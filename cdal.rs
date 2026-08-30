// cd.rs  (formerly cayley_table.rs)
//
// Standalone, dependency-free generator for the Cayley–Dickson multiplication
// table ("累積表") as a colored BMP image, for arbitrary dimension 2^n —
// now with a fully customizable color gradient.
//
// Build:   rustc -O cd.rs -o cd
// Run:     ./cd <n> [cell_px] [out.bmp] [hue_start] [hue_end] [saturation] [light_pos] [light_neg]
//   n           exponent, table dimension = 2^n (e.g. n=7 -> 128x128 table).
//               Any n you like — no artificial cap.
//   cell_px     pixels per cell for upscaling, default 4
//   out.bmp     output filename, default "cayley_table.bmp"
//   hue_start   gradient start, degrees 0-360, default 0
//   hue_end     gradient end,   degrees 0-360, default 300
//   saturation  percent 0-100, default 55
//   light_pos   lightness % for positive-sign cells, default 55
//   light_neg   lightness % for negative-sign cells, default 30
//
// Run with NO arguments at all to get a fully interactive prompt instead —
// it will ask for n, cell size, output filename, and then let you pick a
// color gradient from a small preset menu (rainbow / ocean / fire / mono)
// or enter your own hue/saturation/lightness values, one question at a
// time, the same way Python's input() would.
//
// Example: ./cd 8 4 sedenions128.bmp   -> 256x256 table, 4px cells,
//          so a 1024x1024 image, rainbow gradient (defaults).
//
// -----------------------------------------------------------------------
// Why this is fast even for large dimensions
// -----------------------------------------------------------------------
// The browser version multiplied full 2^n-component vectors recursively
// (the textbook Cayley–Dickson formula (a,b)(c,d) = (ac - d̄b, da + bc̄)
// applied to dense arrays). That costs O(dim^2) per single multiplication,
// so filling the whole dim x dim table costs O(dim^4) — fine up to dim=32,
// unusably slow beyond that (this is exactly the "32x32 is already too big"
// wall you hit).
//
// The key fact that fixes this: the product of two *basis* elements e_i, e_j
// in any Cayley–Dickson algebra is never a combination of several basis
// vectors — it is always exactly ±e_k for some single k. So instead of
// multiplying whole vectors, we only ever need to track a single
// (index, sign) pair and recurse on the SAME formula restricted to that
// pair. That reduces one cell of the table to O(log dim) work, and the
// whole table to O(dim^2 log dim) — instant even for dim in the thousands.
//
// This was derived by case-splitting (a,b)(c,d) = (ac-d̄b, da+bc̄) on which
// half (a-part / b-part) each basis vector occupies, and verified
// exhaustively (all i,j pairs) against the full vector-based multiplication
// for dim = 2,4,8,16,32 before being used here.
//
// The image itself is now streamed row-by-row straight to disk (see
// write_bmp_streaming) rather than being built up in one giant in-memory
// buffer first — so memory use stays proportional to image *width*, not to
// the full width*height, which is what lets n grow arbitrarily large without
// running out of RAM (disk space and time are the only real limits left).

use std::env;
use std::fs::File;
use std::io::{self, BufWriter, Write, BufRead};
use std::time::Instant;

/// Recursively computes e_i * e_j in the 2^k-dimensional Cayley–Dickson
/// algebra ("dim" = 2^k, passed explicitly so the recursion can shrink it).
/// Returns (result_basis_index, sign).
fn mul_basis(dim: usize, i: usize, j: usize) -> (usize, i32) {
    if dim == 1 {
        return (0, 1); // the base case: 1 * 1 = 1
    }
    let half = dim / 2;
    match (i < half, j < half) {
        // both in the "a" (lower) half: (a,0)(c,0) = (ac, 0) — stays put
        (true, true) => mul_basis(half, i, j),

        // x in "a", y in "b": (a,0)(0,d) = (0, da) — lands in upper half
        (true, false) => {
            let jj = j - half;
            let (idx, s) = mul_basis(half, jj, i); // order: d * a
            (half + idx, s)
        }

        // x in "b", y in "a": (0,b)(c,0) = (0, b*conj(c)) — upper half
        (false, true) => {
            let ii = i - half;
            let conj_sign: i32 = if j == 0 { 1 } else { -1 }; // conj(e_0)=e_0, conj(e_k)=-e_k
            let (idx, s) = mul_basis(half, ii, j);
            (half + idx, s * conj_sign)
        }

        // both in "b": (0,b)(0,d) = (-conj(d)*b, 0) — lower half, extra minus
        (false, false) => {
            let ii = i - half;
            let jj = j - half;
            let conj_sign: i32 = if jj == 0 { 1 } else { -1 };
            let (idx, s) = mul_basis(half, jj, ii); // order: conj(d) * b
            (idx, -s * conj_sign)
        }
    }
}

/// HSL -> RGB with h, s, l in [0,1]. Same scheme as the browser artifact:
/// hue encodes which basis index the product landed on, lightness encodes
/// the sign (dark = negative, light = positive).
fn hsl_to_rgb(h: f64, s: f64, l: f64) -> (u8, u8, u8) {
    let h = h.rem_euclid(1.0);
    let c = (1.0 - (2.0 * l - 1.0).abs()) * s;
    let x = c * (1.0 - ((h * 6.0).rem_euclid(2.0) - 1.0).abs());
    let m = l - c / 2.0;
    let (r, g, b) = if h < 1.0 / 6.0 {
        (c, x, 0.0)
    } else if h < 2.0 / 6.0 {
        (x, c, 0.0)
    } else if h < 3.0 / 6.0 {
        (0.0, c, x)
    } else if h < 4.0 / 6.0 {
        (0.0, x, c)
    } else if h < 5.0 / 6.0 {
        (x, 0.0, c)
    } else {
        (c, 0.0, x)
    };
    (
        ((r + m) * 255.0).round() as u8,
        ((g + m) * 255.0).round() as u8,
        ((b + m) * 255.0).round() as u8,
    )
}

/// A customizable color gradient for the table: which hue range basis
/// indices sweep across, how saturated the colors are, and how light
/// positive- vs negative-sign cells are. All angles in degrees (0-360),
/// saturation/lightness as percentages (0-100) — the human-friendly units
/// used at the prompts; converted to HSL's 0..1 range in `color_for_cell`.
#[derive(Clone, Copy)]
struct ColorScheme {
    hue_start: f64,
    hue_end: f64,
    saturation: f64,
    light_pos: f64,
    light_neg: f64,
}
impl ColorScheme {
    fn rainbow() -> Self { Self { hue_start: 0.0,   hue_end: 300.0, saturation: 55.0, light_pos: 55.0, light_neg: 30.0 } }
    fn ocean()   -> Self { Self { hue_start: 180.0, hue_end: 260.0, saturation: 65.0, light_pos: 60.0, light_neg: 25.0 } }
    fn fire()    -> Self { Self { hue_start: 0.0,   hue_end: 50.0,  saturation: 80.0, light_pos: 60.0, light_neg: 25.0 } }
    fn mono()    -> Self { Self { hue_start: 0.0,   hue_end: 0.0,   saturation: 0.0,  light_pos: 72.0, light_neg: 22.0 } }
}

fn color_for_cell(dim: usize, idx: usize, sign: i32, scheme: &ColorScheme) -> (u8, u8, u8) {
    let t = idx as f64 / dim as f64;
    let hue_deg = scheme.hue_start + (scheme.hue_end - scheme.hue_start) * t;
    let light_pct = if sign < 0 { scheme.light_neg } else { scheme.light_pos };
    hsl_to_rgb(hue_deg / 360.0, scheme.saturation / 100.0, light_pct / 100.0)
}

/// Human-readable basis label: 1, i, j, k for quaternions (dim<=4);
/// e1..e_{dim-1} beyond that (octonions, sedenions, ...).
fn basis_label(dim: usize, idx: usize) -> String {
    if idx == 0 {
        return "1".to_string();
    }
    if dim <= 4 {
        match idx {
            1 => "i".to_string(),
            2 => "j".to_string(),
            3 => "k".to_string(),
            _ => format!("e{idx}"),
        }
    } else {
        format!("e{idx}")
    }
}

/// Prints a labeled i*j multiplication table to the terminal — readable for
/// small algebras like quaternions (dim=4) or octonions (dim=8), where a
/// 4x4 or 8x8 color image is too small to read by eye anyway.
fn print_text_table(dim: usize) {
    let labels: Vec<String> = (0..dim).map(|k| basis_label(dim, k)).collect();
    let width = labels.iter().map(|s| s.len()).max().unwrap_or(1).max(2) + 1;

    print!("{:width$}", "", width = width + 1);
    for l in &labels {
        print!("{:>width$}", l, width = width);
    }
    println!();

    for i in 0..dim {
        print!("{:>width$} ", labels[i], width = width);
        for j in 0..dim {
            let (idx, sign) = mul_basis(dim, i, j);
            let cell = if sign < 0 {
                format!("-{}", labels[idx])
            } else {
                labels[idx].clone()
            };
            print!("{:>width$}", cell, width = width);
        }
        println!();
    }
}

/// Writes an uncompressed 24-bit BMP, computing each pixel on the fly via
/// `pixel_fn(x,y)` and writing one row at a time — the whole width*height
/// image is never held in memory at once, only a single row buffer.
fn write_bmp_streaming<F>(path: &str, width: usize, height: usize, mut pixel_fn: F) -> io::Result<()>
where
    F: FnMut(usize, usize) -> (u8, u8, u8),
{
    let row_size = (width * 3 + 3) / 4 * 4; // rows padded to a multiple of 4 bytes
    let pixel_array_size = row_size
        .checked_mul(height)
        .expect("image too large: row_size * height overflows");
    let file_size = 14usize + 40 + pixel_array_size;

    let mut f = BufWriter::with_capacity(1 << 20, File::create(path)?);

    // BITMAPFILEHEADER (14 bytes)
    f.write_all(b"BM")?;
    f.write_all(&(file_size as u32).to_le_bytes())?;
    f.write_all(&0u16.to_le_bytes())?;
    f.write_all(&0u16.to_le_bytes())?;
    f.write_all(&54u32.to_le_bytes())?;

    // BITMAPINFOHEADER (40 bytes)
    f.write_all(&40u32.to_le_bytes())?;
    f.write_all(&(width as i32).to_le_bytes())?;
    f.write_all(&(height as i32).to_le_bytes())?; // positive => bottom-up rows
    f.write_all(&1u16.to_le_bytes())?;
    f.write_all(&24u16.to_le_bytes())?;
    f.write_all(&0u32.to_le_bytes())?;
    f.write_all(&(pixel_array_size as u32).to_le_bytes())?;
    f.write_all(&2835i32.to_le_bytes())?;
    f.write_all(&2835i32.to_le_bytes())?;
    f.write_all(&0u32.to_le_bytes())?;
    f.write_all(&0u32.to_le_bytes())?;

    let mut row_buf = vec![0u8; row_size]; // padding bytes stay 0 automatically
    for y in (0..height).rev() {
        // BMP stores rows bottom-to-top
        for x in 0..width {
            let (r, g, b) = pixel_fn(x, y);
            let o = x * 3;
            row_buf[o] = b;
            row_buf[o + 1] = g;
            row_buf[o + 2] = r;
        }
        f.write_all(&row_buf)?;
    }
    Ok(())
}

fn human_bytes(n: f64) -> String {
    let units = ["B", "KB", "MB", "GB", "TB"];
    let mut v = n;
    let mut u = 0;
    while v >= 1024.0 && u < units.len() - 1 {
        v /= 1024.0;
        u += 1;
    }
    format!("{:.1}{}", v, units[u])
}

fn prompt(msg: &str, default: &str) -> String {
    print!("{msg} [{default}]: ");
    io::stdout().flush().ok();
    let mut line = String::new();
    io::stdin().lock().read_line(&mut line).ok();
    let s = line.trim();
    if s.is_empty() { default.to_string() } else { s.to_string() }
}

fn prompt_f64(msg: &str, default: f64) -> f64 {
    prompt(msg, &default.to_string()).parse().unwrap_or(default)
}

/// Interactive color-scheme picker: a short menu (presets 1-4) plus a
/// "custom" option that asks for hue range / saturation / lightness
/// directly, the same way the rest of the interactive mode works.
fn prompt_color_scheme() -> ColorScheme {
    println!("色のグラデーションを選んでください:");
    println!("  1) 虹色   (色相 0-300度、デフォルト)");
    println!("  2) 海     (青〜シアン系)");
    println!("  3) 炎     (赤〜黄色系)");
    println!("  4) モノクロ (彩度なし、明暗のみ)");
    println!("  5) カスタム (色相・彩度・明度を自分で指定)");
    match prompt("選択", "1").as_str() {
        "2" => ColorScheme::ocean(),
        "3" => ColorScheme::fire(),
        "4" => ColorScheme::mono(),
        "5" => {
            let hue_start = prompt_f64("  色相の開始 (0-360度)", 0.0);
            let hue_end = prompt_f64("  色相の終了 (0-360度)", 300.0);
            let saturation = prompt_f64("  彩度 (0-100%)", 55.0);
            let light_pos = prompt_f64("  明度：符号が正のセル (0-100%)", 55.0);
            let light_neg = prompt_f64("  明度：符号が負のセル (0-100%)", 30.0);
            ColorScheme { hue_start, hue_end, saturation, light_pos, light_neg }
        }
        _ => ColorScheme::rainbow(),
    }
}

fn run(n: u32, cell: usize, out: String, scheme: ColorScheme) {
    // 1usize << n overflows for n >= 64; guard it explicitly instead of panicking.
    let dim: usize = match 1usize.checked_shl(n) {
        Some(d) if n < 63 => d,
        _ => {
            eprintln!("n={n} is too large — 2^{n} does not fit in a 64-bit table dimension.");
            std::process::exit(1);
        }
    };

    let width = dim.checked_mul(cell).expect("width overflow: dim*cell too large");
    let height = dim.checked_mul(cell).expect("height overflow: dim*cell too large");
    let row_size = (width * 3 + 3) / 4 * 4;
    let est_bytes = row_size as f64 * height as f64;

    println!("Cayley–Dickson table: n={n}  dim={dim} ({dim}x{dim} cells)");
    println!(
        "color: hue {:.0}°→{:.0}°  sat {:.0}%  light +{:.0}%/-{:.0}%",
        scheme.hue_start, scheme.hue_end, scheme.saturation, scheme.light_pos, scheme.light_neg
    );

    if dim <= 16 {
        println!();
        print_text_table(dim);
        println!();
    }

    println!("output image: {width}x{height} px  (~{} on disk)", human_bytes(est_bytes));
    if est_bytes > 500_000_000.0 {
        eprintln!(
            "note: that's a large file — this will take a while and use real disk space. \
             Ctrl+C now if that's not what you wanted; proceeding in 3s..."
        );
        std::thread::sleep(std::time::Duration::from_secs(3));
    }

    let start = Instant::now();
    write_bmp_streaming(&out, width, height, |x, y| {
        let i = y / cell;
        let j = x / cell;
        let (idx, sign) = mul_basis(dim, i, j);
        color_for_cell(dim, idx, sign, &scheme)
    })
    .expect("failed to write BMP");
    let elapsed = start.elapsed();

    println!(
        "wrote {out} in {:.3}s ({:.1} Mpixels/s)",
        elapsed.as_secs_f64(),
        (width as f64 * height as f64) / 1e6 / elapsed.as_secs_f64().max(1e-9)
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() == 1 {
        // No arguments at all: friendly interactive mode.
        println!("Cayley–Dickson 乗積表ジェネレータ（対話モード。Enterでデフォルト値）");
        let n: u32 = prompt("n（次元 = 2^n）", "6").parse().unwrap_or(6);
        let cell: usize = prompt("1マスのピクセル数", "4").parse().unwrap_or(4);
        let out = prompt("出力ファイル名", "cayley_table.bmp");
        let scheme = prompt_color_scheme();
        run(n, cell, out, scheme);
        return;
    }

    let n: u32 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(6);
    let cell: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4);
    let out = args.get(3).cloned().unwrap_or_else(|| "cayley_table.bmp".to_string());

    // Optional CLI color overrides, in order, after the output filename:
    //   <n> <cell> <out> [hue_start] [hue_end] [saturation] [light_pos] [light_neg]
    // Any that are omitted fall back to the rainbow preset's values.
    let base = ColorScheme::rainbow();
    let get = |idx: usize, default: f64| args.get(idx).and_then(|s| s.parse().ok()).unwrap_or(default);
    let scheme = ColorScheme {
        hue_start: get(4, base.hue_start),
        hue_end: get(5, base.hue_end),
        saturation: get(6, base.saturation),
        light_pos: get(7, base.light_pos),
        light_neg: get(8, base.light_neg),
    };

    run(n, cell, out, scheme);
}
