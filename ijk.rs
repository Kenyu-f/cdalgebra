// cayley_table.rs
//
// Standalone, dependency-free generator for the Cayley–Dickson multiplication
// table ("累積表") as a colored BMP image, for arbitrary dimension 2^n.
//
// Build:   rustc -O cayley_table.rs -o cayley_table
// Run:     ./cayley_table <n> [cell_px] [out.bmp]
//   n        exponent, table dimension = 2^n (e.g. n=7 -> 128x128 table)
//   cell_px  pixels per cell for upscaling, default 1
//   out.bmp  output filename, default "cayley_table.bmp"
//
// Example: ./cayley_table 8 4 sedenions128.bmp   -> 128x128 table, 4px cells,
//          so a 512x512 image.
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

use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
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

fn color_for_cell(dim: usize, idx: usize, sign: i32) -> (u8, u8, u8) {
    let hue = (idx as f64 / dim as f64) * (300.0 / 360.0); // 0..300deg, matches web version
    let light = if sign < 0 { 0.30 } else { 0.55 };
    hsl_to_rgb(hue, 0.55, light)
}

/// Writes an uncompressed 24-bit BMP. `pixels` is row-major, top row first.
fn write_bmp(path: &str, width: usize, height: usize, pixels: &[(u8, u8, u8)]) -> std::io::Result<()> {
    let row_size = (width * 3 + 3) / 4 * 4; // rows are padded to a multiple of 4 bytes
    let pixel_array_size = row_size * height;
    let file_size = 14 + 40 + pixel_array_size;

    let mut f = BufWriter::new(File::create(path)?);

    // BITMAPFILEHEADER (14 bytes)
    f.write_all(b"BM")?;
    f.write_all(&(file_size as u32).to_le_bytes())?;
    f.write_all(&0u16.to_le_bytes())?; // reserved
    f.write_all(&0u16.to_le_bytes())?; // reserved
    f.write_all(&54u32.to_le_bytes())?; // pixel data offset

    // BITMAPINFOHEADER (40 bytes)
    f.write_all(&40u32.to_le_bytes())?;
    f.write_all(&(width as i32).to_le_bytes())?;
    f.write_all(&(height as i32).to_le_bytes())?; // positive => bottom-up rows
    f.write_all(&1u16.to_le_bytes())?; // color planes
    f.write_all(&24u16.to_le_bytes())?; // bits per pixel
    f.write_all(&0u32.to_le_bytes())?; // no compression
    f.write_all(&(pixel_array_size as u32).to_le_bytes())?;
    f.write_all(&2835i32.to_le_bytes())?; // ~72 DPI
    f.write_all(&2835i32.to_le_bytes())?;
    f.write_all(&0u32.to_le_bytes())?; // colors in palette (0 = default)
    f.write_all(&0u32.to_le_bytes())?; // important colors

    let padding = vec![0u8; row_size - width * 3];
    for y in (0..height).rev() {
        // BMP stores rows bottom-to-top
        for x in 0..width {
            let (r, g, b) = pixels[y * width + x];
            f.write_all(&[b, g, r])?; // BMP wants BGR
        }
        f.write_all(&padding)?;
    }
    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let n: u32 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(6); // default dim=64
    let cell: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4);
    let out = args.get(3).cloned().unwrap_or_else(|| "cayley_table.bmp".to_string());

    if n > 16 {
        eprintln!("n={} would need a {}x{} table — refusing (raise the limit in main() if you really need it).", n, 1usize << n, 1usize << n);
        std::process::exit(1);
    }

    let dim: usize = 1usize << n; // 2^n
    println!("Cayley–Dickson table: n={n}  dim={dim} ({}x{} cells)", dim, dim);

    let start = Instant::now();
    let mut pixels = vec![(0u8, 0u8, 0u8); dim * dim];
    for i in 0..dim {
        for j in 0..dim {
            let (idx, sign) = mul_basis(dim, i, j);
            pixels[i * dim + j] = color_for_cell(dim, idx, sign);
        }
    }
    let elapsed = start.elapsed();
    println!(
        "table computed in {:.3} ms ({} cells, {:.1} ns/cell)",
        elapsed.as_secs_f64() * 1000.0,
        dim * dim,
        elapsed.as_nanos() as f64 / (dim * dim) as f64
    );

    let width = dim * cell;
    let height = dim * cell;
    let mut out_pixels = vec![(0u8, 0u8, 0u8); width * height];
    for y in 0..height {
        for x in 0..width {
            out_pixels[y * width + x] = pixels[(y / cell) * dim + (x / cell)];
        }
    }

    write_bmp(&out, width, height, &out_pixels).expect("failed to write BMP");
    println!("wrote {out} ({width}x{height} px)");
}
