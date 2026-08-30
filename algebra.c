/*
 * cayley_table.c
 *
 * Cayley-Dickson multiplication table generator.
 *
 * Features:
 *   - Recursive multiplication of basis elements
 *   - Colored 24-bit BMP output
 *   - Recursive hierarchy lines
 *   - Radial vignette
 *   - Text multiplication table for small dimensions
 *   - Interactive mode
 *   - Command-line mode
 *   - Streaming BMP output
 *
 * Build:
 *
 *   cc -O3 -std=c11 cayley_table.c -o cayley_table -lm
 *
 * Run:
 *
 *   ./cayley_table
 *   ./cayley_table <n>
 *   ./cayley_table <n> <cell_px>
 *   ./cayley_table <n> <cell_px> <out.bmp>
 *
 * Example:
 *
 *   ./cayley_table 8 4 cayley_256.bmp
 *
 * n = exponent
 * dim = 2^n
 *
 * For n=8:
 *   dim = 256
 *   with cell_px=4:
 *   image = 1024 x 1024 px
 */

/* ============================================================ */
/* Includes                                                     */
/* ============================================================ */

#include <stdio.h>

#include <stdlib.h>

#include <stdint.h>

#include <string.h>

#include <math.h>

#include <time.h>

#include <limits.h>

/* ============================================================ */
/* Types                                                        */
/* ============================================================ */

typedef struct {
  size_t idx;
  int sign;
} MulResult;

/* ============================================================ */
/* Utility                                                      */
/* ============================================================ */

static double clamp01(double x) {
  if (x < 0.0)
    return 0.0;

  if (x > 1.0)
    return 1.0;

  return x;
}

static double mod_positive(double x, double y) {
  double r = fmod(x, y);

  if (r < 0.0)
    r += y;

  return r;
}

/* ============================================================ */
/* Cayley-Dickson multiplication                                 */
/* ============================================================ */

/*
 * Computes
 *
 *     e_i * e_j
 *
 * recursively.
 *
 * The return value is:
 *
 *     (basis index, sign)
 *
 * where sign is +1 or -1.
 *
 * The recursion follows
 *
 *     (a,b)(c,d)
 *       = (ac - conj(d)b,
 *          da + b conj(c))
 *
 * restricted to basis elements.
 */
static MulResult mul_basis(size_t dim, size_t i, size_t j) {
  /* Base case: R */
  if (dim == 1) {
    return (MulResult){.idx = 0, .sign = 1};
  }

  size_t half = dim / 2;

  /*
   * Both elements are in the lower half:
   *
   * (a,0)(c,0) = (ac,0)
   */
  if (i < half && j < half) {

    return mul_basis(half, i, j);
  }

  /*
   * First element is lower,
   * second element is upper:
   *
   * (a,0)(0,d) = (0,da)
   */
  if (i < half && j >= half) {

    size_t jj = j - half;

    MulResult r = mul_basis(half, jj, i);

    r.idx += half;

    return r;
  }

  /*
   * First element is upper,
   * second element is lower:
   *
   * (0,b)(c,0)
   *   = (0,b*conj(c))
   */
  if (i >= half && j < half) {

    size_t ii = i - half;

    /*
     * conj(e_0) = e_0
     * conj(e_k) = -e_k, k > 0
     */
    int conj_sign = (j == 0) ? 1 : -1;

    MulResult r = mul_basis(half, ii, j);

    r.idx += half;
    r.sign *= conj_sign;

    return r;
  }

  /*
   * Both elements are in the upper half:
   *
   * (0,b)(0,d)
   *   = (-conj(d)b,0)
   */
  {
    size_t ii = i - half;

    size_t jj = j - half;

    int conj_sign = (jj == 0) ? 1 : -1;

    MulResult r = mul_basis(half, jj, ii);

    r.sign = -r.sign * conj_sign;

    return r;
  }
}

/* ============================================================ */
/* HSL -> RGB                                                    */
/* ============================================================ */

/*
 * Converts HSL to RGB.
 *
 * h, s, l are in [0,1].
 */
static void hsl_to_rgb(double h, double s, double l, uint8_t *out_r,
                       uint8_t *out_g, uint8_t *out_b) {
  h = mod_positive(h, 1.0);

  double c = (1.0 - fabs(2.0 * l - 1.0)) * s;

  double x = c * (1.0 - fabs(mod_positive(h * 6.0, 2.0) - 1.0));

  double m = l - c / 2.0;

  double r;
  double g;
  double b;

  if (h < 1.0 / 6.0) {

    r = c;
    g = x;
    b = 0.0;

  } else if (h < 2.0 / 6.0) {

    r = x;
    g = c;
    b = 0.0;

  } else if (h < 3.0 / 6.0) {

    r = 0.0;
    g = c;
    b = x;

  } else if (h < 4.0 / 6.0) {

    r = 0.0;
    g = x;
    b = c;

  } else if (h < 5.0 / 6.0) {

    r = x;
    g = 0.0;
    b = c;

  } else {

    r = c;
    g = 0.0;
    b = x;
  }

  *out_r = (uint8_t)lround((r + m) * 255.0);

  *out_g = (uint8_t)lround((g + m) * 255.0);

  *out_b = (uint8_t)lround((b + m) * 255.0);
}

/* ============================================================ */
/* Cell color                                                    */
/* ============================================================ */

static void color_for_cell(size_t dim, size_t idx, int sign, uint8_t *r,
                           uint8_t *g, uint8_t *b) {
  /*
   * Hue:
   *
   * basis index
   * -> color
   */
  double hue = ((double)idx / (double)dim) * (300.0 / 360.0);

  /*
   * Sign:
   *
   * positive -> lighter
   * negative -> darker
   */
  double light = (sign < 0) ? 0.30 : 0.55;

  hsl_to_rgb(hue, 0.55, light, r, g, b);
}

/* ============================================================ */
/* Basis labels                                                   */
/* ============================================================ */

/*
 * Labels:
 *
 * dim <= 4:
 *
 *     1 i j k
 *
 * dim > 4:
 *
 *     1 e1 e2 ...
 */
static void basis_label(size_t dim, size_t idx, char *buffer,
                        size_t buffer_size) {
  if (idx == 0) {

    snprintf(buffer, buffer_size, "1");

    return;
  }

  if (dim <= 4) {

    switch (idx) {

    case 1:

      snprintf(buffer, buffer_size, "i");

      return;

    case 2:

      snprintf(buffer, buffer_size, "j");

      return;

    case 3:

      snprintf(buffer, buffer_size, "k");

      return;

    default:
      break;
    }
  }

  snprintf(buffer, buffer_size, "e%zu", idx);
}

/* ============================================================ */
/* Text multiplication table                                    */
/* ============================================================ */

static void print_text_table(size_t dim) {
  size_t max_width = 2;

  char label[64];

  /*
   * Find largest label width.
   */
  for (size_t k = 0; k < dim; k++) {

    basis_label(dim, k, label, sizeof(label));

    size_t len = strlen(label);

    if (len > max_width)
      max_width = len;
  }

  size_t width = max_width + 1;

  /*
   * Header.
   */
  printf("%*s", (int)(width + 1), "");

  for (size_t j = 0; j < dim; j++) {

    basis_label(dim, j, label, sizeof(label));

    printf("%*s", (int)width, label);
  }

  printf("\n");

  /*
   * Table.
   */
  for (size_t i = 0; i < dim; i++) {

    basis_label(dim, i, label, sizeof(label));

    printf("%*s ", (int)width, label);

    for (size_t j = 0; j < dim; j++) {

      MulResult result = mul_basis(dim, i, j);

      char result_label[64];

      basis_label(dim, result.idx, result_label, sizeof(result_label));

      char cell[80];

      if (result.sign < 0) {

        snprintf(cell, sizeof(cell), "-%s", result_label);

      } else {

        snprintf(cell, sizeof(cell), "%s", result_label);
      }

      printf("%*s", (int)width, cell);
    }

    printf("\n");
  }
}

/* ============================================================ */
/* Pixel manipulation                                            */
/* ============================================================ */

/*
 * Darkens a pixel.
 *
 * factor = 1.0
 *     unchanged
 *
 * factor = 0.0
 *     black
 */
static void darken_pixel(uint8_t *r, uint8_t *g, uint8_t *b, double factor) {
  factor = clamp01(factor);

  *r = (uint8_t)lround((double)(*r) * factor);

  *g = (uint8_t)lround((double)(*g) * factor);

  *b = (uint8_t)lround((double)(*b) * factor);
}

/* ============================================================ */
/* Recursive boundary calculation                               */
/* ============================================================ */

/*
 * Distance from coordinate to the nearest
 * multiple of block_px.
 *
 * Example:
 *
 * block_px = 32
 *
 * coordinates:
 *
 *   0    32    64    96
 *   |-----|-----|-----|
 *
 * return value is distance to one
 * of the boundaries.
 */
static size_t distance_to_boundary(size_t coord, size_t block_px) {
  if (block_px == 0)
    return SIZE_MAX;

  size_t rem = coord % block_px;

  size_t dist_left = rem;

  size_t dist_right = block_px - rem;

  return (dist_left < dist_right) ? dist_left : dist_right;
}

/* ============================================================ */
/* Recursive hierarchy lines                                    */
/* ============================================================ */

/*
 * Draws all recursive Cayley-Dickson
 * subdivision boundaries.
 *
 * For n:
 *
 * level 1:
 *     dim / 2
 *
 * level 2:
 *     dim / 4
 *
 * level 3:
 *     dim / 8
 *
 * ...
 *
 * level n:
 *     1
 */
static void apply_recursive_lines(size_t x, size_t y, size_t dim, size_t cell,
                                  uint32_t n, uint8_t *r, uint8_t *g,
                                  uint8_t *b) {
  if (cell == 0)
    return;

  if (dim == 0)
    return;

  if (n == 0)
    return;

  for (uint32_t level = 1; level <= n; level++) {

    /*
     * Number of cells in one recursive block.
     */
    size_t block_cells = dim >> level;

    if (block_cells == 0)
      break;

    /*
     * Convert from cells to pixels.
     */
    if (block_cells > SIZE_MAX / cell) {
      break;
    }

    size_t block_px = block_cells * cell;

    if (block_px == 0)
      continue;

    /*
     * Width of the line.
     */
    size_t line_width;

    if (block_px >= 64) {

      line_width = 3;

    } else if (block_px >= 16) {

      line_width = 2;

    } else {

      line_width = 1;
    }

    size_t dx = distance_to_boundary(x, block_px);

    size_t dy = distance_to_boundary(y, block_px);

    /*
     * Not close to a boundary.
     */
    if (dx >= line_width && dy >= line_width) {
      continue;
    }

    /*
     * Coarse levels are stronger.
     */
    double strength = 0.18 * pow(0.68, (double)(level - 1));

    double factor = 1.0 - strength;

    darken_pixel(r, g, b, factor);
  }
}

/* ============================================================ */
/* Vignette                                                      */
/* ============================================================ */

/*
 * Applies a radial vignette.
 *
 * Center:
 *     almost unchanged
 *
 * Corners:
 *     darker
 */
static void apply_vignette(size_t x, size_t y, size_t width, size_t height,
                           uint8_t *r, uint8_t *g, uint8_t *b) {
  if (width <= 1 || height <= 1) {
    return;
  }

  double cx = ((double)width - 1.0) / 2.0;

  double cy = ((double)height - 1.0) / 2.0;

  double dx = ((double)x - cx) / cx;

  double dy = ((double)y - cy) / cy;

  /*
   * Radius:
   *
   * 0 = center
   * 1 = corner
   */
  double radius = sqrt(dx * dx + dy * dy) / sqrt(2.0);

  radius = clamp01(radius);

  /*
   * Smooth transition.
   */
  double smooth = radius * radius * (3.0 - 2.0 * radius);

  /*
   * Maximum darkness = 25%.
   */
  double factor = 1.0 - 0.25 * smooth;

  darken_pixel(r, g, b, factor);
}

/* ============================================================ */
/* BMP helpers                                                   */
/* ============================================================ */

static int write_u16_le(FILE *f, uint16_t value) {
  uint8_t bytes[2];

  bytes[0] = (uint8_t)(value & 0xff);

  bytes[1] = (uint8_t)((value >> 8) & 0xff);

  return fwrite(bytes, 1, 2, f) == 2;
}

static int write_u32_le(FILE *f, uint32_t value) {
  uint8_t bytes[4];

  bytes[0] = (uint8_t)(value & 0xff);

  bytes[1] = (uint8_t)((value >> 8) & 0xff);

  bytes[2] = (uint8_t)((value >> 16) & 0xff);

  bytes[3] = (uint8_t)((value >> 24) & 0xff);

  return fwrite(bytes, 1, 4, f) == 4;
}

/* ============================================================ */
/* BMP writer                                                    */
/* ============================================================ */

static int write_bmp_streaming(const char *path, size_t width, size_t height,
                               size_t dim, size_t cell, uint32_t n) {
  /*
   * Standard BMP uses signed 32-bit dimensions.
   */
  if (width > INT32_MAX || height > INT32_MAX) {

    fprintf(stderr, "image dimensions exceed BMP limits\n");

    return 0;
  }

  /*
   * 24-bit BMP:
   *
   * 3 bytes per pixel.
   */
  if (width > (SIZE_MAX - 3) / 3) {

    fprintf(stderr, "row size overflow\n");

    return 0;
  }

  size_t row_size = (width * 3 + 3) / 4 * 4;

  /*
   * Check height multiplication.
   */
  if (height != 0 && row_size > SIZE_MAX / height) {

    fprintf(stderr, "image too large\n");

    return 0;
  }

  size_t pixel_array_size = row_size * height;

  /*
   * Classic BMP file size field
   * is uint32_t.
   */
  if (pixel_array_size > UINT32_MAX - 54) {

    fprintf(stderr, "BMP file exceeds 4 GiB limit\n");

    return 0;
  }

  size_t file_size = 54 + pixel_array_size;

  FILE *f = fopen(path, "wb");

  if (f == NULL) {

    perror("failed to create BMP");

    return 0;
  }

  /* -------------------------------------------------------- */
  /* BITMAPFILEHEADER                                         */
  /* -------------------------------------------------------- */

  if (fwrite("BM", 1, 2, f) != 2 ||

      !write_u32_le(f, (uint32_t)file_size) ||

      !write_u16_le(f, 0) ||

      !write_u16_le(f, 0) ||

      !write_u32_le(f, 54)) {

    fclose(f);

    return 0;
  }

  /* -------------------------------------------------------- */
  /* BITMAPINFOHEADER                                         */
  /* -------------------------------------------------------- */

  /*
   * Height is written positive,
   * therefore BMP rows are bottom-up.
   */
  if (!write_u32_le(f, 40) ||

      !write_u32_le(f, (uint32_t)width) ||

      !write_u32_le(f, (uint32_t)height) ||

      !write_u16_le(f, 1) ||

      !write_u16_le(f, 24) ||

      !write_u32_le(f, 0) ||

      !write_u32_le(f, (uint32_t)pixel_array_size) ||

      !write_u32_le(f, 2835) ||

      !write_u32_le(f, 2835) ||

      !write_u32_le(f, 0) ||

      !write_u32_le(f, 0)) {

    fclose(f);

    return 0;
  }

  /*
   * Only one image row is held in memory.
   */
  uint8_t *row_buf = (uint8_t *)calloc(row_size, 1);

  if (row_buf == NULL) {

    fprintf(stderr, "failed to allocate row buffer\n");

    fclose(f);

    return 0;
  }

  /* -------------------------------------------------------- */
  /* Pixel generation                                         */
  /* -------------------------------------------------------- */

  /*
   * BMP rows are stored bottom -> top.
   */
  for (size_t y = height; y-- > 0;) {

    /*
     * Padding bytes must be zero.
     */
    memset(row_buf, 0, row_size);

    for (size_t x = 0; x < width; x++) {

      /*
       * Map pixel coordinate to table cell.
       */
      size_t i = y / cell;

      size_t j = x / cell;

      /*
       * Compute basis multiplication.
       */
      MulResult result = mul_basis(dim, i, j);

      /*
       * Convert algebraic result
       * to a base color.
       */
      uint8_t r;
      uint8_t g;
      uint8_t b;

      color_for_cell(dim, result.idx, result.sign, &r, &g, &b);

      /*
       * Visual layer 1:
       *
       * recursive hierarchy.
       */
      apply_recursive_lines(x, y, dim, cell, n, &r, &g, &b);

      /*
       * Visual layer 2:
       *
       * radial vignette.
       */
      apply_vignette(x, y, width, height, &r, &g, &b);

      /*
       * BMP uses BGR order.
       */
      size_t offset = x * 3;

      row_buf[offset] = b;

      row_buf[offset + 1] = g;

      row_buf[offset + 2] = r;
    }

    /*
     * Write one row.
     */
    if (fwrite(row_buf, 1, row_size, f) != row_size) {

      fprintf(stderr, "failed to write BMP\n");

      free(row_buf);
      fclose(f);

      return 0;
    }
  }

  free(row_buf);

  if (fclose(f) != 0) {

    perror("failed to close BMP");

    return 0;
  }

  return 1;
}

/* ============================================================ */
/* Human-readable byte size                                     */
/* ============================================================ */

static void human_bytes(double n, char *buffer, size_t buffer_size) {
  const char *units[] = {"B", "KB", "MB", "GB", "TB"};

  size_t unit = 0;

  while (n >= 1024.0 && unit < 4) {

    n /= 1024.0;

    unit++;
  }

  snprintf(buffer, buffer_size, "%.1f%s", n, units[unit]);
}

/* ============================================================ */
/* Input helpers                                                 */
/* ============================================================ */

static void prompt(const char *message, const char *default_value, char *buffer,
                   size_t buffer_size) {
  printf("%s [%s]: ", message, default_value);

  fflush(stdout);

  if (fgets(buffer, (int)buffer_size, stdin) == NULL) {

    snprintf(buffer, buffer_size, "%s", default_value);

    return;
  }

  /*
   * Remove newline.
   */
  buffer[strcspn(buffer, "\r\n")] = '\0';

  /*
   * Empty input -> default.
   */
  if (buffer[0] == '\0') {

    snprintf(buffer, buffer_size, "%s", default_value);
  }
}

/* ============================================================ */
/* Main generation function                                     */
/* ============================================================ */

static int run(uint32_t n, size_t cell, const char *out) {
  /*
   * Avoid shifting into the sign bit
   * or past size_t.
   */
  if (n >= sizeof(size_t) * CHAR_BIT - 1) {

    fprintf(stderr,
            "n=%u is too large: "
            "2^n does not fit in size_t\n",
            n);

    return 1;
  }

  /*
   * dim = 2^n.
   */
  size_t dim = (size_t)1 << n;

  if (cell == 0) {

    fprintf(stderr, "cell size must be positive\n");

    return 1;
  }

  /*
   * width = dim * cell.
   */
  if (dim > SIZE_MAX / cell) {

    fprintf(stderr, "image dimensions overflow\n");

    return 1;
  }

  size_t width = dim * cell;

  size_t height = dim * cell;

  /*
   * Row size for 24-bit BMP.
   */
  if (width > (SIZE_MAX - 3) / 3) {

    fprintf(stderr, "row size overflow\n");

    return 1;
  }

  size_t row_size = (width * 3 + 3) / 4 * 4;

  if (height != 0 && row_size > SIZE_MAX / height) {

    fprintf(stderr, "image size overflow\n");

    return 1;
  }

  double estimated_bytes = (double)row_size * (double)height;

  printf("Cayley-Dickson table: "
         "n=%u  dim=%zu "
         "(%zux%zu cells)\n",
         n, dim, dim, dim);

  /*
   * Small tables are printed as text too.
   */
  if (dim <= 16) {

    printf("\n");

    print_text_table(dim);

    printf("\n");
  }

  char size_text[64];

  human_bytes(estimated_bytes, size_text, sizeof(size_text));

  printf("output image: "
         "%zux%zu px "
         "(~%s on disk)\n",
         width, height, size_text);

  /*
   * Warn about very large outputs.
   */
  if (estimated_bytes > 500000000.0) {

    fprintf(stderr, "note: large file; "
                    "this will use substantial disk space.\n");
  }

  /*
   * CPU time measurement.
   */
  clock_t start = clock();

  int ok = write_bmp_streaming(out, width, height, dim, cell, n);

  clock_t end = clock();

  if (!ok)
    return 1;

  double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

  double megapixels_per_second = ((double)width * (double)height) / 1e6 /
                                 (elapsed > 1e-9 ? elapsed : 1e-9);

  printf("wrote %s in %.3fs "
         "(%.1f Mpixels/s)\n",
         out, elapsed, megapixels_per_second);

  return 0;
}

/* ============================================================ */
/* main                                                          */
/* ============================================================ */

int main(int argc, char **argv) {
  /*
   * No arguments:
   * interactive mode.
   */
  if (argc == 1) {

    char buffer[512];

    printf("Cayley-Dickson "
           "multiplication table generator\n");

    printf("(interactive mode)\n\n");

    /*
     * n
     */
    prompt("n (dimension = 2^n)", "6", buffer, sizeof(buffer));

    char *endptr;

    unsigned long parsed_n = strtoul(buffer, &endptr, 10);

    uint32_t n;

    if (endptr == buffer || *endptr != '\0' || parsed_n > UINT32_MAX) {

      n = 6;

    } else {

      n = (uint32_t)parsed_n;
    }

    /*
     * cell size
     */
    prompt("pixels per cell", "4", buffer, sizeof(buffer));

    unsigned long long parsed_cell = strtoull(buffer, &endptr, 10);

    size_t cell;

    if (endptr == buffer || *endptr != '\0' || parsed_cell == 0 ||
        parsed_cell > SIZE_MAX) {

      cell = 4;

    } else {

      cell = (size_t)parsed_cell;
    }

    /*
     * output filename
     */
    char out[512];

    prompt("output filename", "cayley_table.bmp", out, sizeof(out));

    return run(n, cell, out);
  }

  /*
   * Command-line defaults.
   */
  uint32_t n = 6;

  size_t cell = 4;

  const char *out = "cayley_table.bmp";

  /*
   * Parse n.
   */
  if (argc >= 2) {

    char *endptr;

    unsigned long parsed = strtoul(argv[1], &endptr, 10);

    if (endptr != argv[1] && *endptr == '\0' && parsed <= UINT32_MAX) {

      n = (uint32_t)parsed;
    } else {

      fprintf(stderr, "invalid n: %s\n", argv[1]);

      return 1;
    }
  }

  /*
   * Parse cell size.
   */
  if (argc >= 3) {

    char *endptr;

    unsigned long long parsed = strtoull(argv[2], &endptr, 10);

    if (endptr != argv[2] && *endptr == '\0' && parsed > 0 &&
        parsed <= SIZE_MAX) {

      cell = (size_t)parsed;
    } else {

      fprintf(stderr, "invalid cell size: %s\n", argv[2]);

      return 1;
    }
  }

  /*
   * Output filename.
   */
  if (argc >= 4) {

    out = argv[3];
  }

  /*
   * Ignore additional arguments,
   * matching the simple command-line
   * interface of the original program.
   */
  return run(n, cell, out);
}
