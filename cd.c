/*
 * cayley_table.c
 *
 * Standalone, dependency-free generator for the Cayley–Dickson
 * multiplication table as a colored BMP image.
 *
 * Build:
 *   cc -O3 -std=c11 cayley_table.c -o cayley_table
 *
 * Run:
 *   ./cayley_table <n> [cell_px] [out.bmp]
 *
 * Example:
 *   ./cayley_table 8 4 sedenions128.bmp
 */

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
  size_t idx;
  int sign;
} MulResult;

/* ------------------------------------------------------------ */
/* Cayley–Dickson basis multiplication                          */
/* ------------------------------------------------------------ */

/*
 * Recursively computes e_i * e_j in the dim-dimensional
 * Cayley–Dickson algebra.
 *
 * Returns:
 *   idx  = resulting basis index
 *   sign = +1 or -1
 */
static MulResult mul_basis(size_t dim, size_t i, size_t j) {
  if (dim == 1) {
    return (MulResult){0, 1};
  }

  size_t half = dim / 2;

  if (i < half && j < half) {
    /* (a,0)(c,0) = (ac,0) */
    return mul_basis(half, i, j);
  }

  if (i < half && j >= half) {
    /* (a,0)(0,d) = (0,da) */
    size_t jj = j - half;
    MulResult r = mul_basis(half, jj, i);

    r.idx += half;
    return r;
  }

  if (i >= half && j < half) {
    /* (0,b)(c,0) = (0,b*conj(c)) */

    size_t ii = i - half;
    int conj_sign = (j == 0) ? 1 : -1;

    MulResult r = mul_basis(half, ii, j);

    r.idx += half;
    r.sign *= conj_sign;

    return r;
  }

  /*
   * (0,b)(0,d) = (-conj(d)*b,0)
   */

  size_t ii = i - half;
  size_t jj = j - half;

  int conj_sign = (jj == 0) ? 1 : -1;

  MulResult r = mul_basis(half, jj, ii);

  r.sign = -r.sign * conj_sign;

  return r;
}

/* ------------------------------------------------------------ */
/* Color conversion                                              */
/* ------------------------------------------------------------ */

static double mod_positive(double x, double y) {
  double r = fmod(x, y);

  if (r < 0.0) {
    r += y;
  }

  return r;
}

/*
 * HSL -> RGB
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

static void color_for_cell(size_t dim, size_t idx, int sign, uint8_t *r,
                           uint8_t *g, uint8_t *b) {
  double hue = ((double)idx / (double)dim) * (300.0 / 360.0);

  double light = (sign < 0) ? 0.30 : 0.55;

  hsl_to_rgb(hue, 0.55, light, r, g, b);
}

/* ------------------------------------------------------------ */
/* Basis labels                                                  */
/* ------------------------------------------------------------ */

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

/* ------------------------------------------------------------ */
/* Text table                                                    */
/* ------------------------------------------------------------ */

static void print_text_table(size_t dim) {
  size_t max_width = 2;

  char label[64];

  for (size_t k = 0; k < dim; k++) {

    basis_label(dim, k, label, sizeof(label));

    size_t len = strlen(label);

    if (len > max_width) {
      max_width = len;
    }
  }

  size_t width = max_width + 1;

  printf("%*s", (int)(width + 1), "");

  for (size_t j = 0; j < dim; j++) {

    basis_label(dim, j, label, sizeof(label));

    printf("%*s", (int)width, label);
  }

  printf("\n");

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

/* ------------------------------------------------------------ */
/* Little-endian BMP helpers                                    */
/* ------------------------------------------------------------ */

static int write_u16_le(FILE *f, uint16_t value) {
  uint8_t b[2];

  b[0] = (uint8_t)(value & 0xff);
  b[1] = (uint8_t)((value >> 8) & 0xff);

  return fwrite(b, 1, 2, f) == 2;
}

static int write_u32_le(FILE *f, uint32_t value) {
  uint8_t b[4];

  b[0] = (uint8_t)(value & 0xff);

  b[1] = (uint8_t)((value >> 8) & 0xff);

  b[2] = (uint8_t)((value >> 16) & 0xff);

  b[3] = (uint8_t)((value >> 24) & 0xff);

  return fwrite(b, 1, 4, f) == 4;
}

/* ------------------------------------------------------------ */
/* BMP writer                                                    */
/* ------------------------------------------------------------ */

static int write_bmp_streaming(const char *path, size_t width, size_t height,
                               size_t dim, size_t cell) {
  if (width > INT32_MAX || height > INT32_MAX) {

    fprintf(stderr, "image dimensions exceed BMP limits\n");

    return 0;
  }

  size_t row_size = ((width * 3 + 3) / 4) * 4;

  if (height != 0 && row_size > SIZE_MAX / height) {

    fprintf(stderr, "image too large\n");

    return 0;
  }

  size_t pixel_array_size = row_size * height;

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
  /* BITMAPFILEHEADER                                        */
  /* -------------------------------------------------------- */

  fwrite("BM", 1, 2, f);

  if (!write_u32_le(f, (uint32_t)file_size) ||

      !write_u16_le(f, 0) || !write_u16_le(f, 0) ||

      !write_u32_le(f, 54)) {

    fclose(f);
    return 0;
  }

  /* -------------------------------------------------------- */
  /* BITMAPINFOHEADER                                        */
  /* -------------------------------------------------------- */

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

  uint8_t *row_buf = calloc(row_size, 1);

  if (row_buf == NULL) {

    fprintf(stderr, "failed to allocate row buffer\n");

    fclose(f);

    return 0;
  }

  /*
   * BMP stores positive-height images
   * from bottom to top.
   */

  for (size_t y = height; y-- > 0;) {

    memset(row_buf, 0, row_size);

    for (size_t x = 0; x < width; x++) {

      size_t i = y / cell;

      size_t j = x / cell;

      MulResult result = mul_basis(dim, i, j);

      uint8_t r;
      uint8_t g;
      uint8_t b;

      color_for_cell(dim, result.idx, result.sign, &r, &g, &b);

      size_t offset = x * 3;

      /*
       * BMP pixel order:
       * B, G, R
       */

      row_buf[offset] = b;

      row_buf[offset + 1] = g;

      row_buf[offset + 2] = r;
    }

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

/* ------------------------------------------------------------ */
/* Human-readable byte size                                     */
/* ------------------------------------------------------------ */

static void human_bytes(double n, char *buffer, size_t buffer_size) {
  const char *units[] = {"B", "KB", "MB", "GB", "TB"};

  size_t unit = 0;

  while (n >= 1024.0 && unit < 4) {

    n /= 1024.0;
    unit++;
  }

  snprintf(buffer, buffer_size, "%.1f%s", n, units[unit]);
}

/* ------------------------------------------------------------ */
/* Interactive prompt                                            */
/* ------------------------------------------------------------ */

static void prompt(const char *message, const char *default_value, char *buffer,
                   size_t buffer_size) {
  printf("%s [%s]: ", message, default_value);

  fflush(stdout);

  if (fgets(buffer, (int)buffer_size, stdin) == NULL) {

    snprintf(buffer, buffer_size, "%s", default_value);

    return;
  }

  buffer[strcspn(buffer, "\r\n")] = '\0';

  if (buffer[0] == '\0') {

    snprintf(buffer, buffer_size, "%s", default_value);
  }
}

/* ------------------------------------------------------------ */
/* Main execution                                                */
/* ------------------------------------------------------------ */

static int run(uint32_t n, size_t cell, const char *out) {
  /*
   * dim = 2^n
   */

  if (n >= sizeof(size_t) * 8 - 1) {

    fprintf(stderr,
            "n=%u is too large: "
            "2^n does not fit in size_t\n",
            n);

    return 1;
  }

  size_t dim = (size_t)1 << n;

  if (cell == 0) {

    fprintf(stderr, "cell size must be positive\n");

    return 1;
  }

  if (dim > SIZE_MAX / cell) {

    fprintf(stderr, "image dimensions overflow\n");

    return 1;
  }

  size_t width = dim * cell;

  size_t height = dim * cell;

  size_t row_size = ((width * 3 + 3) / 4) * 4;

  double estimated_bytes = (double)row_size * (double)height;

  printf("Cayley-Dickson table: "
         "n=%u  dim=%zu "
         "(%zux%zu cells)\n",
         n, dim, dim, dim);

  if (dim <= 16) {

    printf("\n");

    print_text_table(dim);

    printf("\n");
  }

  char size_text[64];

  human_bytes(estimated_bytes, size_text, sizeof(size_text));

  printf("output image: "
         "%zux%zu px  "
         "(~%s on disk)\n",
         width, height, size_text);

  if (estimated_bytes > 500000000.0) {

    fprintf(stderr, "note: large file; "
                    "using real disk space. "
                    "Proceeding immediately...\n");
  }

  clock_t start = clock();

  int ok = write_bmp_streaming(out, width, height, dim, cell);

  clock_t end = clock();

  if (!ok) {
    return 1;
  }

  double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

  double megapixels_per_second = ((double)width * (double)height) / 1e6 /
                                 (elapsed > 1e-9 ? elapsed : 1e-9);

  printf("wrote %s in %.3fs "
         "(%.1f Mpixels/s)\n",
         out, elapsed, megapixels_per_second);

  return 0;
}

/* ------------------------------------------------------------ */
/* main                                                          */
/* ------------------------------------------------------------ */

int main(int argc, char **argv) {
  if (argc == 1) {

    char buffer[256];

    printf("Cayley-Dickson "
           "multiplication table generator "
           "(interactive mode)\n");

    prompt("n (dimension = 2^n)", "6", buffer, sizeof(buffer));

    uint32_t n = (uint32_t)strtoul(buffer, NULL, 10);

    prompt("pixels per cell", "4", buffer, sizeof(buffer));

    size_t cell = (size_t)strtoull(buffer, NULL, 10);

    char out[512];

    prompt("output filename", "cayley_table.bmp", out, sizeof(out));

    return run(n, cell, out);
  }

  uint32_t n = 6;

  size_t cell = 4;

  const char *out = "cayley_table.bmp";

  if (argc >= 2) {

    n = (uint32_t)strtoul(argv[1], NULL, 10);
  }

  if (argc >= 3) {

    cell = (size_t)strtoull(argv[2], NULL, 10);
  }

  if (argc >= 4) {

    out = argv[3];
  }

  return run(n, cell, out);
}

/* ------------------------------------------------------------ */
/* Visual effects                                                */
/* ------------------------------------------------------------ */

static double clamp01(double x) {
  if (x < 0.0)
    return 0.0;

  if (x > 1.0)
    return 1.0;

  return x;
}

/*
 * Darkens a pixel.
 */
static void darken_pixel(uint8_t *r, uint8_t *g, uint8_t *b, double factor) {
  factor = clamp01(factor);

  *r = (uint8_t)lround((double)(*r) * factor);
  *g = (uint8_t)lround((double)(*g) * factor);
  *b = (uint8_t)lround((double)(*b) * factor);
}

/*
 * Returns the distance from a coordinate to the nearest
 * recursive boundary at the specified scale.
 *
 * block_px = size of one recursive block in pixels.
 */
static size_t distance_to_boundary(size_t coord, size_t block_px) {
  if (block_px == 0)
    return SIZE_MAX;

  size_t rem = coord % block_px;

  size_t dist_left = rem;
  size_t dist_right = block_px - rem;

  return dist_left < dist_right ? dist_left : dist_right;
}

/*
 * Draws the hierarchy of the Cayley–Dickson recursive
 * decomposition over the generated image.
 *
 * Coarser recursive divisions receive stronger lines.
 * Finer divisions receive weaker lines.
 */
static void apply_recursive_lines(size_t x, size_t y, size_t dim, size_t cell,
                                  uint32_t n, uint8_t *r, uint8_t *g,
                                  uint8_t *b) {
  if (cell == 0 || dim == 0 || n == 0)
    return;

  /*
   * level = 0:
   *   dim / 2
   *
   * level = 1:
   *   dim / 4
   *
   * ...
   *
   * level = n - 1:
   *   1 cell
   */
  for (uint32_t level = 1; level <= n; level++) {

    size_t block_cells = dim >> level;

    if (block_cells == 0)
      break;

    size_t block_px = block_cells * cell;

    if (block_px == 0)
      continue;

    /*
     * The line becomes thinner at deeper levels.
     *
     * Coarse:
     *   2–3 px
     *
     * Fine:
     *   1 px
     */
    size_t line_width = block_px >= 64 ? 3 : block_px >= 16 ? 2 : 1;

    size_t dx = distance_to_boundary(x, block_px);

    size_t dy = distance_to_boundary(y, block_px);

    if (dx >= line_width && dy >= line_width) {

      continue;
    }

    /*
     * Coarser levels are more visible.
     *
     * level 1 ≈ 0.18
     * level 2 ≈ 0.12
     * level 3 ≈ 0.08
     * ...
     */
    double strength = 0.18 * pow(0.68, (double)(level - 1));

    /*
     * Darkening factor.
     */
    double factor = 1.0 - strength;

    darken_pixel(r, g, b, factor);
  }
}

/*
 * Applies a radial vignette.
 *
 * Center:
 *   factor = 1.0
 *
 * Corners:
 *   factor ≈ 0.75
 */
static void apply_vignette(size_t x, size_t y, size_t width, size_t height,
                           uint8_t *r, uint8_t *g, uint8_t *b) {
  if (width <= 1 || height <= 1)
    return;

  double cx = ((double)width - 1.0) / 2.0;

  double cy = ((double)height - 1.0) / 2.0;

  double dx = ((double)x - cx) / cx;

  double dy = ((double)y - cy) / cy;

  /*
   * Normalized radial distance.
   *
   * 0 at center
   * 1 at corners
   */
  double radius = sqrt(dx * dx + dy * dy) / sqrt(2.0);

  radius = clamp01(radius);

  /*
   * Smooth radial attenuation.
   *
   * center -> 1.00
   * edge   -> 0.84
   * corner -> 0.75
   */
  double edge = radius * radius;

  double factor = 1.0 - 0.25 * edge;

  /*
   * Additional smoothstep-like shaping
   * prevents a hard transition around the center.
   */
  double smooth = radius * radius * (3.0 - 2.0 * radius);

  factor = 1.0 - 0.25 * smooth;

  darken_pixel(r, g, b, factor);
}
