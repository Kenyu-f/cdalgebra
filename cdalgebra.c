/*
 * cayley_table.c
 *
 * Cayley-Dickson multiplication table visualizer
 *
 * Features
 *   - Cayley-Dickson basis multiplication
 *   - 24-bit BMP output
 *   - Interactive visualization settings
 *   - Multiple color palettes
 *   - Mathematical quantity -> visual attribute mapping
 *   - Recursive hierarchy lines
 *   - Vignette
 *   - Interactive regeneration
 *
 * Build:
 *
 *   cc -O3 -std=c11 cayley_table.c -o cayley_table -lm
 *
 * Run:
 *
 *   ./cayley_table
 *
 * or:
 *
 *   ./cayley_table <n> [cell_px] [output.bmp]
 *
 * Example:
 *
 *   ./cayley_table 8 4 cayley.bmp
 */

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ============================================================ */
/* Constants                                                     */
/* ============================================================ */

#define PI 3.14159265358979323846

/* ============================================================ */
/* Basic algebraic result                                        */
/* ============================================================ */

typedef struct {
  size_t idx;
  int sign;
} MulResult;

/* ============================================================ */
/* Visualization configuration                                   */
/* ============================================================ */

typedef enum {
  PALETTE_RAINBOW = 1,
  PALETTE_OCEAN,
  PALETTE_FIRE,
  PALETTE_PLASMA,
  PALETTE_MONOCHROME,
  PALETTE_DIVERGING,
  PALETTE_CUSTOM
} Palette;

typedef enum {
  INDEX_TO_HUE = 1,
  INDEX_TO_SATURATION,
  INDEX_TO_LIGHTNESS
} IndexMapping;

typedef enum { SIGN_LIGHTNESS = 1, SIGN_HUE, SIGN_SATURATION } SignMapping;

typedef struct {

  Palette palette;

  IndexMapping index_mapping;

  SignMapping sign_mapping;

  double hue_shift;

  double saturation;

  double lightness;

  double contrast;

  int recursive_lines;

  double recursive_strength;

  int vignette;

  double vignette_strength;

  /*
   * Custom palette parameters.
   */
  double custom_hue_positive;
  double custom_hue_negative;

  double custom_saturation;

  double custom_lightness_positive;
  double custom_lightness_negative;

} VisualConfig;

/* ============================================================ */
/* Utility                                                       */
/* ============================================================ */

static double clamp01(double x) {
  if (x < 0.0)
    return 0.0;

  if (x > 1.0)
    return 1.0;

  return x;
}

static double lerp(double a, double b, double t) { return a + (b - a) * t; }

static double mod_positive(double x, double y) {
  double r = fmod(x, y);

  if (r < 0.0)
    r += y;

  return r;
}

static int read_int(const char *message, int default_value) {
  char buffer[128];

  printf("%s [%d]: ", message, default_value);

  fflush(stdout);

  if (fgets(buffer, sizeof(buffer), stdin) == NULL) {
    return default_value;
  }

  char *endptr;

  long value = strtol(buffer, &endptr, 10);

  if (endptr == buffer) {
    return default_value;
  }

  return (int)value;
}

static double read_double(const char *message, double default_value) {
  char buffer[128];

  printf("%s [%.3f]: ", message, default_value);

  fflush(stdout);

  if (fgets(buffer, sizeof(buffer), stdin) == NULL) {
    return default_value;
  }

  char *endptr;

  double value = strtod(buffer, &endptr);

  if (endptr == buffer) {
    return default_value;
  }

  return value;
}

static void read_string(const char *message, const char *default_value,
                        char *buffer, size_t buffer_size) {
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

/* ============================================================ */
/* Cayley-Dickson multiplication                                 */
/* ============================================================ */

/*
 * Computes
 *
 *     e_i * e_j
 *
 * and returns
 *
 *     +/- e_k
 *
 * as (idx, sign).
 */
static MulResult mul_basis(size_t dim, size_t i, size_t j) {
  if (dim == 1) {
    return (MulResult){0, 1};
  }

  size_t half = dim / 2;

  /*
   * (a,0)(c,0) = (ac,0)
   */
  if (i < half && j < half) {
    return mul_basis(half, i, j);
  }

  /*
   * (a,0)(0,d) = (0,da)
   */
  if (i < half && j >= half) {
    size_t jj = j - half;

    MulResult r = mul_basis(half, jj, i);

    r.idx += half;

    return r;
  }

  /*
   * (0,b)(c,0)
   * = (0,b conj(c))
   */
  if (i >= half && j < half) {
    size_t ii = i - half;

    int conj_sign = (j == 0) ? 1 : -1;

    MulResult r = mul_basis(half, ii, j);

    r.idx += half;
    r.sign *= conj_sign;

    return r;
  }

  /*
   * (0,b)(0,d)
   * = (-conj(d)b,0)
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
/* HSL -> RGB                                                     */
/* ============================================================ */

static void hsl_to_rgb(double h, double s, double l, uint8_t *out_r,
                       uint8_t *out_g, uint8_t *out_b) {
  h = mod_positive(h, 1.0);

  s = clamp01(s);

  l = clamp01(l);

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
/* Palette functions                                             */
/* ============================================================ */

static void palette_rainbow(double index, int sign, const VisualConfig *cfg,
                            uint8_t *r, uint8_t *g, uint8_t *b) {
  double hue = index + cfg->hue_shift;

  double lightness = 0.50;

  /*
   * Sign encoded by brightness.
   */
  if (cfg->sign_mapping == SIGN_LIGHTNESS) {
    lightness = sign < 0 ? 0.30 : 0.58;
  }

  /*
   * Sign encoded by hue.
   */
  if (cfg->sign_mapping == SIGN_HUE) {
    hue += sign < 0 ? 0.5 : 0.0;
  }

  /*
   * Sign encoded by saturation.
   */
  double saturation = cfg->saturation;

  if (cfg->sign_mapping == SIGN_SATURATION) {
    saturation = sign < 0 ? saturation * 0.45 : saturation;
  }

  hsl_to_rgb(hue, saturation, lightness, r, g, b);
}

static void palette_ocean(double index, int sign, const VisualConfig *cfg,
                          uint8_t *r, uint8_t *g, uint8_t *b) {
  /*
   * Blue-cyan range.
   */
  double hue = lerp(0.48, 0.62, index);

  if (cfg->sign_mapping == SIGN_HUE) {
    hue += sign < 0 ? 0.10 : 0.0;
  }

  double lightness = cfg->lightness;

  if (cfg->sign_mapping == SIGN_LIGHTNESS) {
    lightness = sign < 0 ? lightness * 0.55 : lightness;
  }

  double saturation = cfg->saturation;

  if (cfg->sign_mapping == SIGN_SATURATION) {
    saturation = sign < 0 ? saturation * 0.40 : saturation;
  }

  hsl_to_rgb(hue, saturation, lightness, r, g, b);
}

static void palette_fire(double index, int sign, const VisualConfig *cfg,
                         uint8_t *r, uint8_t *g, uint8_t *b) {
  double hue = lerp(0.02, 0.15, index);

  if (cfg->sign_mapping == SIGN_HUE) {
    hue += sign < 0 ? 0.50 : 0.0;
  }

  double lightness = cfg->lightness;

  if (cfg->sign_mapping == SIGN_LIGHTNESS) {
    lightness = sign < 0 ? 0.28 : 0.58;
  }

  hsl_to_rgb(hue, cfg->saturation, lightness, r, g, b);
}

static void palette_plasma(double index, int sign, const VisualConfig *cfg,
                           uint8_t *r, uint8_t *g, uint8_t *b) {
  /*
   * Purple -> pink -> orange.
   */
  double hue = lerp(0.78, 0.02, index);

  if (cfg->sign_mapping == SIGN_HUE) {
    hue += sign < 0 ? 0.50 : 0.0;
  }

  double lightness = cfg->lightness;

  if (cfg->sign_mapping == SIGN_LIGHTNESS) {
    lightness = sign < 0 ? 0.28 : 0.58;
  }

  hsl_to_rgb(hue, cfg->saturation, lightness, r, g, b);
}

static void palette_monochrome(double index, int sign, const VisualConfig *cfg,
                               uint8_t *r, uint8_t *g, uint8_t *b) {
  double l = lerp(0.18, 0.82, index);

  if (cfg->sign_mapping == SIGN_LIGHTNESS) {
    if (sign < 0)
      l *= 0.55;
  }

  uint8_t value = (uint8_t)lround(clamp01(l) * 255.0);

  *r = value;
  *g = value;
  *b = value;
}

static void palette_diverging(double index, int sign, const VisualConfig *cfg,
                              uint8_t *r, uint8_t *g, uint8_t *b) {
  double hue;

  /*
   * Blue for positive,
   * red for negative.
   */
  if (sign > 0) {

    hue = lerp(0.55, 0.65, index);

  } else {

    hue = lerp(0.00, 0.08, index);
  }

  double lightness = cfg->lightness;

  if (cfg->sign_mapping == SIGN_LIGHTNESS) {
    lightness = sign > 0 ? 0.58 : 0.42;
  }

  hsl_to_rgb(hue, cfg->saturation, lightness, r, g, b);
}

static void palette_custom(double index, int sign, const VisualConfig *cfg,
                           uint8_t *r, uint8_t *g, uint8_t *b) {
  double hue = sign > 0 ? cfg->custom_hue_positive : cfg->custom_hue_negative;

  /*
   * Index controls hue if explicitly selected.
   */
  if (cfg->index_mapping == INDEX_TO_HUE) {
    hue += index * 0.15;
  }

  double saturation = cfg->custom_saturation;

  if (cfg->index_mapping == INDEX_TO_SATURATION) {
    saturation = lerp(0.15, cfg->custom_saturation, index);
  }

  double lightness = sign > 0 ? cfg->custom_lightness_positive
                              : cfg->custom_lightness_negative;

  if (cfg->index_mapping == INDEX_TO_LIGHTNESS) {
    lightness = lerp(0.20, 0.80, index);
  }

  hsl_to_rgb(hue, saturation, lightness, r, g, b);
}

/* ============================================================ */
/* Mathematical value -> color                                  */
/* ============================================================ */

static void color_for_cell(size_t dim, size_t idx, int sign,
                           const VisualConfig *cfg, uint8_t *r, uint8_t *g,
                           uint8_t *b) {
  /*
   * Normalize basis index.
   */
  double normalized = dim <= 1 ? 0.0 : (double)idx / (double)(dim - 1);

  /*
   * Apply the selected index mapping.
   *
   * For most palettes Hue is the default.
   *
   * The palette functions additionally
   * interpret this normalized quantity.
   */
  if (cfg->index_mapping == INDEX_TO_SATURATION) {
    /*
     * Temporarily encode it through
     * the global saturation parameter.
     */
    VisualConfig local = *cfg;

    local.saturation = lerp(0.15, cfg->saturation, normalized);

    switch (cfg->palette) {

    case PALETTE_RAINBOW:
      palette_rainbow(normalized, sign, &local, r, g, b);
      break;

    case PALETTE_OCEAN:
      palette_ocean(normalized, sign, &local, r, g, b);
      break;

    case PALETTE_FIRE:
      palette_fire(normalized, sign, &local, r, g, b);
      break;

    case PALETTE_PLASMA:
      palette_plasma(normalized, sign, &local, r, g, b);
      break;

    case PALETTE_MONOCHROME:
      palette_monochrome(normalized, sign, &local, r, g, b);
      break;

    case PALETTE_DIVERGING:
      palette_diverging(normalized, sign, &local, r, g, b);
      break;

    case PALETTE_CUSTOM:
      palette_custom(normalized, sign, &local, r, g, b);
      break;
    }

    return;
  }

  if (cfg->index_mapping == INDEX_TO_LIGHTNESS) {
    VisualConfig local = *cfg;

    local.lightness = lerp(0.20, 0.80, normalized);

    switch (cfg->palette) {

    case PALETTE_RAINBOW:
      palette_rainbow(normalized, sign, &local, r, g, b);
      break;

    case PALETTE_OCEAN:
      palette_ocean(normalized, sign, &local, r, g, b);
      break;

    case PALETTE_FIRE:
      palette_fire(normalized, sign, &local, r, g, b);
      break;

    case PALETTE_PLASMA:
      palette_plasma(normalized, sign, &local, r, g, b);
      break;

    case PALETTE_MONOCHROME:
      palette_monochrome(normalized, sign, &local, r, g, b);
      break;

    case PALETTE_DIVERGING:
      palette_diverging(normalized, sign, &local, r, g, b);
      break;

    case PALETTE_CUSTOM:
      palette_custom(normalized, sign, &local, r, g, b);
      break;
    }

    return;
  }

  /*
   * Default:
   * index -> hue
   */
  switch (cfg->palette) {

  case PALETTE_RAINBOW:

    palette_rainbow(normalized, sign, cfg, r, g, b);

    break;

  case PALETTE_OCEAN:

    palette_ocean(normalized, sign, cfg, r, g, b);

    break;

  case PALETTE_FIRE:

    palette_fire(normalized, sign, cfg, r, g, b);

    break;

  case PALETTE_PLASMA:

    palette_plasma(normalized, sign, cfg, r, g, b);

    break;

  case PALETTE_MONOCHROME:

    palette_monochrome(normalized, sign, cfg, r, g, b);

    break;

  case PALETTE_DIVERGING:

    palette_diverging(normalized, sign, cfg, r, g, b);

    break;

  case PALETTE_CUSTOM:

    palette_custom(normalized, sign, cfg, r, g, b);

    break;
  }
}

/* ============================================================ */
/* Pixel operations                                              */
/* ============================================================ */

static void darken_pixel(uint8_t *r, uint8_t *g, uint8_t *b, double factor) {
  factor = clamp01(factor);

  *r = (uint8_t)lround((double)(*r) * factor);

  *g = (uint8_t)lround((double)(*g) * factor);

  *b = (uint8_t)lround((double)(*b) * factor);
}

/* ============================================================ */
/* Recursive hierarchy                                           */
/* ============================================================ */

static size_t distance_to_boundary(size_t coord, size_t block_px) {
  if (block_px == 0)
    return SIZE_MAX;

  size_t rem = coord % block_px;

  size_t left = rem;

  size_t right = block_px - rem;

  return left < right ? left : right;
}

static void apply_recursive_lines(size_t x, size_t y, size_t dim, size_t cell,
                                  uint32_t n, const VisualConfig *cfg,
                                  uint8_t *r, uint8_t *g, uint8_t *b) {
  if (!cfg->recursive_lines || cell == 0 || dim == 0 || n == 0) {
    return;
  }

  for (uint32_t level = 1; level <= n; level++) {

    size_t block_cells = dim >> level;

    if (block_cells == 0) {
      break;
    }

    if (block_cells > SIZE_MAX / cell) {
      break;
    }

    size_t block_px = block_cells * cell;

    size_t dx = distance_to_boundary(x, block_px);

    size_t dy = distance_to_boundary(y, block_px);

    /*
     * Deep hierarchy:
     * finer line.
     */
    size_t line_width;

    if (block_px >= 64)
      line_width = 3;
    else if (block_px >= 16)
      line_width = 2;
    else
      line_width = 1;

    if (dx >= line_width && dy >= line_width) {
      continue;
    }

    /*
     * Stronger for coarse levels.
     */
    double strength = cfg->recursive_strength * pow(0.68, (double)(level - 1));

    darken_pixel(r, g, b, 1.0 - strength);
  }
}

/* ============================================================ */
/* Vignette                                                      */
/* ============================================================ */

static void apply_vignette(size_t x, size_t y, size_t width, size_t height,
                           const VisualConfig *cfg, uint8_t *r, uint8_t *g,
                           uint8_t *b) {
  if (!cfg->vignette || width <= 1 || height <= 1) {
    return;
  }

  double cx = ((double)width - 1.0) / 2.0;

  double cy = ((double)height - 1.0) / 2.0;

  double dx = ((double)x - cx) / cx;

  double dy = ((double)y - cy) / cy;

  double radius = sqrt(dx * dx + dy * dy) / sqrt(2.0);

  radius = clamp01(radius);

  double smooth = radius * radius * (3.0 - 2.0 * radius);

  double factor = 1.0 - cfg->vignette_strength * smooth;

  darken_pixel(r, g, b, factor);
}

/* ============================================================ */
/* BMP low-level writers                                        */
/* ============================================================ */

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

/* ============================================================ */
/* BMP writer                                                    */
/* ============================================================ */

static int write_bmp(const char *path, size_t width, size_t height, size_t dim,
                     size_t cell, uint32_t n, const VisualConfig *cfg) {
  if (width > INT32_MAX || height > INT32_MAX) {
    fprintf(stderr, "image dimensions exceed BMP limits\n");

    return 0;
  }

  if (width > (SIZE_MAX - 3) / 3) {
    fprintf(stderr, "row size overflow\n");

    return 0;
  }

  size_t row_size = (width * 3 + 3) / 4 * 4;

  if (height != 0 && row_size > SIZE_MAX / height) {
    fprintf(stderr, "image size overflow\n");

    return 0;
  }

  size_t pixel_array_size = row_size * height;

  if (pixel_array_size > UINT32_MAX - 54) {
    fprintf(stderr, "BMP file exceeds 4 GiB\n");

    return 0;
  }

  FILE *f = fopen(path, "wb");

  if (f == NULL) {
    perror("failed to create output");

    return 0;
  }

  /*
   * BITMAPFILEHEADER
   */
  if (fwrite("BM", 1, 2, f) != 2 ||

      !write_u32_le(f, (uint32_t)(54 + pixel_array_size)) ||

      !write_u16_le(f, 0) || !write_u16_le(f, 0) ||

      !write_u32_le(f, 54)) {
    fclose(f);
    return 0;
  }

  /*
   * BITMAPINFOHEADER
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

  uint8_t *row = calloc(row_size, 1);

  if (row == NULL) {
    fprintf(stderr, "failed to allocate row buffer\n");

    fclose(f);

    return 0;
  }

  /*
   * Bottom-up BMP.
   */
  for (size_t y = height; y-- > 0;) {

    memset(row, 0, row_size);

    for (size_t x = 0; x < width; x++) {

      size_t i = y / cell;

      size_t j = x / cell;

      MulResult result = mul_basis(dim, i, j);

      uint8_t r;
      uint8_t g;
      uint8_t b;

      /*
       * Mathematical value
       * -> color.
       */
      color_for_cell(dim, result.idx, result.sign, cfg, &r, &g, &b);

      /*
       * Visual layer 1:
       * recursive structure.
       */
      apply_recursive_lines(x, y, dim, cell, n, cfg, &r, &g, &b);

      /*
       * Visual layer 2:
       * vignette.
       */
      apply_vignette(x, y, width, height, cfg, &r, &g, &b);

      /*
       * BGR order.
       */
      size_t offset = x * 3;

      row[offset] = b;

      row[offset + 1] = g;

      row[offset + 2] = r;
    }

    if (fwrite(row, 1, row_size, f) != row_size) {
      fprintf(stderr, "failed to write row\n");

      free(row);
      fclose(f);

      return 0;
    }
  }

  free(row);

  if (fclose(f) != 0) {
    perror("failed to close file");

    return 0;
  }

  return 1;
}

/* ============================================================ */
/* Size formatting                                               */
/* ============================================================ */

static void human_bytes(double bytes, char *buffer, size_t buffer_size) {
  const char *units[] = {"B", "KB", "MB", "GB", "TB"};

  size_t unit = 0;

  while (bytes >= 1024.0 && unit < 4) {
    bytes /= 1024.0;
    unit++;
  }

  snprintf(buffer, buffer_size, "%.1f%s", bytes, units[unit]);
}

/* ============================================================ */
/* Visualization settings                                       */
/* ============================================================ */

static const char *palette_name(Palette p) {
  switch (p) {

  case PALETTE_RAINBOW:
    return "Rainbow";

  case PALETTE_OCEAN:
    return "Ocean";

  case PALETTE_FIRE:
    return "Fire";

  case PALETTE_PLASMA:
    return "Plasma";

  case PALETTE_MONOCHROME:
    return "Monochrome";

  case PALETTE_DIVERGING:
    return "Diverging";

  case PALETTE_CUSTOM:
    return "Custom";

  default:
    return "Unknown";
  }
}

static const char *index_mapping_name(IndexMapping m) {
  switch (m) {

  case INDEX_TO_HUE:
    return "Hue";

  case INDEX_TO_SATURATION:
    return "Saturation";

  case INDEX_TO_LIGHTNESS:
    return "Lightness";

  default:
    return "Unknown";
  }
}

static const char *sign_mapping_name(SignMapping m) {
  switch (m) {

  case SIGN_LIGHTNESS:
    return "Lightness";

  case SIGN_HUE:
    return "Hue";

  case SIGN_SATURATION:
    return "Saturation";

  default:
    return "Unknown";
  }
}

static void print_config(const VisualConfig *cfg) {
  printf("\n");
  printf("========================================\n");
  printf("           VISUAL CONFIGURATION\n");
  printf("========================================\n");

  printf("Palette            : %s\n", palette_name(cfg->palette));

  printf("Index ->           : %s\n", index_mapping_name(cfg->index_mapping));

  printf("Sign ->             : %s\n", sign_mapping_name(cfg->sign_mapping));

  printf("Hue shift           : %.3f\n", cfg->hue_shift);

  printf("Saturation          : %.3f\n", cfg->saturation);

  printf("Lightness           : %.3f\n", cfg->lightness);

  printf("Recursive lines     : %s\n", cfg->recursive_lines ? "ON" : "OFF");

  printf("Recursive strength  : %.3f\n", cfg->recursive_strength);

  printf("Vignette            : %s\n", cfg->vignette ? "ON" : "OFF");

  printf("Vignette strength   : %.3f\n", cfg->vignette_strength);

  printf("========================================\n");
}

/* ============================================================ */
/* Visualization menu                                            */
/* ============================================================ */

static void configure_visuals(VisualConfig *cfg) {
  for (;;) {

    print_config(cfg);

    printf("\n"
           "1. Palette\n"
           "2. Index mapping\n"
           "3. Sign mapping\n"
           "4. Hue shift\n"
           "5. Saturation\n"
           "6. Lightness\n"
           "7. Recursive lines\n"
           "8. Recursive strength\n"
           "9. Vignette\n"
           "10. Vignette strength\n"
           "11. Custom palette\n"
           "0. Done\n");

    int choice = read_int("Select", 0);

    switch (choice) {

    case 0:
      return;

    case 1: {

      printf("\n"
             "Palette\n"
             "1. Rainbow\n"
             "2. Ocean\n"
             "3. Fire\n"
             "4. Plasma\n"
             "5. Monochrome\n"
             "6. Diverging\n"
             "7. Custom\n");

      int p = read_int("Select", (int)cfg->palette);

      if (p >= 1 && p <= 7) {
        cfg->palette = (Palette)p;
      }

      break;
    }

    case 2: {

      printf("\n"
             "Index mapping\n"
             "1. Hue\n"
             "2. Saturation\n"
             "3. Lightness\n");

      int m = read_int("Select", (int)cfg->index_mapping);

      if (m >= 1 && m <= 3) {
        cfg->index_mapping = (IndexMapping)m;
      }

      break;
    }

    case 3: {

      printf("\n"
             "Sign mapping\n"
             "1. Lightness\n"
             "2. Hue\n"
             "3. Saturation\n");

      int m = read_int("Select", (int)cfg->sign_mapping);

      if (m >= 1 && m <= 3) {
        cfg->sign_mapping = (SignMapping)m;
      }

      break;
    }

    case 4:

      cfg->hue_shift = read_double("Hue shift (0..1)", cfg->hue_shift);

      break;

    case 5:

      cfg->saturation =
          clamp01(read_double("Saturation (0..1)", cfg->saturation));

      break;

    case 6:

      cfg->lightness = clamp01(read_double("Lightness (0..1)", cfg->lightness));

      break;

    case 7:

      cfg->recursive_lines =
          read_int("Recursive lines (0=off, 1=on)", cfg->recursive_lines) ? 1
                                                                          : 0;

      break;

    case 8:

      cfg->recursive_strength = clamp01(
          read_double("Recursive strength (0..1)", cfg->recursive_strength));

      break;

    case 9:

      cfg->vignette = read_int("Vignette (0=off, 1=on)", cfg->vignette) ? 1 : 0;

      break;

    case 10:

      cfg->vignette_strength = clamp01(
          read_double("Vignette strength (0..1)", cfg->vignette_strength));

      break;

    case 11:

      printf("\n"
             "Custom HSL palette\n");

      cfg->custom_hue_positive = mod_positive(
          read_double("Positive hue (0..1)", cfg->custom_hue_positive), 1.0);

      cfg->custom_hue_negative = mod_positive(
          read_double("Negative hue (0..1)", cfg->custom_hue_negative), 1.0);

      cfg->custom_saturation = clamp01(
          read_double("Custom saturation (0..1)", cfg->custom_saturation));

      cfg->custom_lightness_positive = clamp01(read_double(
          "Positive lightness (0..1)", cfg->custom_lightness_positive));

      cfg->custom_lightness_negative = clamp01(read_double(
          "Negative lightness (0..1)", cfg->custom_lightness_negative));

      break;

    default:

      printf("Unknown option.\n");

      break;
    }
  }
}

/* ============================================================ */
/* Default configuration                                         */
/* ============================================================ */

static VisualConfig default_visual_config(void) {
  VisualConfig cfg;

  cfg.palette = PALETTE_RAINBOW;

  cfg.index_mapping = INDEX_TO_HUE;

  cfg.sign_mapping = SIGN_LIGHTNESS;

  cfg.hue_shift = 0.0;

  cfg.saturation = 0.65;

  cfg.lightness = 0.55;

  cfg.contrast = 1.0;

  cfg.recursive_lines = 1;

  cfg.recursive_strength = 0.18;

  cfg.vignette = 1;

  cfg.vignette_strength = 0.25;

  cfg.custom_hue_positive = 0.58;

  cfg.custom_hue_negative = 0.02;

  cfg.custom_saturation = 0.70;

  cfg.custom_lightness_positive = 0.60;

  cfg.custom_lightness_negative = 0.38;

  return cfg;
}

/* ============================================================ */
/* Generate                                                       */
/* ============================================================ */

static int generate(uint32_t n, size_t cell, const char *output,
                    const VisualConfig *cfg) {
  if (n >= sizeof(size_t) * CHAR_BIT - 1) {
    fprintf(stderr, "n=%u is too large for size_t\n", n);

    return 1;
  }

  size_t dim = (size_t)1 << n;

  if (cell == 0) {
    fprintf(stderr, "cell size must be > 0\n");

    return 1;
  }

  if (dim > SIZE_MAX / cell) {
    fprintf(stderr, "image dimensions overflow\n");

    return 1;
  }

  size_t width = dim * cell;

  size_t height = dim * cell;

  if (width > (SIZE_MAX - 3) / 3) {
    fprintf(stderr, "row size overflow\n");

    return 1;
  }

  size_t row_size = (width * 3 + 3) / 4 * 4;

  if (height != 0 && row_size > SIZE_MAX / height) {
    fprintf(stderr, "image size overflow\n");

    return 1;
  }

  double estimated = (double)row_size * (double)height;

  char size_text[64];

  human_bytes(estimated, size_text, sizeof(size_text));

  printf("\n"
         "========================================\n"
         " Cayley-Dickson Visualization\n"
         "========================================\n"
         "n       = %u\n"
         "dimension = %zu\n"
         "cells   = %zux%zu\n"
         "image   = %zux%zu px\n"
         "size    = ~%s\n"
         "palette = %s\n"
         "========================================\n",
         n, dim, dim, dim, width, height, size_text,
         palette_name(cfg->palette));

  if (dim <= 16) {
    printf("\nSmall algebra detected.\n"
           "Text multiplication table:\n\n");

    /*
     * Compact table.
     */
    char labels[16][32];

    size_t label_width = 2;

    for (size_t i = 0; i < dim; i++) {

      if (i == 0) {

        strcpy(labels[i], "1");

      } else if (dim <= 4) {

        const char *s = i == 1 ? "i" : i == 2 ? "j" : "k";

        strcpy(labels[i], s);

      } else {

        snprintf(labels[i], sizeof(labels[i]), "e%zu", i);
      }

      size_t len = strlen(labels[i]);

      if (len > label_width) {
        label_width = len;
      }
    }

    label_width++;

    printf("%*s", (int)(label_width + 1), "");

    for (size_t j = 0; j < dim; j++) {

      printf("%*s", (int)label_width, labels[j]);
    }

    printf("\n");

    for (size_t i = 0; i < dim; i++) {

      printf("%*s ", (int)label_width, labels[i]);

      for (size_t j = 0; j < dim; j++) {

        MulResult r = mul_basis(dim, i, j);

        if (r.sign < 0) {

          printf("%*s", (int)label_width, "-?");

        } else {

          printf("%*s", (int)label_width, labels[r.idx]);
        }
      }

      printf("\n");
    }
  }

  printf("\nGenerating...\n");

  clock_t start = clock();

  int ok = write_bmp(output, width, height, dim, cell, n, cfg);

  clock_t end = clock();

  if (!ok) {
    return 1;
  }

  double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

  if (elapsed < 1e-9) {
    elapsed = 1e-9;
  }

  double mpix = ((double)width * (double)height) / 1e6 / elapsed;

  printf("\n"
         "Generated: %s\n"
         "Time     : %.3f s\n"
         "Speed    : %.2f Mpixels/s\n",
         output, elapsed, mpix);

  return 0;
}

/* ============================================================ */
/* Main                                                          */
/* ============================================================ */

int main(int argc, char **argv) {
  uint32_t n = 6;

  size_t cell = 4;

  char output[512] = "cayley_table.bmp";

  /*
   * Command-line parameters.
   */
  if (argc >= 2) {

    char *end;

    unsigned long value = strtoul(argv[1], &end, 10);

    if (end == argv[1] || *end != '\0' || value > UINT32_MAX) {
      fprintf(stderr, "invalid n: %s\n", argv[1]);

      return 1;
    }

    n = (uint32_t)value;
  }

  if (argc >= 3) {

    char *end;

    unsigned long long value = strtoull(argv[2], &end, 10);

    if (end == argv[2] || *end != '\0' || value == 0 || value > SIZE_MAX) {
      fprintf(stderr, "invalid cell size: %s\n", argv[2]);

      return 1;
    }

    cell = (size_t)value;
  }

  if (argc >= 4) {
    snprintf(output, sizeof(output), "%s", argv[3]);
  }

  VisualConfig cfg = default_visual_config();

  /*
   * No command-line arguments:
   * interactive mode.
   */
  if (argc == 1) {

    printf("Cayley-Dickson Mathematical Visualizer\n");

    char buffer[512];

    printf("\n"
           "Generation parameters\n");

    char n_text[32];

    snprintf(n_text, sizeof(n_text), "%u", n);

    read_string("n (dimension = 2^n)", n_text, buffer, sizeof(buffer));

    char *end;

    unsigned long parsed_n = strtoul(buffer, &end, 10);

    if (end != buffer && *end == '\0' && parsed_n <= UINT32_MAX) {
      n = (uint32_t)parsed_n;
    }

    char cell_text[32];

    snprintf(cell_text, sizeof(cell_text), "%zu", cell);

    read_string("pixels per cell", cell_text, buffer, sizeof(buffer));

    unsigned long long parsed_cell = strtoull(buffer, &end, 10);

    if (end != buffer && *end == '\0' && parsed_cell > 0 &&
        parsed_cell <= SIZE_MAX) {
      cell = (size_t)parsed_cell;
    }

    read_string("output filename", output, output, sizeof(output));

    /*
     * Visualization configuration.
     */
    configure_visuals(&cfg);
  }

  /*
   * CLI mode still allows an interactive
   * visualization configuration after
   * providing the numeric parameters.
   *
   * Ask explicitly so scripts using arguments
   * remain predictable.
   */
  if (argc > 1) {

    int interactive = read_int("Open visualization menu? (0=no, 1=yes)", 1);

    if (interactive) {
      configure_visuals(&cfg);
    }
  }

  /*
   * Generate.
   */
  return generate(n, cell, output, &cfg);
}
