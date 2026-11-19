/* SoX Resampler Library      Copyright (c) 2007-18 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

/* Exercises the packed 24-bit PCM (3 bytes/sample) I/O datatype: full-scale
 * scaling, rounding, saturation/clip counting, interleaved vs. split channel
 * layouts, oneshot vs. streaming use, variable-rate, and interchange with the
 * existing datatypes. */

#include <soxr.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures;

#define CHECK(cond, ...) do { if (!(cond)) { \
  fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); \
  fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); ++failures; } } while (0)

#define IRATE 96000.
#define ORATE 44100.

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static long load_s24(unsigned char const * p)
{
  /* Two's-complement sign extension from 24 bits (host byte order, LE). */
  long x = (long)((unsigned long)p[0] | (unsigned long)p[1] << 8 |
      (unsigned long)p[2] << 16);
  return (x ^ 0x800000L) - 0x800000L;
}

static void store_s24(unsigned char * p, long x)
{
  /* Two's-complement truncation to 24 bits (host byte order, LE). */
  unsigned long u = (unsigned long)(int32_t)x;
  p[0] = (unsigned char)u; p[1] = (unsigned char)(u >> 8);
  p[2] = (unsigned char)(u >> 16);
}

static void gen_s24(unsigned char * buf, size_t n, unsigned ch, int interleaved)
{
  size_t j; unsigned i;
  for (j = 0; j < n; ++j) for (i = 0; i < ch; ++i) {
    long v = (long)(8000000. * sin(2. * M_PI * (j * (i + 1)) / 67.));
    if (interleaved) store_s24(buf + 3 * (j * ch + i), v);
    else store_s24(buf + 3 * (i * n + j), v);
  }
}

static int eq_s24(unsigned char const * a, unsigned char const * b, size_t bytes)
{
  size_t i;
  for (i = 0; i < bytes; ++i) if (a[i] != b[i]) return 0;
  return 1;
}

/* Advance an output buffer by `total' samples per channel, accounting for the
 * split-channel pointer-array layout when needed. */
static void * advance_out(soxr_datatype_t ot, unsigned ch, void * out,
    size_t total)
{
  unsigned i;
  if (ot & SOXR_SPLIT) {
    void ** ptrs = (void **)out;
    for (i = 0; i < ch; ++i)
      ptrs[i] = (unsigned char *)ptrs[i] + 3 * total;
    return ptrs;
  }
  return (unsigned char *)out + 3 * total * ch;
}

static size_t convert(soxr_datatype_t it, soxr_datatype_t ot, unsigned ch,
    void const * in, size_t ilen, void * out, size_t olen,
    soxr_quality_spec_t const * q, unsigned long qflags)
{
  soxr_io_spec_t io = soxr_io_spec(it, ot);
  soxr_quality_spec_t qs;
  soxr_error_t e; soxr_t r;
  size_t idone = 0, odone = 0, total = 0, guard = 0;
  CHECK(!io.e, "soxr_io_spec error");
  if (q) qs = *q; else qs = soxr_quality_spec(SOXR_LQ, qflags);
  r = soxr_create(IRATE, ORATE, ch, &e, &io, &qs, NULL);
  CHECK(!e && r, "soxr_create: %s", soxr_strerror(e));
  if (!r) return 0;
  e = soxr_process(r, (void *)in, ilen, &idone, out, olen, &odone);
  CHECK(!e, "soxr_process: %s", soxr_strerror(e));
  CHECK(idone == ilen, "idone %lu != ilen %lu", (unsigned long)idone,
      (unsigned long)ilen);
  total = odone;
  do {
    void * op = advance_out(ot, ch, out, total);
    e = soxr_process(r, NULL, 0, NULL, op, olen - total, &odone);
    CHECK(!e, "flush: %s", soxr_strerror(e));
    total += odone;
    CHECK(total <= olen, "output overflow %lu > %lu", (unsigned long)total,
        (unsigned long)olen);
  } while (odone && ++guard < 100000);
  CHECK(guard < 100000, "flush did not terminate");
  soxr_delete(r);
  return total;
}

static size_t convert_clips(soxr_datatype_t it, soxr_datatype_t ot,
    unsigned ch, void const * in, size_t ilen, void * out, size_t olen,
    soxr_quality_spec_t const * q, unsigned long qflags, size_t * clips)
{
  soxr_io_spec_t io = soxr_io_spec(it, ot);
  soxr_quality_spec_t qs;
  soxr_error_t e; soxr_t r;
  size_t idone = 0, odone = 0, total = 0, guard = 0;
  if (q) qs = *q; else qs = soxr_quality_spec(SOXR_LQ, qflags);
  r = soxr_create(IRATE, ORATE, ch, &e, &io, &qs, NULL);
  CHECK(!e && r, "soxr_create: %s", soxr_strerror(e));
  if (!r) return 0;
  e = soxr_process(r, (void *)in, ilen, &idone, out, olen, &odone);
  CHECK(!e, "soxr_process: %s", soxr_strerror(e));
  total = odone;
  do {
    void * op = advance_out(ot, ch, out, total);
    e = soxr_process(r, NULL, 0, NULL, op, olen - total, &odone);
    CHECK(!e, "flush: %s", soxr_strerror(e));
    total += odone;
  } while (odone && ++guard < 100000);
  *clips = *soxr_num_clips(r);
  soxr_delete(r);
  return total;
}

static void test_datatype_macros(void)
{
  CHECK(soxr_datatype_size(SOXR_INT24_I) == 3, "INT24_I size");
  CHECK(soxr_datatype_size(SOXR_INT24_S) == 3, "INT24_S size");
  /* Binary compatibility of the existing values: */
  CHECK(soxr_datatype_size(SOXR_FLOAT32_I) == 4, "FLOAT32 size");
  CHECK(soxr_datatype_size(SOXR_FLOAT64_I) == 8, "FLOAT64 size");
  CHECK(soxr_datatype_size(SOXR_INT32_I) == 4, "INT32 size");
  CHECK(soxr_datatype_size(SOXR_INT16_I) == 2, "INT16 size");
  CHECK(SOXR_INT24_I == 8, "SOXR_INT24_I value");
  CHECK(SOXR_INT24_S == 12, "SOXR_INT24_S value");
  CHECK(SOXR_INT16_S == 7, "SOXR_INT16_S value unchanged");
}

static void test_invalid_datatype(void)
{
  soxr_io_spec_t io = soxr_io_spec((soxr_datatype_t)9, SOXR_INT24_I);
  CHECK(io.e != 0, "type 9 must be rejected");
  io = soxr_io_spec(SOXR_INT24_I, (soxr_datatype_t)13);
  CHECK(io.e != 0, "type 13 must be rejected");
  io = soxr_io_spec(SOXR_INT24_I, SOXR_INT24_S);
  CHECK(io.e == 0, "valid 24-bit types rejected");
}

static void test_roundtrip(void)
{
  /* 1:1 rate, 24-bit in/out.  The unity-ratio FIR is designed to be a
   * pass-through; allow at most a 1-LSB rounding wobble on a handful of
   * samples while requiring exact length and zero clipping. */
  size_t n = 4096, i;
  unsigned char * in = malloc(3 * n), * out = malloc(3 * (n + 256));
  soxr_io_spec_t io = soxr_io_spec(SOXR_INT24_I, SOXR_INT24_I);
  soxr_error_t e; soxr_t r; size_t idone = 0, odone = 0, total = 0, guard = 0;
  size_t wobble = 0;
  for (i = 0; i < n; ++i)
    store_s24(in + 3 * i, (long)(8000000. * sin(2. * M_PI * i / 31.)));
  r = soxr_create(1, 1, 1, &e, &io, NULL, NULL);
  CHECK(!e && r, "roundtrip create: %s", soxr_strerror(e));
  if (r) {
    e = soxr_process(r, in, n, &idone, out, n + 256, &odone);
    CHECK(!e, "roundtrip process: %s", soxr_strerror(e));
    CHECK(idone == n, "roundtrip idone %lu != %lu", (unsigned long)idone,
        (unsigned long)n);
    total = odone;
    do {
      e = soxr_process(r, NULL, 0, NULL, out + 3 * total,
          n + 256 - total, &odone);
      CHECK(!e, "roundtrip flush: %s", soxr_strerror(e));
      total += odone;
    } while (odone && ++guard < 100000);
    CHECK(total == n, "roundtrip total %lu != %lu", (unsigned long)total,
        (unsigned long)n);
    for (i = 0; i < n; ++i) {
      long diff = load_s24(out + 3 * i) - load_s24(in + 3 * i);
      if (diff < 0) diff = -diff;
      if (diff == 1) ++wobble;
      else CHECK(diff == 0, "roundtrip sample %lu differs by %ld",
          (unsigned long)i, diff);
    }
    CHECK(wobble <= 8, "too many 1-LSB wobble samples: %lu",
        (unsigned long)wobble);
    CHECK(*soxr_num_clips(r) == 0, "roundtrip clips");
    soxr_delete(r);
  }
  free(in); free(out);
}

static void test_rounding_and_clipping(void)
{
  /* Float32 -> packed 24-bit at unity rate.  Each level is held for a long
   * plateau so the resampling filter settles to the exact level (transients at
   * plateau edges are ignored).  The engine scales by 2^23, hence the /2^23
   * normalization; output is rounded to nearest (ties away from zero) and
   * out-of-range levels saturate and count as clips. */
  struct { float level; long code; int clip; } const table[] = {
    { 0.25f / 8388608.f,             0, 0},  /* below half LSB   -> 0    */
    { 0.5f  / 8388608.f,             1, 0},  /* +0.5 LSB tie     -> +1   */
    {-0.5f  / 8388608.f,            -1, 0},  /* -0.5 LSB tie     -> -1   */
    { 100.5f / 8388608.f,          101, 0},  /* interior tie     -> +101 */
    { 8388607.0f / 8388608.f,  8388607, 0},  /* +max code        */
    {-8388608.0f / 8388608.f, -8388608, 0},  /* -max code        */
    { 8388608.0f / 8388608.f,  8388607, 1},  /* +2^23 over-range -> clip */
    {-8388609.0f / 8388608.f, -8388608, 1},  /* below -2^23      -> clip */
  };
  enum {PLATEAU = 4096, GUARD = 512};
  size_t nlevels = sizeof(table) / sizeof(table[0]);
  size_t ilen = nlevels * PLATEAU, i, m, total;
  size_t cap = ilen + 256;
  float * in = malloc(4 * ilen);
  unsigned char * out = malloc(3 * cap);
  size_t clips = 0; int bad = 0;
  for (m = 0; m < nlevels; ++m)
    for (i = 0; i < PLATEAU; ++i) in[m * PLATEAU + i] = table[m].level;

  total = convert_clips(SOXR_FLOAT32_I, SOXR_INT24_I, 1, in, ilen, out, cap,
      NULL, 0, &clips);

  for (m = 0; m < nlevels && !bad; ++m) {
    long code = table[m].code;
    /* Sample the settled middle of each plateau (skip edges on both sides). */
    size_t lo = m * PLATEAU + GUARD, hi = (m + 1) * PLATEAU - GUARD;
    for (i = lo; i < hi && i < total; ++i) {
      long got = load_s24(out + 3 * i);
      if (got != code) {
        fprintf(stderr,
            "FAIL level %lu (%.8g): sample %lu got %ld, want %ld\n",
            (unsigned long)m, table[m].level, (unsigned long)i, got, code);
        bad = 1; break;
      }
    }
  }
  CHECK(!bad, "rounding/saturation level mismatch");
  /* soxr_num_clips counts individual saturated samples: two plateaus, each of
   * PLATEAU samples (filter transients at the edges may clip as well). */
  CHECK(clips >= 2 * (PLATEAU - 2 * GUARD) && clips <= 2 * PLATEAU,
      "clip count %lu outside expected saturation range",
      (unsigned long)clips);
  free(in); free(out);
}

static void test_layouts_and_engines(unsigned long recipe, unsigned long flags,
    char const * tag)
{
  enum {CH = 2, ILEN = 8000};
  size_t olen = (size_t)(ILEN * ORATE / IRATE) + 256;
  unsigned char * in_i = malloc(3 * (size_t)ILEN * CH);
  unsigned char * in_sb[CH];
  void * in_s_ptrs[CH];
  unsigned char * out_ii = malloc(3 * olen * CH),
      * out_si = malloc(3 * olen * CH);
  unsigned char * out_isb[CH], * out_ssb[CH];
  void * out_is_ptrs[CH], * out_ss_ptrs[CH];
  soxr_quality_spec_t q = soxr_quality_spec(recipe, flags);
  size_t c_ii, c_ss, n_ii, n_ss, n_is, n_si, k; int got = 0; unsigned c;
  in_sb[0] = malloc(3 * ILEN); in_sb[1] = malloc(3 * ILEN);
  out_isb[0] = malloc(3 * olen); out_isb[1] = malloc(3 * olen);
  out_ssb[0] = malloc(3 * olen); out_ssb[1] = malloc(3 * olen);
  for (c = 0; c < CH; ++c) {
    in_s_ptrs[c] = in_sb[c];
    out_is_ptrs[c] = out_isb[c];
    out_ss_ptrs[c] = out_ssb[c];
  }
  /* One canonical per-channel signal; both input layouts derive from it. */
  for (k = 0; k < ILEN; ++k) for (c = 0; c < CH; ++c) {
    long v = (long)(8000000. * sin(2. * M_PI * (k * (c + 1)) / 67.));
    store_s24(in_i + 3 * (k * CH + c), v);
    store_s24(in_sb[c] + 3 * k, v);
  }

  /* All four i/o layout combinations must agree: */
  n_ii = convert_clips(SOXR_INT24_I, SOXR_INT24_I, CH, in_i, ILEN,
      out_ii, olen, &q, 0, &c_ii);
  n_is = convert(SOXR_INT24_I, SOXR_INT24_S, CH, in_i, ILEN,
      out_is_ptrs, olen, &q, 0);
  n_si = convert(SOXR_INT24_S, SOXR_INT24_I, CH, in_s_ptrs, ILEN,
      out_si, olen, &q, 0);
  n_ss = convert_clips(SOXR_INT24_S, SOXR_INT24_S, CH, in_s_ptrs, ILEN,
      out_ss_ptrs, olen, &q, 0, &c_ss);
  CHECK(n_ii == n_ss && n_ii == n_is && n_ii == n_si,
      "[%s] output length mismatch %lu/%lu/%lu/%lu", tag,
      (unsigned long)n_ii, (unsigned long)n_is,
      (unsigned long)n_si, (unsigned long)n_ss);
  CHECK(c_ii == c_ss, "[%s] clip count mismatch %lu/%lu", tag,
      (unsigned long)c_ii, (unsigned long)c_ss);
  CHECK(eq_s24(out_ii, out_si, 3 * n_ii * CH),
      "[%s] interleaved outputs differ (input layout)", tag);

  /* Compare split output buffers channel-by-channel against interleaved: */
  for (k = 0; k < n_ii && !got; ++k)
    for (c = 0; c < CH; ++c) {
      long a = load_s24(out_ii + 3 * (k * CH + c));
      long b = load_s24(out_isb[c] + 3 * k);
      long d = load_s24(out_ssb[c] + 3 * k);
      if (a != b || a != d) { got = 1; break; }
    }
  CHECK(!got, "[%s] split/interleaved output mismatch at %lu", tag,
      (unsigned long)k);
  free(in_i); free(in_sb[0]); free(in_sb[1]);
  free(out_ii); free(out_si);
  free(out_isb[0]); free(out_isb[1]); free(out_ssb[0]); free(out_ssb[1]);
}

static void test_streaming(void)
{
  /* Same signal, oneshot vs. streamed in small blocks must be identical. */
  enum {CH = 2, ILEN = 6000};
  size_t olen = (size_t)(ILEN * ORATE / IRATE) + 256;
  unsigned char * in = malloc(3 * (size_t)ILEN * CH),
      * out_one = malloc(3 * olen * CH), * out_str = malloc(3 * olen * CH);
  soxr_io_spec_t io = soxr_io_spec(SOXR_INT24_I, SOXR_INT24_I);
  soxr_error_t e; soxr_t r;
  size_t idone = 0, odone = 0, total = 0, pos = 0, n_one, guard = 0;
  gen_s24(in, ILEN, CH, 1);

  n_one = convert(SOXR_INT24_I, SOXR_INT24_I, CH, in, ILEN,
      out_one, olen, NULL, 0);

  r = soxr_create(IRATE, ORATE, CH, &e, &io, NULL, NULL);
  CHECK(!e && r, "streaming create: %s", soxr_strerror(e));
  if (r) {
    while (pos < ILEN) {
      size_t block = ILEN - pos; if (block > 373) block = 373;
      e = soxr_process(r, in + 3 * pos * CH, block, &idone,
          out_str + 3 * total * CH, olen - total, &odone);
      CHECK(!e, "streaming process: %s", soxr_strerror(e));
      pos += idone; total += odone;
    }
    do {
      e = soxr_process(r, NULL, 0, NULL,
          out_str + 3 * total * CH, olen - total, &odone);
      CHECK(!e, "streaming flush: %s", soxr_strerror(e));
      total += odone;
    } while (odone && ++guard < 100000);
    CHECK(total == n_one, "streaming length %lu != oneshot %lu",
        (unsigned long)total, (unsigned long)n_one);
    CHECK(eq_s24(out_one, out_str, 3 * total * CH),
        "streaming output differs from oneshot (%lu samples)",
        (unsigned long)total);
    soxr_delete(r);
  }
  free(in); free(out_one); free(out_str);
}

static void test_interchange(void)
{
  /* The int32 and packed int24 full-scale conventions must agree.  The same
   * physical level is represented by int32 code c*256 and by int24 code c
   * (full-scale 2^31 vs 2^23); feeding both must yield identical 24-bit
   * output, proving the full-scale mapping in soxr_datatype_full_scale(). */
  enum {ILEN = 4000};
  size_t olen = (size_t)(ILEN * ORATE / IRATE) + 256, i, n1, n2;
  size_t wobble = 0;
  int32_t * in32 = malloc(4 * ILEN);
  unsigned char * in24 = malloc(3 * ILEN),
      * out_from32 = malloc(3 * olen), * out_from24 = malloc(3 * olen);
  for (i = 0; i < ILEN; ++i) {
    long c = (long)(2000000. * sin(2. * M_PI * i / 53.)); /* 24-bit code */
    store_s24(in24 + 3 * i, c);
    in32[i] = (int32_t)(c * 256); /* same level expressed as an int32 code. */
  }
  n1 = convert(SOXR_INT32_I, SOXR_INT24_I, 1, in32, ILEN,
      out_from32, olen, NULL, 0);
  n2 = convert(SOXR_INT24_I, SOXR_INT24_I, 1, in24, ILEN,
      out_from24, olen, NULL, 0);
  CHECK(n1 == n2, "interchange length mismatch %lu/%lu",
      (unsigned long)n1, (unsigned long)n2);
  for (i = 0; i < n1; ++i) {
    long diff = load_s24(out_from32 + 3 * i) - load_s24(out_from24 + 3 * i);
    if (diff < 0) diff = -diff;
    if (diff == 1) ++wobble;
    else CHECK(diff == 0, "interchange sample %lu differs by %ld",
        (unsigned long)i, diff);
  }
  CHECK(wobble <= n1 / 100, "excessive interchange wobble: %lu",
      (unsigned long)wobble);
  free(in32); free(in24); free(out_from32); free(out_from24);
}

static void test_variable_rate(void)
{
  /* VR engine with 24-bit I/O: exercise the SOXR_VR path end-to-end, varying
   * the ratio mid-stream; output must stay inside the 24-bit range. */
  enum {ILEN = 4000};
  double ratio = 1.5;
  size_t olen = (size_t)(ILEN / ratio) + 512, i, total = 0, pos = 0,
      idone = 0, odone = 0, guard = 0;
  unsigned char * in = malloc(3 * ILEN), * out = malloc(3 * olen);
  soxr_quality_spec_t q = soxr_quality_spec(SOXR_HQ, SOXR_VR);
  soxr_io_spec_t io = soxr_io_spec(SOXR_INT24_I, SOXR_INT24_I);
  soxr_error_t e; soxr_t r;
  gen_s24(in, ILEN, 1, 1);

  r = soxr_create(ratio, 1, 1, &e, &io, &q, NULL); /* max ratio */
  CHECK(!e && r, "vr create: %s", soxr_strerror(e));
  if (r) {
    soxr_set_io_ratio(r, ratio, 0);
    while (pos < ILEN) {
      size_t block = ILEN - pos; if (block > 512) block = 512;
      if (pos) soxr_set_io_ratio(r, ratio * (pos % 2? 1.1 : 0.9), 64);
      e = soxr_process(r, in + 3 * pos, block, &idone,
          out + 3 * total, olen - total, &odone);
      CHECK(!e, "vr process: %s", soxr_strerror(e));
      pos += idone; total += odone;
    }
    do {
      e = soxr_process(r, NULL, 0, NULL, out + 3 * total,
          olen - total, &odone);
      CHECK(!e, "vr flush: %s", soxr_strerror(e));
      total += odone;
    } while (odone && ++guard < 100000);
    CHECK(total > 0, "vr produced no output");
    for (i = 0; i < total; ++i) {
      long v = load_s24(out + 3 * i);
      CHECK(v >= -8388608L && v <= 8388607L, "vr out of range at %lu: %ld",
          (unsigned long)i, v);
    }
    soxr_delete(r);
  }
  free(in); free(out);
}

int main(void)
{
  test_datatype_macros();
  test_invalid_datatype();
  test_roundtrip();
  test_rounding_and_clipping();
  test_layouts_and_engines(SOXR_LQ, 0, "LQ/32-bit engine");
  test_layouts_and_engines(SOXR_HQ, 0, "HQ/32-bit engine");
  test_layouts_and_engines(SOXR_VHQ, 0, "VHQ/64-bit engine");
  test_layouts_and_engines(SOXR_VHQ, SOXR_DOUBLE_PRECISION,
      "VHQ/explicit-double");
  test_streaming();
  test_interchange();
  test_variable_rate();
  if (failures) {
    fprintf(stderr, "io24-test: %d failure(s)\n", failures);
    return 1;
  }
  printf("io24-test: all tests passed\n");
  return 0;
}
