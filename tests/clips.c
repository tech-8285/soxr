/* SoX Resampler Library      Copyright (c) 2007-18 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

/* Tests of clipping statistics for integer output.
 *
 * The reconstruction filter can overshoot: a signal whose input samples are
 * all within the legal [-1, 1] range can produce filtered peaks beyond +/-1,
 * and output gain (io_spec.scale) can amplify those further.  Converting such
 * a peak to an integer type must limit the value to the range of that type
 * and count the limited output sample in soxr_num_clips().
 *
 * Checked here:
 *
 *  1. A sample that merely rounds to an extreme value is not counted as a
 *     clip; only a sample that would round *outside* the type's range is
 *     limited and counted.  The same output value (e.g. 32767) can therefore
 *     occur both with and without a clip having taken place.
 *
 *  2. The two's complement ranges are asymmetric, so positive and negative
 *     peaks have different saturation boundaries (int16: +32767 / -32768;
 *     int32: +2147483647 / -2147483648); and each output width has its own
 *     boundary, so a peak can be clipped for int16 whilst being representable
 *     (and therefore unclipped) for int32.
 *
 *  3. The counter counts the number of *output samples actually limited*,
 *     not the number of input peaks: a band-limited signal can have clipped
 *     output samples even when no input sample reaches full scale.
 *
 *  4. Float output keeps the full amplitude (including the overshoot) and
 *     never contributes to the clip count, whilst integer output adds only
 *     quantisation error to samples that are not clipped; that error is
 *     checked separately from the clipping.
 *
 *  5. The dither switch (int16) moves a sample by less than one LSB: values
 *     must not wrap around, and dithering can only change the clip status of
 *     samples lying within about one LSB of a boundary.  (The dither sequence
 *     is seeded from the clock, so these are properties that hold for any
 *     sequence, rather than exact expected outputs.)  int32 output is
 *     unaffected by the switch.
 *
 *  6. Interleaved, split-channel and per-channel-mono output agree, and the
 *     per-channel clip counts add up.
 *
 *  7. The clip count accumulates over successive soxr_process() calls.
 *
 * All signals used here have magnitude <= 1.  Where a boundary has to be
 * probed to within a fraction of an LSB, output gain is used to place the
 * filtered value at the required offset (the 1:1 rate ratio being an exact
 * gain stage); expected clip counts are otherwise derived from a float64
 * reference run of the *same* resampler configuration, so they do not depend
 * on the filter shape of any particular engine.  Boundary probes sit a
 * quarter of an LSB away from the boundary, so no result depends on a
 * rounding tie. */

#ifdef NDEBUG /* N.B. assert used with active statements so enable always. */
#undef NDEBUG /* Must undef above assert.h or other that might include it. */
#endif

#include <soxr.h>
#include "std-types.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#define NMAX 16000              /* Max. samples per channel in a buffer. */
#define limit16  32768.         /* Value of +1.0, in int16 LSBs. */
#define limit32  2147483648.    /* Value of +1.0, in int32 LSBs. */

#define I16_MAX  32767
#define I16_MIN  (-32767 - 1)
#define I32_MAX  2147483647L
#define I32_MIN  (-2147483647L - 1)

/* A pre-quantisation value d, in LSBs, is clipped iff rounding it would leave
 * the type's range, i.e. iff it is at, or beyond, half an LSB outside a
 * limit: */
#define clips_pos(d, limit)  ((d) >= (limit) - .5)
#define clips_neg(d, limit)  ((d) <= -(limit) - .5)

static int failures;

static void check(int ok, char const * what)
{
  if (!ok) {
    printf("FAIL: %s\n", what);
    ++failures;
  }
}

static void check_count(char const * what, size_t got, size_t want)
{
  if (got != want) {
    printf("FAIL: %s: got %u, expected %u\n",
        what, (unsigned)got, (unsigned)want);
    ++failures;
  }
}

static void check_le(char const * what, double got, double limit)
{
  if (!(got <= limit)) {
    printf("FAIL: %s: got %.17g, limit %.17g\n", what, got, limit);
    ++failures;
  }
}

static void note(char const * what)
{
  printf("  %s\n", what);
}



/* ------------------------------- Test setup ------------------------------- */



typedef struct {
  double irate, orate;
  double gain;                /* As io_spec.scale. */
  soxr_datatype_t otype;
  unsigned long flags;        /* As io_spec.flags (e.g. SOXR_NO_DITHER). */
  unsigned chans;
  unsigned long recipe;       /* Quality recipe. */
  unsigned long q_flags;      /* Quality flags. */
} cfg_t;

static void cfg_init(cfg_t * c, soxr_datatype_t otype, unsigned chans)
{
  memset(c, 0, sizeof(*c));
  c->irate = 44100., c->orate = 96000.;
  c->gain = 1.;
  c->otype = otype;
  c->chans = chans;
  c->recipe = SOXR_HQ;
}

/* Use the double-precision engine: gains and offsets used here are chosen to
 * be exactly representable, and this engine keeps them so (the higher
 * precision path is also what makes sub-LSB int32 probes meaningful).  If the
 * library was built without the double-precision engine, the flag is ignored
 * and the single-precision path is used; see double_precision_engine(). */
static void use_double_precision(cfg_t * c)
{
  c->recipe = SOXR_VHQ;
  c->q_flags = SOXR_DOUBLE_PRECISION;
}

/* Resample `ilen' samples per channel in one (flushing) call.  Returns the
 * number of output samples per channel; `clips', if given, receives the final
 * value of the clip counter. */
static size_t resample(cfg_t const * c, void const * in, size_t ilen,
    soxr_out_t out, size_t olen, size_t * clips)
{
  soxr_error_t error = NULL;
  soxr_io_spec_t io_spec = soxr_io_spec(SOXR_FLOAT64_I, c->otype);
  soxr_quality_spec_t q_spec = soxr_quality_spec(c->recipe, c->q_flags);
  soxr_t p;
  size_t idone = 0, odone = 0;

  io_spec.scale = c->gain;
  io_spec.flags = c->flags;
  p = soxr_create(c->irate, c->orate, c->chans,
      &error, &io_spec, &q_spec, NULL);
  if (!error)
    error = soxr_process(p, in, ~ilen, &idone, out, olen, &odone);
  if (error) {
    printf("FAIL: resample: %s\n", soxr_strerror(error));
    ++failures;
    odone = 0;
  }
  if (clips)
    *clips = p? *soxr_num_clips(p) : 0;
  soxr_delete(p);
  return odone;
}

/* As resample(), but with split (one buffer per channel) output. */
static size_t resample_split(cfg_t const * c, void const * in, size_t ilen,
    soxr_bufs_t out, size_t olen, size_t * clips)
{
  soxr_error_t error = NULL;
  soxr_io_spec_t io_spec = soxr_io_spec(SOXR_FLOAT64_I,
      (soxr_datatype_t)(c->otype | SOXR_SPLIT));
  soxr_quality_spec_t q_spec = soxr_quality_spec(c->recipe, c->q_flags);
  soxr_t p;
  size_t idone = 0, odone = 0;

  io_spec.scale = c->gain;
  io_spec.flags = c->flags;
  p = soxr_create(c->irate, c->orate, c->chans,
      &error, &io_spec, &q_spec, NULL);
  if (!error)
    error = soxr_process(p, in, ~ilen, &idone, (soxr_out_t)out, olen, &odone);
  if (error) {
    printf("FAIL: resample_split: %s\n", soxr_strerror(error));
    ++failures;
    odone = 0;
  }
  if (clips)
    *clips = p? *soxr_num_clips(p) : 0;
  soxr_delete(p);
  return odone;
}

/* Reference run: the same resampling and gain, with float64 output scaled to
 * `full_scale' LSBs per unit.  The integer output conversion reads exactly
 * these (pre-quantisation) values, so the reference states what the integer
 * converter should be seeing. */
static size_t reference(cfg_t const * c, void const * in, size_t ilen,
    double * out, size_t olen, double full_scale)
{
  cfg_t rc = *c;
  rc.otype = SOXR_FLOAT64_I;
  rc.flags = 0;
  rc.gain = c->gain * full_scale;
  return resample(&rc, in, ilen, out, olen, NULL);
}

/* Whether the library actually provides the double-precision engine.  A
 * build without it (WITH_CR64/WITH_CR64S off) uses single precision, whose
 * 24-bit significand cannot represent LSB offsets of a few units at int32
 * magnitude; the sub-LSB int32 boundary probes are therefore only meaningful
 * when this returns true.  The 1:1 ratio is an exact gain stage, so a value
 * needing more than float precision tests the engine directly. */
static int double_precision_engine(void)
{
  static int known = -1;
  if (known < 0) {
    double const in = 1. + ldexp(1., -30); /* Not representable as float. */
    double out = 0;
    cfg_t c;
    cfg_init(&c, SOXR_FLOAT64_I, 1);
    c.orate = c.irate;
    use_double_precision(&c);
    known = resample(&c, &in, 1, &out, 1, NULL) == 1 && out == in;
  }
  return known;
}



/* Deterministic test signals; all have magnitude <= 1. */

static void make_dc(double * buf, size_t n, double level)
{
  size_t i;
  for (i = 0; i < n; ++i)
    buf[i] = level;
}

/* Square wave of the given period, in input samples. */
static void make_square(double * buf, size_t n, double period, double amp)
{
  size_t i;
  for (i = 0; i < n; ++i)
    buf[i] = ((size_t)((double)i / (.5 * period)) & 1)? -amp : amp;
}

static void make_sine(double * buf, size_t n, double freq, double amp)
{
  size_t i;
  for (i = 0; i < n; ++i)
    buf[i] = amp * sin(2 * M_PI * freq * (double)i / 44100.);
}



/* --------------------------- Clip-count prediction ------------------------- */



/* Count the reference samples that must be clipped.  Samples within
 * `epsilon' of a boundary are reported separately in `ambiguous', because
 * their rounding (and whether they just cross the boundary) depends on
 * details -- such as tie-breaking, or the precision of the engine that
 * produced the reference -- that are not part of the clipping contract.
 * The epsilon is chosen far below the quarter-LSB offsets of the boundary
 * probes, yet well above the reference's own rounding error. */
static size_t predicted_clips(double const * ref, size_t n, double limit,
    size_t * ambiguous)
{
  size_t i, clips = 0;
  /* For int16, a float (single-precision) engine's reference value is
   * accurate to about 1e-2 LSB; for int32 it is accurate to about 1e3 LSB
   * (a relative 1e-6 of full scale either way).  int32 boundary probes are
   * only run with the double-precision engine present, hence 1e-3. */
  double epsilon;
  if (limit == limit16)
    epsilon = .05;
  else
    epsilon = double_precision_engine()? 1e-3 : 1e3;
  *ambiguous = 0;
  for (i = 0; i < n; ++i) {
    double const d = ref[i];
    if (d > limit - .5 + epsilon || d < -limit - .5 - epsilon)
      ++clips;
    else if (fabs(d - (limit - .5)) <= epsilon ||
             fabs(d - (-limit - .5)) <= epsilon)
      ++*ambiguous;
  }
  return clips;
}

static void check_predicted(char const * what, size_t got, size_t predicted,
    size_t ambiguous)
{
  if (ambiguous) {
    if (!(got >= predicted && got <= predicted + ambiguous)) {
      printf("FAIL: %s: got %u, predicted %u (+%u ambiguous)\n", what,
          (unsigned)got, (unsigned)predicted, (unsigned)ambiguous);
      ++failures;
    }
  }
  else
    check_count(what, got, predicted);
}

/* Check a run of one integer output width against its reference.  Clipped
 * samples must be at the (asymmetric) saturated extreme; the rest must show
 * no more than half an LSB of quantisation error (plus the reference's own
 * error, which is negligible with the double-precision engine).  Returns the
 * largest such error seen amongst the unclipped samples. */
static double check_integer16(int16_t const * out, double const * ref,
    size_t n)
{
  size_t i;
  double max_err = 0;
  for (i = 0; i < n; ++i) {
    double const d = ref[i];
    if (clips_pos(d, limit16))
      check(out[i] == I16_MAX, "int16: clipped sample is at +32767");
    else if (clips_neg(d, limit16))
      check(out[i] == I16_MIN, "int16: clipped sample is at -32768");
    else {
      double const err = fabs((double)out[i] - d);
      check_le("int16: quantisation error <= half an LSB", err, .5 + 1e-9);
      if (err > max_err)
        max_err = err;
    }
  }
  return max_err;
}

static double check_integer32(int32_t const * out, double const * ref,
    size_t n)
{
  size_t i;
  double max_err = 0;
  double const ref_error = double_precision_engine()? 0. : 1e3;
  for (i = 0; i < n; ++i) {
    double const d = ref[i];
    if (clips_pos(d, limit32))
      check(out[i] == I32_MAX, "int32: clipped sample is at +2147483647");
    else if (clips_neg(d, limit32))
      check(out[i] == I32_MIN, "int32: clipped sample is at -2147483648");
    else {
      double const err = fabs((double)out[i] - d);
      /* A single-precision engine's reference value is only accurate to a
       * few thousand LSB at this magnitude. */
      check_le("int32: quantisation error <= half an LSB", err,
          .5 + ref_error);
      if (err > max_err)
        max_err = err;
    }
  }
  return max_err;
}



/* --------------------------------- Tests ---------------------------------- */



/* The 1:1 rate ratio is an exact gain stage.  That is the prerequisite which
 * lets the boundary tests place filtered values at precise LSB offsets, so
 * verify it here. */
static void test_unity_gain_stage(void)
{
  enum {n = 64};
  static double in[n], out[n];
  static double const vals[] = {0., .125, -.125, .5, -.5, 1., -1., .75, -.75};
  size_t const nvals = sizeof(vals) / sizeof(*vals);
  cfg_t c;
  size_t i, nout;

  for (i = 0; i < n; ++i)
    in[i] = vals[i % nvals];
  cfg_init(&c, SOXR_FLOAT64_I, 1);
  c.orate = c.irate;
  use_double_precision(&c);

  nout = resample(&c, in, n, out, n, NULL);
  check_count("1:1: output length", nout, n);
  for (i = 0; i < nout; ++i)
    check(out[i] == in[i], "1:1: pass-through is exact");

  c.gain = 2.;
  nout = resample(&c, in, n, out, n, NULL);
  check_count("1:1 (gain 2): output length", nout, n);
  for (i = 0; i < nout; ++i)
    check(out[i] == 2. * in[i], "1:1: gain is applied exactly");

  note("1:1 ratio is an exact gain stage (used for boundary probes below)");
}



/* Round away from zero (the ISO C float-to-int conversion rule), expressed
 * for a double; values beyond +/-limit are not passed to this. */
static double round_away(double d)
{
  return d < 0? ceil(d - .5) : floor(d + .5);
}

/* Boundary ladder.  Each case places the (exactly known) filtered value at a
 * quarter LSB offset from a limit, using a one-sample stream.  Offsets that
 * leave the value within the range are merely rounded, and may still land on
 * an extreme; offsets that cross the boundary are clipped.  The output value
 * alone does not reveal which happened, so only the clip counter can
 * distinguish them.
 *
 * Positive side (for int16: max 32767, boundary 32767.5, limit 32768):
 *     limit-1.25 -> 32766.75 -> 32767, not clipped
 *     limit-0.75 -> 32767.25 -> 32767, not clipped
 *     limit-0.25 -> 32767.75 -> 32767, clipped
 *     limit+0.25 -> 32768.25 -> 32767, clipped
 * Negative side (min -32768, boundary -32768.5):
 *     -limit+0.75 -> -32767.25 -> -32767, not clipped
 *     -limit-0.25 -> -32768.25 -> -32768, not clipped (and is the extreme)
 *     -limit-0.75 -> -32768.75 -> -32768, clipped
 *     -limit-1.25 -> -32769.25 -> -32768, clipped   */
static void test_boundary_ladder(double limit, int is32)
{
  enum {n = 8};
  static double const offsets[n] = {
    -1.25, -.75, -.25, .25,    /* From +limit. */
    .75, -.25, -.75, -1.25};   /* From -limit (sign inverted below). */
  static int const clipped[n] = {0, 0, 1, 1, 0, 0, 1, 1};
  double in, ref;
  int16_t out16;
  int32_t out32;
  cfg_t c;
  size_t i, nout, case_clips, ambiguous, predicted;
  size_t total = 0;

  cfg_init(&c, is32? SOXR_INT32_I : SOXR_INT16_I, 1);
  c.orate = c.irate;
  c.flags = SOXR_NO_DITHER;   /* Dither would shift values by up to 1 LSB. */
  use_double_precision(&c);

  for (i = 0; i < n; ++i) {
    double const d = (i < 4? limit : -limit) + offsets[i];
    in = i < 4? 1. : -1.;
    c.gain = fabs(d) / limit;
    nout = reference(&c, &in, 1, &ref, 1, limit);
    check_count("ladder: reference output length", nout, 1);
    nout = resample(&c, &in, 1, is32? (soxr_out_t)&out32 : (soxr_out_t)&out16,
        1, &case_clips);
    check_count("ladder: output length", nout, 1);
    /* The probe is at the intended offset (no filter transient at 1:1, and
     * a DC input gives the same value at every output sample). */
    check(ref == d, "ladder: reference value is at the intended offset");
    predicted = predicted_clips(&ref, nout, limit, &ambiguous);
    check_count("ladder: predicted clip count", predicted, clipped[i]? 1 : 0);
    check_predicted("ladder: clip count", case_clips, predicted, ambiguous);
    total += case_clips;
    /* Positive unclipped offsets still round onto the positive extreme, so
     * the same output value occurs both clipped and unclipped. */
    if (!is32 && i < 4)
      check(out16 == I16_MAX,
          "ladder: unclipped positive value rounds to the extreme");

    /* The output is the (away-from-zero) rounding of the reference value,
     * clamped to the range; the saturated cases are the clipped ones. */
    if (is32)
      check(out32 == (clipped[i]?
            (i < 4? I32_MAX : I32_MIN) : (int32_t)round_away(ref)),
          "ladder: int32 output is the clamped (or rounded) value");
    else
      check(out16 == (clipped[i]?
            (i < 4? I16_MAX : I16_MIN) : (int16_t)round_away(ref)),
          "ladder: int16 output is the clamped (or rounded) value");
  }
  check_count("ladder: total clip count", total, 4);
  note(is32?
      "int32 ladder: 4 of 8 extreme-valued outputs were not clipped" :
      "int16 ladder: 4 of 8 extreme-valued outputs were not clipped");
}



/* Exact full-scale input at unity gain.  +1.0 would have to become +limit,
 * which is outside the range, so it is clipped; -1.0 is exactly -limit, which
 * is inside the range, so it is not clipped.  This is the asymmetry of the
 * two's complement range, at amplitudes that are legal inputs. */
static void test_full_scale_asymmetry(void)
{
  enum {n = 64};
  static double in[n], ref[n];
  static int16_t out16[n];
  static int32_t out32[n];
  cfg_t c;
  size_t i, nref, nout, clips;

  cfg_init(&c, SOXR_INT16_I, 1);
  c.orate = c.irate;
  c.flags = SOXR_NO_DITHER;
  use_double_precision(&c);

  make_dc(in, n, 1.);
  nref = reference(&c, in, n, ref, n, limit16);
  nout = resample(&c, in, n, out16, n, &clips);
  check_count("+1.0: output length", nout, nref);
  check_count("+1.0: int16 clip count (one per output sample)", clips, nout);
  for (i = 0; i < nout; ++i)
    check(out16[i] == I16_MAX, "+1.0: int16 saturates to +32767");

  make_dc(in, n, -1.);
  nref = reference(&c, in, n, ref, n, limit16);
  nout = resample(&c, in, n, out16, n, &clips);
  check_count("-1.0: output length", nout, nref);
  check_count("-1.0: int16 clip count (exact, so none)", clips, 0);
  for (i = 0; i < nout; ++i)
    check(out16[i] == I16_MIN, "-1.0: int16 is exactly -32768");

  cfg_init(&c, SOXR_INT32_I, 1);
  c.orate = c.irate;
  use_double_precision(&c);
  make_dc(in, n, 1.);
  nref = reference(&c, in, n, ref, n, limit32);
  nout = resample(&c, in, n, out32, n, &clips);
  check_count("+1.0: int32 output length", nout, nref);
  check_count("+1.0: int32 clip count (one per output sample)", clips, nout);
  for (i = 0; i < nout; ++i)
    check(out32[i] == I32_MAX, "+1.0: int32 saturates to +2147483647");

  make_dc(in, n, -1.);
  nref = reference(&c, in, n, ref, n, limit32);
  nout = resample(&c, in, n, out32, n, &clips);
  check_count("-1.0: int32 clip count (exact, so none)", clips, 0);
  for (i = 0; i < nout; ++i)
    check(out32[i] == I32_MIN, "-1.0: int32 is exactly -2147483648");

  note("full scale: +1.0 clips, -1.0 does not (asymmetric range)");
}



/* The same filtered value can be clipped for int16 yet be representable, and
 * hence unclipped, for int32: each width obeys its own boundary.  The value
 * chosen is 32767.75 int16 LSBs: a quarter LSB beyond int16's positive
 * boundary, but still below int32's (the int16 LSB is 65536 int32 LSBs, so
 * int32 has plenty of room there). */
static void test_width_dependent_boundary(void)
{
  enum {n = 8};
  static double in[n], ref[NMAX];
  static int16_t out16[n];
  static int32_t out32[n];
  cfg_t c;
  size_t n16, n32, clips16, clips32;
  double const d = limit16 - .25;       /* 32767.75, in int16 LSBs. */
  double const in32 = d * 65536.;       /* The same value, in int32 LSBs. */

  check(in32 < limit32 - .5, "width: test value is within int32's range");
  make_dc(in, n, 1.);

  cfg_init(&c, SOXR_INT16_I, 1);
  c.orate = c.irate;
  c.flags = SOXR_NO_DITHER;
  use_double_precision(&c);
  c.gain = d / limit16;
  n16 = reference(&c, in, n, ref, n, limit16);
  n16 = resample(&c, in, n, out16, n, &clips16);
  check_count("width: int16 output length", n16, n);
  check(ref[0] == d, "width: value is at the intended int16 offset");
  check(clips_pos(d, limit16), "width: value is beyond int16's boundary");
  check_count("width: int16 clip count (every sample clipped)", clips16, n16);
  check(out16[0] == I16_MAX, "width: int16 output is the positive extreme");

  cfg_init(&c, SOXR_INT32_I, 1);
  c.orate = c.irate;
  use_double_precision(&c);
  c.gain = d / limit16;                 /* The same pre-quantisation value. */
  n32 = reference(&c, in, n, ref, n, limit32);
  check(ref[0] == in32, "width: same value, expressed in int32 LSBs");
  n32 = resample(&c, in, n, out32, n, &clips32);
  check_count("width: int32 output length", n32, n);
  check(!clips_pos(in32, limit32) && !clips_neg(in32, limit32),
      "width: value is within int32's range");
  check_count("width: int32 clip count (value is representable)", clips32, 0);
  check(out32[0] == (int32_t)in32,
      "width: int32 output is the exact value, not an extreme");

  note("width: a peak clipped by int16 is representable (unclipped) by int32");
}



/* Filter overshoot, and its interaction with output gain, drive the clip
 * count.  Float output keeps the overshoot and is never counted; integer
 * output saturates, and the count agrees with a prediction made from the
 * float reference. */
static void test_overshoot_and_gain(void)
{
  enum {n = 3000};
  static double in[n], ref[NMAX];
  static int16_t out16[NMAX];
  static int32_t out32[NMAX];
  static double const gains[] = {.25, .5, 1., 2.};
  size_t const ngains = sizeof(gains) / sizeof(*gains);
  cfg_t c;
  size_t i, nref, nout, clips, ambiguous, predicted, at_full_scale = 0;
  size_t previous_clips = 0;
  double peak;

  make_square(in, n, 11.025, 1.);
  for (i = 0; i < n; ++i)
    if (fabs(in[i]) >= 1.)
      ++at_full_scale;

  cfg_init(&c, SOXR_FLOAT64_I, 1);
  use_double_precision(&c);

  /* Float output: the overshoot is retained, and nothing is clipped. */
  nref = resample(&c, in, n, ref, NMAX, &clips);
  peak = 0;
  for (i = 0; i < nref; ++i)
    if (fabs(ref[i]) > peak)
      peak = fabs(ref[i]);
  check(peak > 1.1, "overshoot: the filtered peak exceeds full scale");
  check_count("overshoot: float output is never clipped", clips, 0);
  note("overshoot: input within [-1,1] yields a filtered peak > 1");

  cfg_init(&c, SOXR_INT16_I, 1);
  c.flags = SOXR_NO_DITHER;
  use_double_precision(&c);
  for (i = 0; i < ngains; ++i) {
    c.gain = gains[i];
    nref = reference(&c, in, n, ref, NMAX, limit16);
    nout = resample(&c, in, n, out16, NMAX, &clips);
    predicted = predicted_clips(ref, nref, limit16, &ambiguous);
    check_predicted("gain: int16 clip count matches prediction", clips,
        predicted, ambiguous);
    check_integer16(out16, ref, nout);
    if (gains[i] * peak < 1. - .5 / limit16)
      check_count("gain: no clipping while the peak is within range", clips, 0);
    else {
      check(clips > 0, "gain: overshoot clips the integer output");
      if (previous_clips > 0)
        check(clips > previous_clips, "gain: more gain, more clipped samples");
    }
    previous_clips = clips;
  }
  /* The count is of output samples: at gain 2 it exceeds the number of
   * (full-scale) input samples that produced them. */
  check_count("count: test input is full scale throughout", at_full_scale, n);
  check(previous_clips > at_full_scale,
      "count: clipped output samples exceed the input sample count");

  /* The same run for int32, whose finer boundary gives its own count. */
  cfg_init(&c, SOXR_INT32_I, 1);
  use_double_precision(&c);
  nref = reference(&c, in, n, ref, NMAX, limit32);
  nout = resample(&c, in, n, out32, NMAX, &clips);
  predicted = predicted_clips(ref, nref, limit32, &ambiguous);
  check_predicted("int32: clip count matches prediction", clips, predicted,
      ambiguous);
  check_integer32(out32, ref, nout);
  check(clips > 0, "int32: overshoot clips the integer output");

  note("gain: output gain amplifies the overshoot, raising the clip count");
}



/* A signal that never reaches full scale is still clipped, because the
 * filter overshoots: clipping is a property of the filtered output, not of
 * the input. */
static void test_no_full_scale_input(void)
{
  enum {n = 3000};
  static double in[n], ref[NMAX];
  static int16_t out16[NMAX];
  cfg_t c;
  size_t i, nref, nout, clips, ambiguous, predicted, at_full_scale = 0;
  double peak;

  make_square(in, n, 5.5, .9);
  for (i = 0; i < n; ++i)
    if (fabs(in[i]) >= 1.)
      ++at_full_scale;
  check_count("quiet input: no sample at full scale", at_full_scale, 0);

  cfg_init(&c, SOXR_INT16_I, 1);
  c.flags = SOXR_NO_DITHER;
  use_double_precision(&c);

  nref = reference(&c, in, n, ref, NMAX, limit16);
  nout = resample(&c, in, n, out16, NMAX, &clips);
  peak = 0;
  for (i = 0; i < nref; ++i)
    if (fabs(ref[i]) > peak)
      peak = fabs(ref[i]);
  check(peak > limit16 - .5, "quiet input: filtered peak exceeds the range");
  predicted = predicted_clips(ref, nref, limit16, &ambiguous);
  check_predicted("quiet input: clip count matches prediction", clips,
      predicted, ambiguous);
  check_integer16(out16, ref, nout);
  check(clips > 0, "quiet input: clipping without any full-scale input");
  note("count: clipping follows the filtered output, not the input peaks");
}



/* Float output retains the amplitude; integer output adds only quantisation
 * error to samples it does not clip.  The two effects are checked
 * separately. */
static void test_float_versus_integer(void)
{
  enum {n = 3000};
  static double in[n], reff[NMAX], ref[NMAX];
  static int16_t out16[NMAX];
  static int32_t out32[NMAX];
  cfg_t c;
  size_t i, nf, nref, nout, clips, ambiguous, predicted;
  double peak_f = 0, err16, err32;

  make_square(in, n, 8., 1.);

  cfg_init(&c, SOXR_FLOAT64_I, 1);
  use_double_precision(&c);
  nf = resample(&c, in, n, reff, NMAX, &clips);
  check_count("float: nothing is clipped", clips, 0);
  for (i = 0; i < nf; ++i)
    if (fabs(reff[i]) > peak_f)
      peak_f = fabs(reff[i]);
  check(peak_f > 1.1, "float: the overshoot is retained");
  check(peak_f < 4., "float: the overshoot is finite and modest");

  cfg_init(&c, SOXR_INT16_I, 1);
  c.flags = SOXR_NO_DITHER;
  use_double_precision(&c);
  nref = reference(&c, in, n, ref, NMAX, limit16);
  nout = resample(&c, in, n, out16, NMAX, &clips);
  predicted = predicted_clips(ref, nref, limit16, &ambiguous);
  check_predicted("integer: clip count matches prediction", clips, predicted,
      ambiguous);
  err16 = check_integer16(out16, ref, nout);
  check(err16 <= .5 + 1e-9, "integer: quantisation error is bounded");
  check(err16 > .4, "integer: quantisation error is indeed present");
  note("float keeps the amplitude; integer adds <= 1/2 LSB quantisation");

  cfg_init(&c, SOXR_INT32_I, 1);
  use_double_precision(&c);
  nref = reference(&c, in, n, ref, NMAX, limit32);
  nout = resample(&c, in, n, out32, NMAX, &clips);
  predicted = predicted_clips(ref, nref, limit32, &ambiguous);
  check_predicted("int32: clip count matches prediction", clips, predicted,
      ambiguous);
  err32 = check_integer32(out32, ref, nout);
  check(err32 <= .5 + 1e-9, "int32: quantisation error is bounded");
  /* (A single-precision engine quantises exactly, so no rounding error will
   * be observed amongst its unclipped samples.) */
  if (double_precision_engine())
    check(err32 > .4, "int32: quantisation error is indeed present");
}



/* The dither switch (int16).  Dither adds a value smaller than one LSB before
 * rounding, so each dithered sample must lie within one step of the unclipped
 * quantisation (or exactly on the extreme when the result is clamped), must
 * stay within the type's range, and must not wrap around.  The clip count may
 * change only for samples within one LSB of a boundary.  The dither sequence
 * is seeded from the clock, so these are properties that hold for any
 * sequence.  int32 output has no dither, so its output must be identical with
 * and without the flag. */
static void test_dither(void)
{
  enum {n = 3000};
  static double in[n], ref[NMAX];
  static int16_t plain[NMAX], dithered[NMAX];
  static int32_t plain32[NMAX], flagged32[NMAX];
  cfg_t c;
  size_t i, nref, nplain, ndith, clips_plain, clips_dith, n32, clips32;
  size_t near_boundary = 0, differing = 0, clips;

  make_square(in, n, 11.025, 1.);

  cfg_init(&c, SOXR_INT16_I, 1);
  use_double_precision(&c);
  c.flags = SOXR_NO_DITHER;
  nref = reference(&c, in, n, ref, NMAX, limit16);
  nplain = resample(&c, in, n, plain, NMAX, &clips_plain);
  c.flags = 0;                                  /* Dithering enabled. */
  ndith = resample(&c, in, n, dithered, NMAX, &clips_dith);

  check_count("dither: output length unchanged", ndith, nplain);
  check_count("dither: output length unchanged", nplain, nref);
  for (i = 0; i < ndith; ++i) {
    double const d = ref[i];
    double ideal = round_away(d);
    int const near = fabs(d - (limit16 - .5)) <= 1. ||
                     fabs(d - (-limit16 - .5)) <= 1.;
    /* The value the converter is aiming at: the rounded reference, clamped
     * into the type's range. */
    if (ideal > I16_MAX)
      ideal = I16_MAX;
    if (ideal < I16_MIN)
      ideal = I16_MIN;
    check(dithered[i] >= I16_MIN && dithered[i] <= I16_MAX,
        "dither: output stays within the type's range");
    check(fabs((double)dithered[i] - ideal) <= 1.,
        "dither: within one step of the clamped quantisation");
    /* The undithered output is the reference rounded to the nearest integer,
     * so it lies within half an LSB of it (ties may go either way); clipped
     * samples instead sit at the extreme. */
    if (!clips_pos(d, limit16) && !clips_neg(d, limit16))
      check_le("dither off: within half an LSB of the reference",
          fabs((double)plain[i] - d), .5 + 1e-9);
    if (near)
      ++near_boundary;
    if (plain[i] != dithered[i])
      ++differing;
  }
  check(differing > 0, "dither: samples are in fact perturbed");
  /* Only samples near a boundary can change clip status.  (Samples further
   * than one LSB from a boundary cannot be moved across it by dither.) */
  check(clips_plain == clips_dith ||
      (clips_plain < clips_dith? clips_dith - clips_plain <= near_boundary :
       clips_plain - clips_dith <= near_boundary),
      "dither: clip count changes only close to a boundary");
  /* int32 output is not dithered: the flag must have no effect. */
  cfg_init(&c, SOXR_INT32_I, 1);
  use_double_precision(&c);
  c.flags = SOXR_NO_DITHER;
  n32 = resample(&c, in, n, plain32, NMAX, &clips32);
  c.flags = 0;
  {
    size_t const n32b = resample(&c, in, n, flagged32, NMAX, &clips);
    check_count("dither: int32 output length unchanged by the flag", n32b, n32);
    check_count("dither: int32 clip count unchanged by the flag", clips,
        clips32);
    for (i = 0; i < n32b; ++i)
      check(plain32[i] == flagged32[i], "dither: int32 output is undithered");
  }
  note("dither: perturbation < 1 LSB, no wraparound, clip count stable");
}



/* Interleaved, split and per-channel-mono output must agree, and the
 * per-channel clip counts must add up. */
static void test_multichannel(void)
{
  enum {n = 3000, chans = 2};
  static double ch0[n], ch1[n], in[chans * n];
  static int16_t interleaved[chans * NMAX];
  static int16_t split0[NMAX], split1[NMAX], mono[NMAX];
  soxr_buf_t split[chans];
  cfg_t c;
  size_t i, ch, nout, clips, clips_split, clips_mono, clips_sum = 0, nmono;

  make_square(ch0, n, 11.025, 1.);          /* Overshoots: clips. */
  make_sine(ch1, n, 440., .25);             /* Quiet: never clips. */
  for (i = 0; i < n; ++i) {
    in[chans * i] = ch0[i];
    in[chans * i + 1] = ch1[i];
  }

  cfg_init(&c, SOXR_INT16_I, chans);
  c.flags = SOXR_NO_DITHER;
  use_double_precision(&c);
  /* olen, like ilen, counts samples per channel; the interleaved buffer holds
   * `chans' times as many samples. */
  nout = resample(&c, in, n, interleaved, NMAX, &clips);
  check(nout > 0, "multi: interleaved run produces output");
  check(clips > 0, "multi: the loud channel clips");

  cfg_init(&c, SOXR_INT16_I, chans);
  c.flags = SOXR_NO_DITHER;
  use_double_precision(&c);
  split[0] = split0, split[1] = split1;
  {
    size_t const nsplit = resample_split(&c, in, n, split,
        NMAX, &clips_split);
    check_count("multi: split output length equals interleaved", nsplit, nout);
    check_count("multi: split clip count equals interleaved", clips_split,
        clips);
    for (i = 0; i < nsplit; ++i) {
      check(split0[i] == interleaved[chans * i],
          "multi: split channel 0 equals interleaved");
      check(split1[i] == interleaved[chans * i + 1],
          "multi: split channel 1 equals interleaved");
    }
  }

  /* Per-channel mono runs: each channel's contribution to the clip count. */
  for (ch = 0; ch < chans; ++ch) {
    cfg_t mc;
    cfg_init(&mc, SOXR_INT16_I, 1);
    mc.flags = SOXR_NO_DITHER;
    use_double_precision(&mc);
    nmono = resample(&mc, ch? ch1 : ch0, n, mono, NMAX, &clips_mono);
    check_count("multi: mono output length equals interleaved", nmono, nout);
    for (i = 0; i < nmono; ++i)
      check(mono[i] == interleaved[chans * i + ch],
          "multi: mono channel equals its interleaved channel");
    if (ch == 0)
      check(clips_mono > 0, "multi: channel 0 clips");
    else
      check_count("multi: the quiet channel does not clip", clips_mono, 0);
    clips_sum += clips_mono;
  }
  check_count("multi: channel clip counts add up to the total", clips,
      clips_sum);
  note("multi: interleaved, split and per-channel-mono agree; counts add up");
}



/* The clip count accumulates over successive soxr_process() calls, and a
 * chunked run gives the same result as a single call. */
static void test_streaming(void)
{
  enum {n = 3000, chunk = 377};
  static double in[n], ref[NMAX];
  static int16_t out[NMAX], whole[NMAX];
  cfg_t c;
  soxr_error_t error = NULL;
  soxr_io_spec_t io_spec;
  soxr_quality_spec_t q_spec;
  soxr_t p;
  size_t i, total = 0, odone, idone, nwhole, clips_whole, clips_stream;
  size_t previous = 0;

  make_square(in, n, 11.025, 1.);
  cfg_init(&c, SOXR_INT16_I, 1);
  c.flags = SOXR_NO_DITHER;
  use_double_precision(&c);

  nwhole = resample(&c, in, n, whole, NMAX, &clips_whole);
  check_count("stream: single call output length matches reference",
      nwhole, reference(&c, in, n, ref, NMAX, limit16));
  check(clips_whole > 0, "stream: single call clips");

  io_spec = soxr_io_spec(SOXR_FLOAT64_I, SOXR_INT16_I);
  io_spec.flags = SOXR_NO_DITHER;
  q_spec = soxr_quality_spec(c.recipe, c.q_flags);
  p = soxr_create(c.irate, c.orate, 1, &error, &io_spec, &q_spec, NULL);
  check(!error, "stream: create");
  for (i = 0; i < n && !error; i += chunk) {
    size_t const ilen = n - i < chunk? n - i : chunk;
    error = soxr_process(p, in + i, ilen, &idone, out + total,
        NMAX - total, &odone);
    check(!error, "stream: process");
    check(idone == ilen, "stream: all input consumed");
    total += odone;
    check(*soxr_num_clips(p) >= previous, "stream: clip count is monotonic");
    previous = *soxr_num_clips(p);
  }
  while (!error) {
    error = soxr_process(p, NULL, 0, NULL, out + total, NMAX - total, &odone);
    if (error || !odone)
      break;
    total += odone;
  }
  check(!error, "stream: flush");
  clips_stream = *soxr_num_clips(p);
  soxr_delete(p);

  check_count("stream: total output length", total, nwhole);
  check_count("stream: clip count equals the single-call count", clips_stream,
      clips_whole);
  check(clips_stream > (size_t)0, "stream: the chunked run clips");
  for (i = 0; i < total && i < nwhole; ++i)
    check(out[i] == whole[i], "stream: output equals the single-call output");
  note("stream: clip counts accumulate over process() calls");
}



int main(int argc, char const * argv[])
{
  (void)argc, (void)argv;
  printf("clips: test of clipping statistics for integer output\n");

  test_unity_gain_stage();
  test_boundary_ladder(limit16, 0);
  /* A quarter-LSB offset at int32 magnitude is below the resolution of a
   * single-precision engine, so that ladder needs the double-precision one. */
  if (double_precision_engine())
    test_boundary_ladder(limit32, 1);
  else
    note("int32 ladder: skipped (no double-precision engine in this build)");
  test_full_scale_asymmetry();
  test_width_dependent_boundary();
  test_overshoot_and_gain();
  test_no_full_scale_input();
  test_float_versus_integer();
  test_dither();
  test_multichannel();
  test_streaming();

  if (failures) {
    printf("clips: %d failure(s)\n", failures);
    return 1;
  }
  printf("clips: all tests passed\n");
  return 0;
}
