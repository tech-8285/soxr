/* SoX Resampler Library      Copyright (c) 2007-15 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

/* Test 2: the streaming contract of soxr_process around a finite output space.
 *
 * A streaming resampler has to keep three books consistent across arbitrary
 * chunk boundaries:
 *
 *   - consumption : input samples accepted from the application,
 *   - production  : output samples handed back to the application,
 *   - delay       : output-sample equivalents still held by the filter.
 *
 * and it has to distinguish three very different situations:
 *
 *   NOT-READY : no output can be produced right now (e.g. olen == 0) but the
 *               stream is alive and more input will be accepted;
 *   DRAINING  : all input has been consumed (in == NULL announced the end) but
 *               the filter tail is still being emptied, delay > 0;
 *   FINISHED  : the tail is empty, delay == 0 and every further call is a
 *               stable no-op.
 *
 * The test checks that
 *   - a zero-length, non-NULL input call is never mistaken for end-of-input;
 *   - end-of-input is indicated solely by in == NULL;
 *   - effective audio is neither duplicated nor lost: chunked streaming is bit
 *     identical to a soxr_oneshot reference, under many chunk schedules and
 *     both idone-pointer calling conventions;
 *   - consumption, production and soxr_delay satisfy one time-length mapping
 *     after every single call, independently of the chunk boundaries;
 *   - the oneshot (~ilen) and the continuous-streaming return semantics agree;
 *   - soxr_clear resets the engine for a fresh signal without changing the
 *     ownership (contents) of the samples already produced.
 */

#ifdef NDEBUG /* N.B. assert used with active statements so enable always. */
#undef NDEBUG /* Must undef above assert.h or other that might include it. */
#endif

#include <soxr.h>
#include "../examples/examples-common.h"

#define ranqd1(x) ((x) = 1664525 * (x) + 1013904223) /* int32_t x */
#define franqd1(x) (float)(ranqd1(x) * (1. / (65536. * 32768.))) /* [-1,1) */

#define CHECK(cond) do { \
  if (!(cond)) { \
    fprintf(stderr, "CHECK failed at %s:%d: %s\n", __FILE__, __LINE__, #cond); \
    abort(); \
  } \
} while (0)

struct scenario {
  double irate, orate;
  unsigned channels;
  size_t n;                 /* Input length, in samples per channel. */
  char const * name;
};

struct scenario const scenarios[] = {
  { 9600.,  4410.,  1, 300,  "non-integer down-rate"},
  { 4410.,  9600.,  1, 300,  "non-integer up-rate"},
  {48000., 48000.,  1, 257,  "unity rate"},
  {16000., 48000.,  1,  99,  "integer x3 up-rate"},
  {48000., 16000.,  1, 999,  "integer /3 down-rate"},
  {96000., 44100.,  1, 1000, "classic down-rate"},
  {44100., 96000.,  1, 777,  "classic up-rate"},
  { 9600.,  4410.,  2, 333,  "2-channel interleaved"},
  {48000., 16000.,  1,   1,  "single input sample"},
};

/* Chunk schedules: cyclic tables of input chunk and output-cap lengths. */
static size_t const ic_extreme[] = {1};
static size_t const oc_extreme[] = {1};
static size_t const ic_mixed[]   = {3, 7, 13, 29};
static size_t const oc_mixed[]   = {2, 5, 11, 1};
static size_t const ic_bigin[]   = {128};
static size_t const oc_bigin[]   = {1, 3};
static size_t const ic_smallin[] = {1, 2};
static size_t const oc_smallin[] = {64, 128};
static size_t const ic_nospace[] = {5, 11};
static size_t const oc_nospace[] = {0, 1, 7};       /* Interleaved zero space. */
static size_t const ic_pow2[]    = {64, 32};
static size_t const oc_pow2[]    = {32, 64};

struct schedule {
  size_t const * ic; size_t nic;
  size_t const * oc; size_t noc;
  char const * name;
};

static struct schedule const schedules[] = {
  {ic_extreme, 1, oc_extreme, 1, "1-in/1-out"},
  {ic_mixed,   4, oc_mixed,   4, "mixed coprime chunks"},
  {ic_bigin,   1, oc_bigin,   2, "large input, 1..3 output space"},
  {ic_smallin, 2, oc_smallin, 2, "1..2 input, large output space"},
  {ic_nospace, 2, oc_nospace, 3, "interleaved zero output space"},
  {ic_pow2,    2, oc_pow2,    2, "power-of-two blocks"},
};

static void fill_signal(float * buf, size_t frames, unsigned channels, int32_t seed)
{
  size_t i;
  int32_t ran = seed;
  for (i = 0; i < frames * channels; ++i)
    buf[i] = franqd1(ran);
}

static void check_audio_identical(char const * tag,
    float const * got, float const * ref, size_t frames, unsigned channels)
{
  size_t i;
  for (i = 0; i < frames * channels; ++i)
    if (got[i] != ref[i]) {
      fprintf(stderr, "%s: audio diverges at sample %u of %u (got %.9g, ref %.9g)\n",
          tag, (unsigned)i, (unsigned)(frames * channels), got[i], ref[i]);
      abort();
    }
}

/* Streaming-state invariant: delay == consumed / io_ratio - produced.
 * For interpolating (low-latency) engines the idealised RHS can transiently go
 * slightly negative - the filter briefly "owes" a fractional output sample;
 * only the flush phase guarantees a non-negative delay.  The equality itself
 * must hold exactly, regardless of sign. */
static void check_streaming_delay(soxr_t s, double io_ratio,
    size_t consumed, size_t produced)
{
  double expected = (double)consumed / io_ratio - (double)produced;
  double delay = soxr_delay(s);
  CHECK(fabs(delay - expected) <= 1. + 1e-9 * fabs(expected));
}

/* Flushing-state invariant.  At the start of flushing the engine resets its
 * output counter to -round(consumed/io_ratio), so the reported delay is
 * round(consumed/io_ratio) - produced: its rounded value is exactly the number
 * of samples still to be emitted.  Once a flush call has actually been made
 * and the tail is empty, delay is exactly zero. */
static void check_flush_delay(soxr_t s, size_t remaining)
{
  double delay = soxr_delay(s);
  CHECK((size_t)(delay + .5) == remaining);
}

static void check_call_bounds(size_t idone, size_t ilen,
    size_t odone, size_t olen, double io_ratio, int use_idone)
{
  CHECK(idone <= ilen);
  CHECK(odone <= olen);
  if (use_idone)
    /* With an idone pointer, consumption is bounded by the output space. */
    CHECK(idone <= (size_t)ceil((double)olen * io_ratio));
  if (olen == 0)
    CHECK(odone == 0);
  if (use_idone && olen == 0)
    CHECK(idone == 0);
}

/* Reference: everything in one soxr_oneshot call. */
static void run_reference(struct scenario const * sc,
    soxr_quality_spec_t const * q_spec, float const * in,
    float * ref, size_t cap, size_t * ref_n)
{
  size_t idone = 0, odone = 0;
  soxr_error_t error;
  double io_ratio = sc->irate / sc->orate;

  error = soxr_oneshot(sc->irate, sc->orate, sc->channels,
      (void *)in, sc->n, &idone, ref, cap, &odone, NULL, q_spec, NULL);
  CHECK(!error);
  CHECK(idone == sc->n);                       /* All input is consumed once. */
  CHECK(odone <= cap);
  /* Final length mapping: the whole filter delay is returned at the tail. */
  CHECK(odone == (size_t)((double)sc->n / io_ratio + .5));
  *ref_n = odone;
}

/* Continuous streaming over one chunk schedule, on an existing resampler. */
static void run_streaming(struct scenario const * sc,
    soxr_quality_spec_t const * q_spec, struct schedule const * sch,
    int use_idone, float const * in, float const * ref, size_t ref_n,
    float * out, size_t cap, size_t * produced_total)
{
  soxr_error_t error;
  soxr_t soxr;
  double io_ratio = sc->irate / sc->orate;
  size_t consumed = 0, produced = 0;
  size_t ii = 0, oo = 0, calls = 0, idone, odone, ilen, olen, k;
  float non_null_dummy = 0.f;

  soxr = soxr_create(sc->irate, sc->orate, sc->channels,
      &error, NULL, q_spec, NULL);
  CHECK(!error);

  /* Phase 1: feed the stream, honouring the finite output space every call. */
  while (consumed < sc->n) {
    ilen = min(sch->ic[ii % sch->nic], sc->n - consumed);
    olen = sch->oc[oo % sch->noc];
    idone = odone = 0;
    error = soxr_process(soxr,
        (void *)(in + consumed * sc->channels), ilen, use_idone ? &idone : NULL,
        out + produced * sc->channels, olen, &odone);
    CHECK(!error);
    check_call_bounds(idone, ilen, odone, olen, io_ratio, use_idone);
    if (use_idone)
      consumed += idone;
    else
      consumed += ilen;                 /* Without idone, all input is queued. */
    CHECK(consumed <= sc->n);
    produced += odone;
    CHECK(produced <= cap);
    check_streaming_delay(soxr, io_ratio, consumed, produced);
    ++ii; ++oo; ++calls;
    CHECK(calls < sc->n * 200 + 100);  /* A zero-space table must cycle past. */
  }

  /* Input is exhausted but the stream is NOT ended: a zero-length call with a
   * non-NULL input pointer is an ordinary call, not end-of-input. */
  idone = odone = 0;
  error = soxr_process(soxr, &non_null_dummy, 0, &idone,
      out + produced * sc->channels, 16, &odone);
  CHECK(!error);
  CHECK(idone == 0);
  CHECK(odone <= 16);
  produced += odone;
  CHECK(produced <= cap);
  check_streaming_delay(soxr, io_ratio, consumed, produced);

  /* Phase 2: in == NULL announces the end; drain the filter tail. */
  while (produced < ref_n) {
    olen = sch->oc[oo % sch->noc];
    idone = odone = 0;
    error = soxr_process(soxr, NULL, 0, &idone,
        out + produced * sc->channels, olen, &odone);
    CHECK(!error);
    CHECK(idone == 0);
    CHECK(odone <= olen);
    if (olen == 0)
      CHECK(odone == 0);                /* NOT-READY again, but still draining. */
    else
      CHECK(odone > 0);                 /* DRAINING: pending tail flows out. */
    produced += odone;
    check_flush_delay(soxr, ref_n - produced);
    ++oo;
  }

  /* Phase 3: FINISHED. The empty tail is stable; repeated end calls are no-ops. */
  for (k = 0; k < 3; ++k) {
    odone = 0;
    error = soxr_process(soxr, NULL, 0, NULL,
        out + produced * sc->channels, 16, &odone);
    CHECK(!error);
    CHECK(odone == 0);
    check_flush_delay(soxr, 0);
  }

  CHECK(consumed == sc->n);
  CHECK(produced == ref_n);
  check_audio_identical(sch->name, out, ref, ref_n, sc->channels);
  *produced_total = produced;
  soxr_delete(soxr);
}

/* Explicit three-state narrative on one engine: NOT-READY, DRAINING, FINISHED,
 * plus the empty-input-call and output-ownership rules. */
static void test_three_states(struct scenario const * sc,
    soxr_quality_spec_t const * q_spec, float const * in,
    float const * ref, size_t ref_n, float * out, size_t cap)
{
  soxr_error_t error;
  soxr_t soxr;
  double io_ratio = sc->irate / sc->orate;
  double delay0;
  size_t consumed = 0, produced = 0, idone, odone, k;
  float non_null_dummy = 0.f;

  (void)cap;
  soxr = soxr_create(sc->irate, sc->orate, sc->channels,
      &error, NULL, q_spec, NULL);
  CHECK(!error);
  CHECK(soxr_delay(soxr) == 0.);

  /* NOT-READY: zero output space can neither consume nor produce, and does
   * not end the stream. */
  idone = odone = 99;
  error = soxr_process(soxr, in, 5, &idone, out, 0, &odone);
  CHECK(!error);
  CHECK(idone == 0);
  CHECK(odone == 0);
  CHECK(soxr_delay(soxr) == 0.);

  /* A zero-length, non-NULL input call (even with no output buffer) is likewise
   * a no-op and must not be mistaken for end-of-input. */
  idone = odone = 99;
  error = soxr_process(soxr, &non_null_dummy, 0, &idone, NULL, 0, &odone);
  CHECK(!error);
  CHECK(idone == 0);
  CHECK(odone == 0);
  CHECK(soxr_delay(soxr) == 0.);

  /* Prove the stream is alive: feed the whole signal in odd little chunks. */
  while (consumed < sc->n) {
    size_t ilen = min((size_t)7, sc->n - consumed);
    idone = odone = 0;
    error = soxr_process(soxr,
        (void *)(in + consumed * sc->channels), ilen, &idone,
        out + produced * sc->channels, 3, &odone);
    CHECK(!error);
    CHECK(idone <= ilen);
    CHECK(odone <= 3);
    consumed += idone;
    produced += odone;
    check_streaming_delay(soxr, io_ratio, consumed, produced);
  }
  CHECK(consumed == sc->n);

  /* Exhausted, but not ended: one more non-NULL empty call keeps streaming. */
  delay0 = soxr_delay(soxr);
  idone = odone = 0;
  error = soxr_process(soxr, &non_null_dummy, 0, &idone,
      out + produced * sc->channels, 8, &odone);
  CHECK(!error);
  CHECK(idone == 0);
  produced += odone;
  check_streaming_delay(soxr, io_ratio, consumed, produced);
  /* An empty call cannot advance the filter: delay may only stay or shrink by
   * the samples just produced, never grow. */
  CHECK(soxr_delay(soxr) <= delay0 + 1e-12);

  /* DRAINING: NULL announces the end.  The rounded delay is the count of
   * pending tail samples; drain it with one-sample output buffers, tolerating
   * a transient call that emits nothing before the FIFO is primed.  Some
   * configurations (e.g. unity-rate cubic) have no tail at all, in which case
   * the resampler moves straight to FINISHED. */
  {
    size_t idle = 0;
    while ((size_t)(soxr_delay(soxr) + .5) > 0) {
      odone = 99;
      error = soxr_process(soxr, NULL, 0, &idone,
          out + produced * sc->channels, 1, &odone);
      CHECK(!error);
      CHECK(idone == 0);
      CHECK(odone <= 1);
      produced += odone;
      check_flush_delay(soxr, ref_n - produced);
      ++idle;
      CHECK(idle <= ref_n + 2);         /* At most one idle call per sample. */
    }
  }

  /* FINISHED: delay is zero and stays zero under further end calls. */
  for (k = 0; k < 4; ++k) {
    odone = 99;
    error = soxr_process(soxr, NULL, 0, &idone,
        out + produced * sc->channels, 16, &odone);
    CHECK(!error);
    CHECK(idone == 0);
    CHECK(odone == 0);
    check_flush_delay(soxr, 0);
  }
  CHECK(produced == ref_n);
  check_audio_identical("three-state", out, ref, ref_n, sc->channels);
  soxr_delete(soxr);
}

/* Empty (but non-NULL) input calls interspersed in a stream must not change
 * any accounting, and the completed stream must still match the reference. */
static void test_empty_input_ignored(struct scenario const * sc,
    soxr_quality_spec_t const * q_spec, float const * in,
    float const * ref, size_t ref_n, float * out, size_t cap)
{
  soxr_error_t error;
  soxr_t soxr;
  double io_ratio = sc->irate / sc->orate;
  size_t consumed = 0, produced = 0, idone, odone;
  float non_null_dummy = 0.f;

  soxr = soxr_create(sc->irate, sc->orate, sc->channels,
      &error, NULL, q_spec, NULL);
  CHECK(!error);

  while (consumed < sc->n) {
    size_t ilen = min((size_t)5, sc->n - consumed);

    idone = odone = 0;
    error = soxr_process(soxr, &non_null_dummy, 0, &idone,
        out + produced * sc->channels, 4, &odone);
    CHECK(!error);
    CHECK(idone == 0);                  /* Nothing consumed. */
    CHECK(odone <= 4);
    produced += odone;
    check_streaming_delay(soxr, io_ratio, consumed, produced);

    idone = odone = 0;
    error = soxr_process(soxr,
        (void *)(in + consumed * sc->channels), ilen, &idone,
        out + produced * sc->channels, 4, &odone);
    CHECK(!error);
    consumed += idone;
    produced += odone;
    check_streaming_delay(soxr, io_ratio, consumed, produced);
  }

  for (;;) {
    odone = 0;
    error = soxr_process(soxr, NULL, 0, &idone,
        out + produced * sc->channels, 9, &odone);
    CHECK(!error);
    produced += odone;
    check_flush_delay(soxr, ref_n - produced);
    if (!odone)
      break;
  }
  check_flush_delay(soxr, 0);
  CHECK(produced == ref_n);
  check_audio_identical("empty-input", out, ref, ref_n, sc->channels);
  (void)cap;
  soxr_delete(soxr);
}

/* soxr_clear prepares a fresh signal with the same configuration; already
 * produced samples keep their ownership (contents) and are not re-emitted. */
static void test_clear(struct scenario const * sc,
    soxr_quality_spec_t const * q_spec, float const * in,
    float const * ref, size_t ref_n, float * out, size_t cap,
    float * saved, size_t saved_cap)
{
  soxr_error_t error;
  soxr_t soxr;
  double io_ratio = sc->irate / sc->orate;
  size_t consumed = 0, produced = 0, idone, odone, k, saved_n;

  soxr = soxr_create(sc->irate, sc->orate, sc->channels,
      &error, NULL, q_spec, NULL);
  CHECK(!error);

  /* Stream roughly half of the signal, then clear mid-stream. */
  while (consumed < sc->n / 2 + 1 && consumed < sc->n) {
    size_t ilen = min((size_t)11, sc->n - consumed);
    idone = odone = 0;
    error = soxr_process(soxr,
        (void *)(in + consumed * sc->channels), ilen, &idone,
        out + produced * sc->channels, 7, &odone);
    CHECK(!error);
    consumed += idone;
    produced += odone;
    check_streaming_delay(soxr, io_ratio, consumed, produced);
  }
  CHECK(consumed > 0);                     /* n >= 1 and input was accepted. */
  saved_n = produced;
  CHECK(saved_n <= saved_cap);
  memcpy(saved, out, saved_n * sc->channels * sizeof(float));

  error = soxr_clear(soxr);
  CHECK(!error);
  check_flush_delay(soxr, 0);          /* All internal state is discarded. */

  /* The samples already handed to the application are untouched by clear. */
  for (k = 0; k < saved_n * sc->channels; ++k)
    CHECK(out[k] == saved[k]);

  /* A second clear is harmless. */
  CHECK(!soxr_clear(soxr));
  check_flush_delay(soxr, 0);

  /* The same engine now restarts the same signal from sample 0: the result
   * must equal the reference exactly (no duplicated pre-clear audio). */
  consumed = produced = 0;
  while (consumed < sc->n) {
    size_t ilen = min((size_t)13, sc->n - consumed);
    idone = odone = 0;
    error = soxr_process(soxr,
        (void *)(in + consumed * sc->channels), ilen, &idone,
        out + produced * sc->channels, 5, &odone);
    CHECK(!error);
    consumed += idone;
    produced += odone;
    check_streaming_delay(soxr, io_ratio, consumed, produced);
  }
  while (produced < ref_n) {
    odone = 0;
    error = soxr_process(soxr, NULL, 0, &idone,
        out + produced * sc->channels, 6, &odone);
    CHECK(!error);
    produced += odone;
    check_flush_delay(soxr, ref_n - produced);
  }
  check_flush_delay(soxr, 0);
  CHECK(produced == ref_n);
  check_audio_identical("after-clear", out, ref, ref_n, sc->channels);
  soxr_delete(soxr);

  /* Clearing immediately after creation must also leave a usable engine. */
  soxr = soxr_create(sc->irate, sc->orate, sc->channels,
      &error, NULL, q_spec, NULL);
  CHECK(!error);
  CHECK(!soxr_clear(soxr));
  check_flush_delay(soxr, 0);
  consumed = produced = 0;
  while (consumed < sc->n) {
    idone = odone = 0;
    error = soxr_process(soxr,
        (void *)(in + consumed * sc->channels), sc->n - consumed, &idone,
        out + produced * sc->channels, cap - produced, &odone);
    CHECK(!error);
    consumed += idone;
    produced += odone;
    check_streaming_delay(soxr, io_ratio, consumed, produced);
  }
  while (produced < ref_n) {
    odone = 0;
    error = soxr_process(soxr, NULL, 0, &idone,
        out + produced * sc->channels, cap - produced, &odone);
    CHECK(!error);
    produced += odone;
    check_flush_delay(soxr, ref_n - produced);
  }
  CHECK(produced == ref_n);
  check_audio_identical("fresh-clear", out, ref, ref_n, sc->channels);
  soxr_delete(soxr);
}

/* The oneshot return convention (~ilen on soxr_process) must agree with
 * soxr_oneshot and with the chunked streaming results. */
static void test_tilde_ilen(struct scenario const * sc,
    soxr_quality_spec_t const * q_spec, float const * in,
    float const * ref, size_t ref_n, float * out, size_t cap)
{
  soxr_error_t error;
  soxr_t soxr;
  size_t idone = 0, odone = 0;

  soxr = soxr_create(sc->irate, sc->orate, sc->channels,
      &error, NULL, q_spec, NULL);
  CHECK(!error);
  error = soxr_process(soxr, (void *)in, ~sc->n, &idone,
      out, cap, &odone);
  CHECK(!error);
  CHECK(idone == sc->n);
  CHECK(odone == ref_n);
  check_audio_identical("tilde-ilen", out, ref, ref_n, sc->channels);
  check_flush_delay(soxr, 0);

  odone = 99;
  CHECK(!soxr_process(soxr, NULL, 0, &idone, out, 16, &odone));
  CHECK(odone == 0);
  soxr_delete(soxr);
}

int main(int argc, char const * arg[])
{
  static unsigned long const recipes[] = {SOXR_QQ, SOXR_HQ};
  static size_t const frame_cap = 4096;
  size_t const sample_cap = frame_cap * 2;   /* Scenarios use <= 2 channels. */
  float * in = malloc(sample_cap * sizeof(float));
  float * ref = malloc(sample_cap * sizeof(float));
  float * out = malloc(sample_cap * sizeof(float));
  float * saved = malloc(sample_cap * sizeof(float));
  size_t si, qi, ti, ref_n, produced;
  int32_t seed;

  (void)argc; (void)arg;
  CHECK(in && ref && out && saved);

  for (si = 0; si < AL(scenarios); ++si) {
    struct scenario const * sc = &scenarios[si];
    size_t cap = (size_t)((double)sc->n * sc->orate / sc->irate) + 256;
    CHECK(cap <= frame_cap);

    for (qi = 0; qi < AL(recipes); ++qi) {
      soxr_quality_spec_t q_spec = soxr_quality_spec(recipes[qi], 0);
      CHECK(!q_spec.e);
      seed = (int32_t)(si * 1000 + qi * 17 + 1);
      fill_signal(in, sc->n, sc->channels, seed);

      fprintf(stderr, "scenario %u (%s, ch=%u, n=%u, quality=%lu)\n",
          (unsigned)si, sc->name, sc->channels, (unsigned)sc->n, recipes[qi]);

      run_reference(sc, &q_spec, in, ref, cap, &ref_n);
      /* A lone input sample on a heavy down-rate can yield zero output. */
      CHECK(ref_n <= cap);

      /* Every chunk boundary must preserve the same time-length mapping and
       * produce bit-identical audio: nothing duplicated, nothing lost. */
      for (ti = 0; ti < AL(schedules); ++ti) {
        memset(out, 0xA5, sample_cap * sizeof(float));
        run_streaming(sc, &q_spec, &schedules[ti], 1, in, ref, ref_n,
            out, cap, &produced);
        CHECK(produced == ref_n);
      }

      /* The NULL-idone calling convention queues all offered input; the
       * resulting audio must still match exactly. */
      memset(out, 0x5A, sample_cap * sizeof(float));
      run_streaming(sc, &q_spec, &schedules[1], 0, in, ref, ref_n,
          out, cap, &produced);
      CHECK(produced == ref_n);
      memset(out, 0x3C, sample_cap * sizeof(float));
      run_streaming(sc, &q_spec, &schedules[4], 0, in, ref, ref_n,
          out, cap, &produced);
      CHECK(produced == ref_n);

      test_three_states(sc, &q_spec, in, ref, ref_n, out, cap);
      test_empty_input_ignored(sc, &q_spec, in, ref, ref_n, out, cap);
      test_clear(sc, &q_spec, in, ref, ref_n, out, cap, saved, sample_cap);
      test_tilde_ilen(sc, &q_spec, in, ref, ref_n, out, cap);
    }
  }

  free(in); free(ref); free(out); free(saved);
  fprintf(stderr, "2-stream-contract: all streaming contract checks passed\n");
  return 0;
}
