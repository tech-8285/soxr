/* SoX Resampler Library      Copyright (c) 2007-18 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

/* Test 2: multi-channel resampling contracts.
 *
 * The same input audio is resampled through interleaved and split
 * (deinterleaved) I/O layouts, one-shot and streamed with various input
 * chunkings and output capacities.  Failures are tagged so that the cause
 * can be located to a specific contract:
 *
 *   [layout]    Layout conversion alone must never change the signal: all
 *               four I/O layout combinations must give bit-identical output.
 *   [crosstalk] Energy placed in one channel must never leak into another.
 *   [phase]     Identical channel signals must stay bit-identical (no
 *               cross-channel time skew), including the end-of-stream tail.
 *   [stream]    Streaming must consume the same input lengths (idone) and
 *               reproduce the one-shot reference samples exactly, whatever
 *               the output capacity; the per-call (idone, odone) sequence
 *               must be layout-independent.
 *   [delay]     At end-of-input, the residual filter delay reported by
 *               soxr_delay must equal the remaining playback duration,
 *               which is common to all channels. */

#ifdef NDEBUG /* N.B. assert used with active statements so enable always. */
#undef NDEBUG /* Must undef above assert.h or other that might include it. */
#endif

#include <soxr.h>
#include "../examples/examples-common.h"

#define IRATE 9600.
#define ORATE 4410.
#define ILEN  4000                              /* Samples per channel. */
#define OCAP  ((size_t)(ILEN * ORATE / IRATE + .5) + 128)
#define MAX_CH 8
#define MAX_STEPS (ILEN + OCAP + 64)

static int failures;

#define CHECK(cat, cond, ...) do { \
    if (!(cond)) { ++failures; \
      fprintf(stderr, "FAIL [%s] line %d: ", cat, __LINE__); \
      fprintf(stderr, __VA_ARGS__); \
      fputc('\n', stderr); \
    } \
  } while (0)

/* Deterministic, in-band, channel-specific test signal: */
static float test_signal(size_t i, unsigned ch)
{
  double t = (double)i / IRATE, f = 97. * (ch + 1);
  return (float)(.4 * sin(2 * M_PI * f * t) +
                 .2 * sin(2 * M_PI * (f + 131) * t));
}

static float ibuf    [MAX_CH][ILEN];   /* Input, one array per channel. */
static float ilv_in  [MAX_CH * ILEN];  /* The same, channel-interleaved. */
static float ref_out [MAX_CH * OCAP];  /* One-shot reference, interleaved. */
static float test_out[MAX_CH * OCAP];  /* Interleaved output under test. */
static float split_out[MAX_CH][OCAP];  /* Split output under test. */
static size_t ref_len;                 /* Reference samples per channel. */

static size_t log_i1[MAX_STEPS], log_o1[MAX_STEPS]; /* Streaming step logs */
static size_t log_i2[MAX_STEPS], log_o2[MAX_STEPS]; /* for two layouts.    */

enum {FILL_DISTINCT, FILL_IDENTICAL, FILL_SINGLE};

static void fill_input(unsigned ch, int mode, unsigned single)
{
  unsigned c;
  size_t i;
  for (c = 0; c < ch; ++c)
    for (i = 0; i < ILEN; ++i) {
      float v = mode == FILL_SINGLE && c != single? 0 :
          test_signal(i, mode == FILL_IDENTICAL? 0 : c);
      ibuf[c][i] = v;
      ilv_in[i * ch + c] = v;
    }
}

static void check_eq_ref(char const * cat, float const * got, size_t n,
    unsigned ch, char const * what)
{
  size_t i;
  unsigned c;
  for (i = 0; i < n; ++i)
    for (c = 0; c < ch; ++c)
      if (got[i * ch + c] != ref_out[i * ch + c]) {
        CHECK(cat, 0, "%s: ch=%u sample=%u: got %.9g, reference %.9g",
            what, c, (unsigned)i, got[i * ch + c], ref_out[i * ch + c]);
        return;
      }
}

static void check_split_eq_ref(char const * cat, size_t n, unsigned ch,
    char const * what)
{
  size_t i;
  unsigned c;
  for (i = 0; i < n; ++i)
    for (c = 0; c < ch; ++c)
      if (split_out[c][i] != ref_out[i * ch + c]) {
        CHECK(cat, 0, "%s: ch=%u sample=%u: got %.9g, reference %.9g",
            what, c, (unsigned)i, split_out[c][i], ref_out[i * ch + c]);
        return;
      }
}

static size_t oneshot(unsigned ch, soxr_io_spec_t const * io,
    soxr_in_t in, soxr_out_t out, char const * cat)
{
  size_t idone = 0, odone = 0;
  soxr_quality_spec_t q_spec = soxr_quality_spec(SOXR_HQ, 0);
  soxr_error_t error = soxr_oneshot(IRATE, ORATE, ch,
      in, ILEN, &idone, out, OCAP, &odone, io, &q_spec, NULL);
  CHECK(cat, !error, "soxr_oneshot: %s", soxr_strerror(error));
  CHECK(cat, idone == ILEN, "oneshot consumed %u of %u samples/channel",
      (unsigned)idone, (unsigned)ILEN);
  return odone;
}

/* [layout] All four I/O layout combinations must agree bit-exactly: */
static void test_layouts(unsigned ch)
{
  soxr_io_spec_t io_ii = soxr_io_spec(SOXR_FLOAT32_I, SOXR_FLOAT32_I);
  soxr_io_spec_t io_ss = soxr_io_spec(SOXR_FLOAT32_S, SOXR_FLOAT32_S);
  soxr_io_spec_t io_is = soxr_io_spec(SOXR_FLOAT32_I, SOXR_FLOAT32_S);
  soxr_io_spec_t io_si = soxr_io_spec(SOXR_FLOAT32_S, SOXR_FLOAT32_I);
  void const * ip[MAX_CH];
  void       * op[MAX_CH];
  size_t n;
  unsigned c;

  ref_len = oneshot(ch, &io_ii, ilv_in, ref_out, "layout");

  for (c = 0; c < ch; ++c) ip[c] = ibuf[c], op[c] = split_out[c];

  n = oneshot(ch, &io_ss, ip, op, "layout");
  CHECK("layout", n == ref_len, "split->split: %u samples/channel, "
      "reference %u", (unsigned)n, (unsigned)ref_len);
  check_split_eq_ref("layout", min(n, ref_len), ch, "split->split");

  n = oneshot(ch, &io_is, ilv_in, op, "layout");
  CHECK("layout", n == ref_len, "interleaved->split: %u samples/channel, "
      "reference %u", (unsigned)n, (unsigned)ref_len);
  check_split_eq_ref("layout", min(n, ref_len), ch, "interleaved->split");

  n = oneshot(ch, &io_si, ip, test_out, "layout");
  CHECK("layout", n == ref_len, "split->interleaved: %u samples/channel, "
      "reference %u", (unsigned)n, (unsigned)ref_len);
  check_eq_ref("layout", test_out, min(n, ref_len), ch, "split->interleaved");
}

/* [crosstalk] A signal in one channel must not appear in any other: */
static void test_crosstalk(unsigned ch)
{
  soxr_io_spec_t io_ii = soxr_io_spec(SOXR_FLOAT32_I, SOXR_FLOAT32_I);
  soxr_io_spec_t io_ss = soxr_io_spec(SOXR_FLOAT32_S, SOXR_FLOAT32_S);
  void const * ip[MAX_CH];
  void       * op[MAX_CH];
  unsigned single, c;
  size_t i, n;

  for (c = 0; c < ch; ++c) ip[c] = ibuf[c], op[c] = split_out[c];

  for (single = 0; single < ch; ++single) {
    int layout;
    fill_input(ch, FILL_SINGLE, single);
    for (layout = 0; layout < 2; ++layout) {
      size_t peak_i = 0;
      double peak = 0;
      if (layout) n = oneshot(ch, &io_ss, ip, op, "crosstalk");
      else        n = oneshot(ch, &io_ii, ilv_in, test_out, "crosstalk");
      for (i = 0; i < n; ++i)
        for (c = 0; c < ch; ++c) {
          float v = layout? split_out[c][i] : test_out[i * ch + c];
          if (c == single) {
            if (fabs(v) > peak) peak = fabs(v), peak_i = i;
          } else CHECK("crosstalk", v == 0, "%s: signal from ch=%u leaked "
              "%g into ch=%u at sample=%u", layout? "split" : "interleaved",
              single, v, c, (unsigned)i);
        }
      CHECK("crosstalk", peak > .1, "%s: source ch=%u silent (peak=%g at "
          "sample=%u)", layout? "split" : "interleaved", single, peak,
          (unsigned)peak_i);
    }
  }
}

/* Stream all of the global input through one resampler, in chunks of
 * `ichunk' input samples and with output capacity `olen' samples/channel
 * per call; output is collected into test_out (interleaved) or split_out
 * (split).  Returns the total samples/channel produced.  Verifies the
 * [delay] residual-delay contract at every step. */
static size_t run_stream(unsigned ch, size_t ichunk, size_t olen_param,
    int split, size_t * log_i, size_t * log_o, size_t * log_n)
{
  soxr_io_spec_t io_spec = soxr_io_spec(
      split? SOXR_FLOAT32_S : SOXR_FLOAT32_I,
      split? SOXR_FLOAT32_S : SOXR_FLOAT32_I);
  soxr_quality_spec_t q_spec = soxr_quality_spec(SOXR_HQ, 0);
  soxr_error_t error;
  soxr_t soxr = soxr_create(IRATE, ORATE, ch, &error, &io_spec, &q_spec, NULL);
  size_t isent = 0, ototal = 0, steps = 0;
  size_t bad_step = (size_t)-1, bad_o = 0, bad_rem = 0;
  double delay_at_flush = 0, bad_delay = 0;
  size_t ototal_at_flush = 0, flushed = 0;

  CHECK("stream", soxr != 0, "soxr_create: %s", soxr_strerror(error));
  if (!soxr) return 0;

  for (;;) {
    void const * ip[MAX_CH];
    void       * op[MAX_CH];
    soxr_in_t   in = NULL;
    soxr_out_t  out;
    size_t ilen = 0, idone = 0, odone = 0, olen;
    unsigned c;

    if (isent < ILEN) {                 /* Offer the next input chunk. */
      ilen = min(ichunk, ILEN - isent);
      if (split) {
        for (c = 0; c < ch; ++c) ip[c] = ibuf[c] + isent;
        in = ip;
      } else in = ilv_in + isent * ch;
    } else if (!flushed++) {      /* First end-of-input indication. */
      delay_at_flush = soxr_delay(soxr);
      ototal_at_flush = ototal;
    }

    olen = min(olen_param, OCAP - ototal);
    if (split) {
      for (c = 0; c < ch; ++c) op[c] = split_out[c] + ototal;
      out = op;
    } else out = test_out + ototal * ch;

    error = soxr_process(soxr, in, ilen, &idone, out, olen, &odone);
    CHECK("stream", !error, "soxr_process: %s", soxr_strerror(error));
    if (error) break;

    isent += idone;
    ototal += odone;
    if (log_i) log_i[steps] = idone, log_o[steps] = odone;
    ++steps;

    /* [delay] Once all input has been consumed, the output still to come
     * must equal the filter delay line content (an exact conservation law,
     * as in test 1); mid-feed this holds only to within a sample, since
     * the output count of the unsubmitted input is phase-dependent. */
    {
      double delay = soxr_delay(soxr);
      size_t remaining = (size_t)(delay + .5);
      if (bad_step == (size_t)-1 &&
          (delay < 0 || (isent == ILEN && ototal + remaining != ref_len)))
        bad_step = steps, bad_o = ototal, bad_rem = remaining,
        bad_delay = delay;
    }

    if (!in && !odone) break;                        /* Drained. */
    if (steps >= MAX_STEPS) {
      CHECK("stream", 0, "no termination after %u calls (isent=%u/%u)",
          (unsigned)steps, (unsigned)isent, (unsigned)ILEN);
      break;
    }
  }

  CHECK("delay", bad_step == (size_t)-1, "ch=%u ichunk=%u olen=%u %s: "
      "at step %u, produced %u + remaining %u != reference %u (delay=%g)",
      ch, (unsigned)ichunk, (unsigned)olen_param, split? "split":"interleaved",
      (unsigned)bad_step, (unsigned)bad_o, (unsigned)bad_rem,
      (unsigned)ref_len, bad_delay);

  /* [delay] The residual filter delay at end-of-input must account for
   * exactly the playback duration still to be drained (common to all
   * channels, since odone is a per-channel count). */
  CHECK("delay", ototal - ototal_at_flush == (size_t)(delay_at_flush + .5),
      "ch=%u ichunk=%u olen=%u %s: drained %u samples/channel after "
      "end-of-input but residual delay was %g", ch, (unsigned)ichunk,
      (unsigned)olen_param, split? "split" : "interleaved",
      (unsigned)(ototal - ototal_at_flush), delay_at_flush);

  CHECK("stream", isent == ILEN, "ch=%u: consumed %u of %u samples/channel",
      ch, (unsigned)isent, (unsigned)ILEN);

  if (log_n) *log_n = steps;
  soxr_delete(soxr);
  return ototal;
}

/* [stream] Any chunking/capacity, in either layout, must reproduce the
 * one-shot reference exactly, with a layout-independent call sequence: */
static void test_streaming(unsigned ch)
{
  static size_t const ichunks[] = {977, ILEN};     /* Chunked & one-shot. */
  static size_t const olens[] = {1, 23, 257, OCAP};
  unsigned i, j;
  size_t s;

  for (i = 0; i < AL(ichunks); ++i)
    for (j = 0; j < AL(olens); ++j) {
      size_t n1, n2, o1, o2;
      o1 = run_stream(ch, ichunks[i], olens[j], 0, log_i1, log_o1, &n1);
      CHECK("stream", o1 == ref_len, "ch=%u ichunk=%u olen=%u interleaved: "
          "%u samples/channel, reference %u", ch, (unsigned)ichunks[i],
          (unsigned)olens[j], (unsigned)o1, (unsigned)ref_len);
      check_eq_ref("stream", test_out, min(o1, ref_len), ch,
          "streamed interleaved vs one-shot");

      o2 = run_stream(ch, ichunks[i], olens[j], 1, log_i2, log_o2, &n2);
      CHECK("stream", o2 == ref_len, "ch=%u ichunk=%u olen=%u split: "
          "%u samples/channel, reference %u", ch, (unsigned)ichunks[i],
          (unsigned)olens[j], (unsigned)o2, (unsigned)ref_len);
      check_split_eq_ref("stream", min(o2, ref_len), ch,
          "streamed split vs one-shot");

      /* The streaming contract must not depend on the layout: */
      CHECK("stream", n1 == n2, "ch=%u ichunk=%u olen=%u: %u calls "
          "interleaved vs %u split", ch, (unsigned)ichunks[i],
          (unsigned)olens[j], (unsigned)n1, (unsigned)n2);
      for (s = 0; s < min(n1, n2); ++s)
        CHECK("stream", log_i1[s] == log_i2[s] && log_o1[s] == log_o2[s],
            "ch=%u ichunk=%u olen=%u step %u: interleaved idone=%u odone=%u,"
            " split idone=%u odone=%u", ch, (unsigned)ichunks[i],
            (unsigned)olens[j], (unsigned)s, (unsigned)log_i1[s],
            (unsigned)log_o1[s], (unsigned)log_i2[s], (unsigned)log_o2[s]);
    }
}

/* [phase] Identical channel signals must remain bit-identical per channel
 * (no cross-channel time skew), through a small-capacity stream tail: */
static void test_phase(unsigned ch)
{
  size_t n, i;
  unsigned c;
  fill_input(ch, FILL_IDENTICAL, 0);
  n = run_stream(ch, 977, 13, 1, NULL, NULL, NULL);
  CHECK("phase", n == ref_len, "ch=%u: %u samples/channel, reference %u",
      ch, (unsigned)n, (unsigned)ref_len);
  for (c = 1; c < ch; ++c)
    for (i = 0; i < n; ++i)
      if (split_out[c][i] != split_out[0][i]) {
        CHECK("phase", 0, "ch=%u vs ch=0 differ by %g at sample=%u "
            "(cross-channel time skew)", c, split_out[c][i] - split_out[0][i],
            (unsigned)i);
        return;
      }
}

int main(int argc, char const * arg[])
{
  static unsigned const chans[] = {2, 3, MAX_CH};
  unsigned i;

  for (i = 0; i < AL(chans); ++i) {
    unsigned ch = chans[i];
    fill_input(ch, FILL_DISTINCT, 0);
    test_layouts(ch);       /* Also establishes ref_out / ref_len. */
    test_streaming(ch);
    test_crosstalk(ch);
    test_phase(ch);
  }

  fprintf(stderr, "2-multichannel-layout: %s (%d failure%s)\n",
      failures? "FAIL" : "PASS", failures, failures == 1? "" : "s");
  return (0 * argc * !arg) || failures;
}
