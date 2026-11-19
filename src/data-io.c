/* SoX Resampler Library      Copyright (c) 2007-16 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

#include <limits.h>
#include <math.h>
#include <string.h>

#include "data-io.h"
#include "internal.h"



/* Packed, 3-byte-per-sample 24-bit PCM in two's-complement (the same signed
 * convention as the 16- and 32-bit integer types), in the host's byte order.
 * LOAD_S24 and STORE_S24 are exact inverses over [-2^23, 2^23-1]. */

#if HAVE_BIGENDIAN
  #define LOAD_S24(p) ((long)(((unsigned long)((p)[0] & 0xff) << 16 | \
      (unsigned long)((p)[1] & 0xff) << 8 | \
      (unsigned long)((p)[2] & 0xff)) ^ 0x800000UL) - 0x800000L)
  #define STORE_S24(p, x) do { \
    unsigned long u = (unsigned long)(int32_t)(x); \
    (p)[0] = (unsigned char)(u >> 16); \
    (p)[1] = (unsigned char)(u >> 8); \
    (p)[2] = (unsigned char)(u); \
  } while (0)
#else
  #define LOAD_S24(p) ((long)(((unsigned long)((p)[0] & 0xff) | \
      (unsigned long)((p)[1] & 0xff) << 8 | \
      (unsigned long)((p)[2] & 0xff) << 16) ^ 0x800000UL) - 0x800000L)
  #define STORE_S24(p, x) do { \
    unsigned long u = (unsigned long)(int32_t)(x); \
    (p)[0] = (unsigned char)(u); \
    (p)[1] = (unsigned char)(u >> 8); \
    (p)[2] = (unsigned char)(u >> 16); \
  } while (0)
#endif



#define DEINTERLEAVE_FROM(T,flag) do { \
  unsigned i; \
  size_t j; \
  T const * src = *src0; \
  if (ch > 1) for (j = 0; j < n; ++j) \
    for (i = 0; i < ch; ++i) dest[i][j] = (DEINTERLEAVE_TO)*src++; \
  else if (flag) memcpy(dest[0], src, n * sizeof(T)), src = &src[n]; \
  else for (j = 0; j < n; dest[0][j++] = (DEINTERLEAVE_TO)*src++); \
  *src0 = src; \
} while (0)

#define DEINTERLEAVE_FROM_24 do { \
  unsigned i; \
  size_t j; \
  unsigned char const * src = *src0; \
  if (ch > 1) for (j = 0; j < n; ++j) \
    for (i = 0; i < ch; ++i) \
      dest[i][j] = (DEINTERLEAVE_TO)LOAD_S24(src), src += 3; \
  else for (j = 0; j < n; ++j) \
    dest[0][j] = (DEINTERLEAVE_TO)LOAD_S24(src), src += 3; \
  *src0 = src; \
} while (0)



#if WITH_CR64 || WITH_CR64S
void _soxr_deinterleave(double * * dest, /* Round/clipping not needed here */
    soxr_datatype_t data_type, void const * * src0, size_t n, unsigned ch)
{
#define DEINTERLEAVE_TO double
  switch (data_type & ~SOXR_SPLIT) {
    case SOXR_FLOAT32: DEINTERLEAVE_FROM(float, 0); break;
    case SOXR_FLOAT64: DEINTERLEAVE_FROM(double, 1); break;
    case SOXR_INT32:   DEINTERLEAVE_FROM(int32_t, 0); break;
    case SOXR_INT16:   DEINTERLEAVE_FROM(int16_t, 0); break;
    case SOXR_INT24:   DEINTERLEAVE_FROM_24; break;
    default: break;
  }
}
#endif



#if WITH_CR32 || WITH_CR32S || WITH_VR32
void _soxr_deinterleave_f(float * * dest, /* Round/clipping not needed here */
    soxr_datatype_t data_type, void const * * src0, size_t n, unsigned ch)
{
#undef DEINTERLEAVE_TO
#define DEINTERLEAVE_TO float
  switch (data_type & ~SOXR_SPLIT) {
    case SOXR_FLOAT32: DEINTERLEAVE_FROM(float, 1); break;
    case SOXR_FLOAT64: DEINTERLEAVE_FROM(double, 0); break;
    case SOXR_INT32:   DEINTERLEAVE_FROM(int32_t, 0); break;
    case SOXR_INT16:   DEINTERLEAVE_FROM(int16_t, 0); break;
    case SOXR_INT24:   DEINTERLEAVE_FROM_24; break;
    default: break;
  }
}
#endif



#include "rint.h"



#if defined FE_INVALID && defined FPU_RINT32 && defined __STDC_VERSION__
  #if __STDC_VERSION__ >= 199901L
    #pragma STDC FENV_ACCESS ON
  #endif
#endif

#if WITH_CR64 || WITH_CR64S
#define FLOATX double

#define LSX_RINT_CLIP_2 lsx_rint32_clip_2
#define LSX_RINT_CLIP lsx_rint32_clip
#define RINT_CLIP rint32_clip
#define RINT rint32D
#if defined FPU_RINT32
  #define FPU_RINT
#endif
#define RINT_T int32_t
#define RINT_MAX 2147483647L
#include "rint-clip.h"

#define LSX_RINT_CLIP_2 lsx_rint16_clip_2
#define LSX_RINT_CLIP lsx_rint16_clip
#define RINT_CLIP rint16_clip
#define RINT rint16D
#if defined FPU_RINT16
  #define FPU_RINT
#endif
#define RINT_T int16_t
#define RINT_MAX 32767
#include "rint-clip.h"

#define LSX_RINT_CLIP_2 lsx_rint16_clip_2_dither
#define LSX_RINT_CLIP lsx_rint16_clip_dither
#define RINT_CLIP rint16_clip_dither
#define RINT rint16D
#if defined FPU_RINT16
  #define FPU_RINT
#endif
#define RINT_T int16_t
#define RINT_MAX 32767
#define DITHER
#include "rint-clip.h"

#undef FLOATX
#endif



#if WITH_CR32 || WITH_CR32S || WITH_VR32
#define FLOATX float

#define LSX_RINT_CLIP_2 lsx_rint32_clip_2_f
#define LSX_RINT_CLIP lsx_rint32_clip_f
#define RINT_CLIP rint32_clip_f
#define RINT rint32F
#if defined FPU_RINT32
  #define FPU_RINT
#endif
#define RINT_T int32_t
#define RINT_MAX 2147483647L
#include "rint-clip.h"

#define LSX_RINT_CLIP_2 lsx_rint16_clip_2_f
#define LSX_RINT_CLIP lsx_rint16_clip_f
#define RINT_CLIP rint16_clip_f
#define RINT rint16F
#if defined FPU_RINT16
  #define FPU_RINT
#endif
#define RINT_T int16_t
#define RINT_MAX 32767
#include "rint-clip.h"

#define LSX_RINT_CLIP_2 lsx_rint16_clip_2_dither_f
#define LSX_RINT_CLIP lsx_rint16_clip_dither_f
#define RINT_CLIP rint16_clip_dither_f
#define RINT rint16D
#if defined FPU_RINT16
  #define FPU_RINT
#endif
#define RINT_T int16_t
#define RINT_MAX 32767
#define DITHER
#include "rint-clip.h"

#undef FLOATX
#endif

#if defined FE_INVALID && defined FPU_RINT32 && defined __STDC_VERSION__
  #if __STDC_VERSION__ >= 199901L
    #pragma STDC FENV_ACCESS OFF
  #endif
#endif



/* Packed 24-bit output: round to nearest (ties away from zero), saturate to
 * [-2^23, 2^23-1] and count saturated samples as clips.  No dithering.  The
 * same portable arithmetic is used for both the single- and double-precision
 * engines so that either gives identical 24-bit results.                  */

#define RINT24_CLIP_FUNCS(SUFFIX, FLOATX_T) \
static size_t lsx_rint24_clip##SUFFIX(void * * const dest0, \
    FLOATX_T const * const src, size_t const n) \
{ \
  static const double N = 8388608.; /* 2^23 */ \
  size_t i, clips = 0; \
  unsigned char * dest = *dest0; \
  double d; \
  for (i = 0; i < n; ++i) { \
    d = (double)src[i]; \
    d = d > 0? d + .5 >= N? ++clips, N - 1 : d + .5 \
             : d - .5 <= -N - 1? ++clips, -N : d - .5; \
    STORE_S24(dest, (long)d); \
    dest += 3; \
  } \
  *dest0 = dest; \
  return clips; \
} \
\
static size_t lsx_rint24_clip_2##SUFFIX(void * * dest0, \
    FLOATX_T const * const * srcs, unsigned const stride, size_t const n) \
{ \
  static const double N = 8388608.; /* 2^23 */ \
  unsigned j; \
  size_t i, clips = 0; \
  unsigned char * dest = *dest0; \
  double d; \
  if (n) { \
    for (j = 0; j < stride; ++j) { \
      FLOATX_T const * const src = srcs[j]; \
      for (i = 0; i < n; ++i) { \
        d = (double)src[i]; \
        d = d > 0? d + .5 >= N? ++clips, N - 1 : d + .5 \
                 : d - .5 <= -N - 1? ++clips, -N : d - .5; \
        STORE_S24(dest + 3 * (stride * i + j), (long)d); \
      } \
    } \
    *dest0 = dest + 3 * stride * n; \
  } \
  return clips; \
}

#if WITH_CR64 || WITH_CR64S
RINT24_CLIP_FUNCS(_d, double)
#endif
#if WITH_CR32 || WITH_CR32S || WITH_VR32
RINT24_CLIP_FUNCS(_f, float)
#endif

#undef RINT24_CLIP_FUNCS



#define INTERLEAVE_TO(T,flag) do { \
  unsigned i; \
  size_t j; \
  T * dest = *dest0; \
  if (ch > 1) \
  for (j = 0; j < n; ++j) for (i = 0; i < ch; ++i) *dest++ = (T)src[i][j]; \
  else if (flag) memcpy(dest, src[0], n * sizeof(T)), dest = &dest[n]; \
  else for (j = 0; j < n; *dest++ = (T)src[0][j++]); \
  *dest0 = dest; \
  return 0; \
} while (0)

#if WITH_CR64 || WITH_CR64S
size_t /* clips */ _soxr_interleave(soxr_datatype_t data_type, void * * dest0,
  double const * const * src, size_t n, unsigned ch, unsigned long * seed)
{
  switch (data_type & ~SOXR_SPLIT) {
    case SOXR_FLOAT32: INTERLEAVE_TO(float, 0);
    case SOXR_FLOAT64: INTERLEAVE_TO(double, 1);

    case SOXR_INT32: if (ch == 1)
        return lsx_rint32_clip(dest0, src[0], n);
      return lsx_rint32_clip_2(dest0, src, ch, n);

    case SOXR_INT16: if (seed) {
      if (ch == 1)
        return lsx_rint16_clip_dither(dest0, src[0], n, seed);
      return lsx_rint16_clip_2_dither(dest0, src, ch, n, seed);
    }
    if (ch == 1)
        return lsx_rint16_clip(dest0, src[0], n);
      return lsx_rint16_clip_2(dest0, src, ch, n);

    case SOXR_INT24: if (ch == 1)
        return lsx_rint24_clip_d(dest0, src[0], n);
      return lsx_rint24_clip_2_d(dest0, src, ch, n);
    default: break;
  }
  return 0;
}
#endif

#if WITH_CR32 || WITH_CR32S || WITH_VR32
size_t /* clips */ _soxr_interleave_f(soxr_datatype_t data_type, void * * dest0,
  float const * const * src, size_t n, unsigned ch, unsigned long * seed)
{
  switch (data_type & ~SOXR_SPLIT) {
    case SOXR_FLOAT32: INTERLEAVE_TO(float, 1);
    case SOXR_FLOAT64: INTERLEAVE_TO(double, 0);

    case SOXR_INT32: if (ch == 1)
        return lsx_rint32_clip_f(dest0, src[0], n);
      return lsx_rint32_clip_2_f(dest0, src, ch, n);

    case SOXR_INT16: if (seed) {
      if (ch == 1)
        return lsx_rint16_clip_dither_f(dest0, src[0], n, seed);
      return lsx_rint16_clip_2_dither_f(dest0, src, ch, n, seed);
    }
    if (ch == 1)
        return lsx_rint16_clip_f(dest0, src[0], n);
      return lsx_rint16_clip_2_f(dest0, src, ch, n);

    case SOXR_INT24: if (ch == 1)
        return lsx_rint24_clip_f(dest0, src[0], n);
      return lsx_rint24_clip_2_f(dest0, src, ch, n);
    default: break;
  }
  return 0;
}
#endif
