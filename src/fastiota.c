#include "fastiota.h"
#ifndef FASTIOTA_C_H__
#define FASTIOTA_C_H__
#include <x86intrin.h>
#include <assert.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef signed long long signed_size_t;

#if __AVX512F__
#define init_start(func, sv)   func(sv + 15, sv + 14, sv + 13, sv + 12, sv + 11, sv + 10, sv + 9, sv + 8, sv + 7, sv + 6, sv + 5, sv + 4, sv + 3, sv + 2, sv + 1, sv)
#define init_start64(func, sv) func(sv + 7, sv + 6, sv + 5, sv + 4, sv + 3, sv + 2, sv + 1, sv)
#define init_start_val(func, sv, val)  func(sv + val * 15, sv + val * 14, sv + val * 13, sv + val * 12, sv + val * 11, sv + val * 10, sv + val * 9, sv + val * 8, sv + val * 7, sv + val * 6, sv + val * 5, sv + val * 4, sv + val * 3, sv + val * 2, sv + val * 1, sv)
#define init_start_val64(func, sv, val)  func(sv + val * 7, sv + val * 6, sv + val * 5, sv + val * 4, sv + val * 3, sv + val * 2, sv + val * 1, sv)
#define VECSIZE sizeof(__m512i)
#elif __AVX2__
#define init_start(func, sv)   func(sv + 7, sv + 6, sv + 5, sv + 4, sv + 3, sv + 2, sv + 1, sv)
#define init_start64(func, sv) func(sv + 3, sv + 2, sv + 1, sv)
#define init_start_val(func, sv, val)  func(sv + val * 7, sv + val * 6, sv + val * 5, sv + val * 4, sv + val * 3, sv + val * 2, sv + val * 1, sv)
#define init_start_val64(func, sv, val)  func(sv + val * 3, sv + val * 2, sv + val * 1, sv)
#define VECSIZE sizeof(__m256i)
#elif __SSE2__
#define init_start(func, sv)   func(sv + 3, sv + 2, sv + 1, sv)
#define init_start64(func, sv) func(sv + 1, sv)
#define init_start_val(func, sv, val)  func(sv + val * 3, sv + val * 2, sv + val * 1, sv)
#define init_start_val64(func, sv, val)  func(sv + val * 1, sv)
#define VECSIZE sizeof(__m128i)
#endif

#define UNROLL4 _Pragma("GCC unroll 4")

FASTIOTA_EXPORT uint32_t *fastiota32(uint32_t *start, size_t nelem, uint32_t sv) {
#define CASEN *start++ = sv++; __attribute__((fallthrough));
    {
        const size_t leftover = VECSIZE / sizeof(uint32_t) - (uint64_t)start % (VECSIZE / sizeof(uint32_t));
        if(leftover && leftover != VECSIZE / sizeof(uint32_t)) {
            nelem -= leftover;
            switch(leftover) {
                case 15: CASEN case 14: CASEN case 13: CASEN case 12: CASEN case 11: CASEN case 10: CASEN case 9: CASEN case 8: CASEN case 7: CASEN case 6: CASEN case 5: CASEN case 4: CASEN case 3: CASEN case 2: CASEN case 1: CASEN
                default: ;
            }
        }
    }
    uint32_t *const end = start + nelem, *const origstart = start;
#if __AVX512F__
    __m512i m = init_start(_mm512_set_epi32, sv);
    const __m512i inc = _mm512_set1_epi32(sizeof(m) / sizeof(uint32_t));
    assert((uint64_t)start % sizeof(m) == 0);
    if((uint64_t)start % sizeof(m) == 0) {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint32_t))) {
            _mm512_store_si512((void *)start, m);
            m = _mm512_add_epi32(m, inc);
            start += sizeof(m) / sizeof(uint32_t);
        }
    } else {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint32_t))) {
            _mm512_storeu_si512((void *)start, m);
            m = _mm512_add_epi32(m, inc);
            start += sizeof(m) / sizeof(uint32_t);
        }
    }
    sv += (start - origstart);
#elif __AVX2__
    __m256i m = init_start(_mm256_set_epi32, sv);
    assert((uint64_t)start % sizeof(m) == 0);
    const __m256i inc = _mm256_set1_epi32(sizeof(m) / sizeof(uint32_t));
    if((uint64_t)start % sizeof(m) == 0) {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint32_t))) {
            _mm256_store_si256((__m256i *)start, m);
            m = _mm256_add_epi32(m, inc);
            start += sizeof(m) / sizeof(uint32_t);
        }
    } else {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint32_t))) {
            _mm256_storeu_si256((__m256i *)start, m);
            m = _mm256_add_epi32(m, inc);
            start += sizeof(m) / sizeof(uint32_t);
        }
    }
    sv += (start - origstart);
#elif __SSE2__
    __m128i m = init_start(_mm_set_epi32, sv);
    assert((uint64_t)start % sizeof(m) == 0);
    const __m128i inc = _mm_set1_epi32(sizeof(m) / sizeof(uint32_t));
    if((uint64_t)start % sizeof(m) == 0) {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint32_t))) {
            _mm_store_si128((__m128i *)start, m);
            m = _mm_add_epi32(m, inc);
            start += sizeof(m) / sizeof(uint32_t);
        }
    } else {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint32_t))) {
            _mm_storeu_si128((__m128i *)start, m);
            m = _mm_add_epi32(m, inc);
            start += sizeof(m) / sizeof(uint32_t);
        }
    }
    sv += (start - origstart);
#endif
    while(start < end) *start++ = sv++;
    return origstart;
}

FASTIOTA_EXPORT int32_t *fastiota32i(int32_t *start, size_t nelem, int32_t sv)
{
    return (int32_t *)fastiota32((uint32_t *)start, nelem, sv);
}
FASTIOTA_EXPORT uint64_t *fastiota64(uint64_t *start, size_t nelem, uint64_t sv) {
    {
        const size_t leftover = VECSIZE / sizeof(uint64_t) - (uint64_t)start % (VECSIZE / sizeof(uint64_t));
        if(leftover && leftover != VECSIZE / sizeof(uint64_t)) {
            nelem -= leftover;
            switch(leftover) {
                case 7: CASEN case 6: CASEN case 5: CASEN case 4: CASEN case 3: CASEN case 2: CASEN case 1: CASEN
                default: ;
            }
        }
    }
#undef CASEN
    uint64_t *const end = start + nelem, *const origstart = start;
#if __AVX512F__
    __m512i m = init_start64(_mm512_set_epi64, sv);
    const __m512i inc = _mm512_set1_epi64((signed_size_t)(sizeof(m) / sizeof(uint64_t)));
    if((uint64_t)start % sizeof(m) == 0) {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint64_t))) {
            _mm512_store_si512((void *)start, m);
            m = _mm512_add_epi64(m, inc);
            start += (signed_size_t)(sizeof(m) / sizeof(uint64_t));
        }
    } else {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint64_t))) {
            _mm512_storeu_si512((void *)start, m);
            m = _mm512_add_epi64(m, inc);
            start += (signed_size_t)(sizeof(m) / sizeof(uint64_t));
        }
    }
    sv += (start - origstart);
#elif __AVX2__
    __m256i m = init_start64(_mm256_set_epi64x, sv);
    const __m256i inc = _mm256_set1_epi64x((signed_size_t)(sizeof(m) / sizeof(uint64_t)));
    if((uint64_t)start % sizeof(m) == 0) {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint64_t))) {
            _mm256_store_si256((__m256i *)start, m);
            m = _mm256_add_epi64(m, inc);
            start += (signed_size_t)(sizeof(m) / sizeof(uint64_t));
        }
    } else {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint64_t))) {
            _mm256_storeu_si256((__m256i *)start, m);
            m = _mm256_add_epi64(m, inc);
            start += (signed_size_t)(sizeof(m) / sizeof(uint64_t));
        }
    }
    sv += (start - origstart);
#elif __SSE2__
    __m128i m = init_start64(_mm_set_epi64x, sv);
    const __m128i inc = _mm_set1_epi64x((signed_size_t)(sizeof(m) / sizeof(uint64_t)));
    if((uint64_t)start % sizeof(m) == 0) {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint64_t))) {
            _mm_store_si128((__m128i *)start, m);
            m = _mm_add_epi64(m, inc);
            start += (signed_size_t)(sizeof(m) / sizeof(uint64_t));
        }
    } else {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint64_t))) {
            _mm_storeu_si128((__m128i *)start, m);
            m = _mm_add_epi64(m, inc);
            start += (signed_size_t)(sizeof(m) / sizeof(uint64_t));
        }
    }
    sv += (start - origstart);
#endif
    while(start < end) *start++ = sv++;
    return origstart;
}

FASTIOTA_EXPORT int64_t *fastiota64i(int64_t *start, size_t nelem, int64_t sv)
{
    return (int64_t *)fastiota64((uint64_t *)start, nelem, sv);
}

FASTIOTA_EXPORT uint32_t *fastiota32_inc(uint32_t *start, size_t nelem, uint32_t sv, uint32_t incv) {
    if(incv == 1) return fastiota32(start, nelem, sv);
    uint32_t *const end = start + nelem, *const origstart = start;
#if __AVX512F__
    __m512i m = init_start_val(_mm512_set_epi32, sv, incv);
    const __m512i inc = _mm512_set1_epi32(sizeof(m) / sizeof(uint32_t));
    if((uint64_t)start % sizeof(m) == 0) {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint32_t))) {
            _mm512_store_si512((void *)start, m);
            m = _mm512_add_epi32(m, inc);
            start += sizeof(m) / sizeof(uint32_t);
        }
    } else {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint32_t))) {
            _mm512_storeu_si512((void *)start, m);
            m = _mm512_add_epi32(m, inc);
            start += sizeof(m) / sizeof(uint32_t);
        }
    }
    sv += (start - origstart);
#elif __AVX2__
    __m256i m = init_start_val(_mm256_set_epi32, sv, incv);
    const __m256i inc = _mm256_set1_epi32((signed_size_t)(sizeof(m) / sizeof(uint32_t)));
    if((uint64_t)start % sizeof(m) == 0) {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint32_t))) {
            _mm256_store_si256((__m256i *)start, m);
            m = _mm256_add_epi32(m, inc);
            start += sizeof(m) / sizeof(uint32_t);
        }
    } else {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint32_t))) {
            _mm256_storeu_si256((__m256i *)start, m);
            m = _mm256_add_epi32(m, inc);
            start += sizeof(m) / sizeof(uint32_t);
        }
    }
    sv += (start - origstart);
#elif __SSE2__
    __m128i m = init_start_val(_mm_set_epi32, sv, incv);
    const __m128i inc = _mm_set1_epi32(sizeof(m) / sizeof(uint32_t));
    if((uint64_t)start % sizeof(m) == 0) {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint32_t))) {
            _mm_store_si128((__m128i *)start, m);
            m = _mm_add_epi32(m, inc);
            start += sizeof(m) / sizeof(uint32_t);
        }
    } else {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint32_t))) {
            _mm_storeu_si128((__m128i *)start, m);
            m = _mm_add_epi32(m, inc);
            start += sizeof(m) / sizeof(uint32_t);
        }
    }
    sv += (start - origstart);
#endif
    while(start < end) *start++ = sv++;
    return origstart;
}

FASTIOTA_EXPORT uint64_t *fastiota64_inc(uint64_t *start, size_t nelem, uint64_t sv, uint64_t incv) {
    if(incv == 1) return fastiota64(start, nelem, sv);
    uint64_t *const end = start + nelem, *const origstart = start;
#if __AVX512F__
    __m512i m = init_start_val64(_mm512_set_epi64, sv, incv);
    const __m512i inc = _mm512_set1_epi64((signed_size_t)(sizeof(m) / sizeof(uint64_t)));
    if((uint64_t)start % sizeof(m) == 0) {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint64_t))) {
            _mm512_store_si512((void *)start, m);
            m = _mm512_add_epi64(m, inc);
            start += (signed_size_t)(sizeof(m) / sizeof(uint64_t));
        }
    } else {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint64_t))) {
            _mm512_storeu_si512((void *)start, m);
            m = _mm512_add_epi64(m, inc);
            start += (signed_size_t)(sizeof(m) / sizeof(uint64_t));
        }
    }
    sv += (start - origstart);
#elif __AVX2__
    __m256i m = init_start_val64(_mm256_set_epi64x, sv, incv);
    const __m256i inc = _mm256_set1_epi64x((signed_size_t)(sizeof(m) / sizeof(uint64_t)));
    if((uint64_t)start % sizeof(m) == 0) {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint64_t))) {
            _mm256_store_si256((__m256i *)start, m);
            m = _mm256_add_epi64(m, inc);
            start += (signed_size_t)(sizeof(m) / sizeof(uint64_t));
        }
    } else {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint64_t))) {
            _mm256_storeu_si256((__m256i *)start, m);
            m = _mm256_add_epi64(m, inc);
            start += (signed_size_t)(sizeof(m) / sizeof(uint64_t));
        }
    }
    sv += (start - origstart);
#elif __SSE2__
    __m128i m = init_start_val64(_mm_set_epi64x, sv, incv);
    const __m128i inc = _mm_set1_epi64x((signed_size_t)(sizeof(m) / sizeof(uint64_t)));
    if((uint64_t)start % sizeof(m) == 0) {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint64_t))) {
            _mm_store_si128((__m128i *)start, m);
            m = _mm_add_epi64(m, inc);
            start += (signed_size_t)(sizeof(m) / sizeof(uint64_t));
        }
    } else {
        UNROLL4
        while(end - start >= (signed_size_t)(sizeof(m) / sizeof(uint64_t))) {
            _mm_storeu_si128((__m128i *)start, m);
            m = _mm_add_epi64(m, inc);
            start += (signed_size_t)(sizeof(m) / sizeof(uint64_t));
        }
    }
    sv += (start - origstart);
#endif
    while(start < end) *start++ = sv++;
    return origstart;
}


FASTIOTA_EXPORT int64_t *fastiota64i_inc(int64_t *start, size_t nelem, int64_t sv, int64_t inc)
{
    return (int64_t *)fastiota64_inc((uint64_t *)start, nelem, sv, inc);
}

FASTIOTA_EXPORT int32_t *fastiota32i_inc(int32_t *start, size_t nelem, int32_t sv, int32_t inc)
{
    return (int32_t *)fastiota32_inc((uint32_t *)start, nelem, sv, inc);
}

#ifdef __cplusplus
}
#endif
#undef init_start
#undef init_start64
#undef init_start_val
#undef init_start_val64

#endif /* FASTIOTA_C_H__ */
