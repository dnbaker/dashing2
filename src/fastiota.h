#ifndef FAST_IOTA_H__
#define FAST_IOTA_H__
#include <stdint.h>
#include <stdlib.h>
#ifdef SINGLE_HEADER_FAST_IOTA
#define FASTIOTA_EXPORT static inline
#else
#define FASTIOTA_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif
FASTIOTA_EXPORT uint32_t *fastiota32(uint32_t *start, size_t nelem, uint32_t sv);
FASTIOTA_EXPORT int32_t *fastiota32i(int32_t *start, size_t nelem, int32_t sv);
FASTIOTA_EXPORT uint64_t *fastiota64(uint64_t *start, size_t nelem, uint64_t sv);
FASTIOTA_EXPORT int64_t *fastiota64i(int64_t *start, size_t nelem, int64_t sv);
FASTIOTA_EXPORT uint32_t *fastiota32_inc(uint32_t *start, size_t nelem, uint32_t sv, uint32_t inc);
FASTIOTA_EXPORT int32_t *fastiota32i_inc(int32_t *start, size_t nelem, int32_t sv, int32_t inc);
FASTIOTA_EXPORT uint64_t *fastiota64_inc(uint64_t *start, size_t nelem, uint64_t sv, uint64_t inc);
FASTIOTA_EXPORT int64_t *fastiota64i_inc(int64_t *start, size_t nelem, int64_t sv, int64_t inc);
#ifdef __cplusplus
}
#endif

#ifdef SINGLE_HEADER_FAST_IOTA
#include "./fastiota.c"
#endif

#ifdef __cplusplus
#include <array>

namespace fastiota {
template<typename IT, typename IT2=IT, typename IT3=IT2>
IT *iota(IT *start, size_t nelem, IT2 firstv=0, IT3 inc=1) {
    if(sizeof(IT) == 4) {
        return reinterpret_cast<IT *>(fastiota32_inc(reinterpret_cast<uint32_t *>(start), nelem, firstv, inc));
    }
    if(sizeof(IT) == 8) {
        return reinterpret_cast<IT *>(fastiota64_inc(reinterpret_cast<uint64_t *>(start), nelem, firstv, inc));
    }
    IT *ptr = start;
    const IT *const end = start + nelem;
    if(inc == 1) {
        while(ptr + 8 <= end) {
            ptr[0] = firstv;
            ptr[1] = firstv + 1;
            ptr[2] = firstv + 2;
            ptr[3] = firstv + 3;
            ptr[4] = firstv + 4;
            ptr[5] = firstv + 5;
            ptr[6] = firstv + 6;
            ptr[7] = firstv + 7;
            ptr += 8;
            firstv += 8;
        }
        while(ptr < end) *ptr++ = firstv++;
    } else {
        const std::array<const IT, 8> tmp {IT(0), IT(1) * inc, IT(2) * inc, IT(3) * inc, IT(4) * inc, IT(5) * inc, IT(6) * inc, IT(7) * inc};
        for(const IT blockinc = 8 * inc;ptr + 8 <= end;ptr += 8, firstv += blockinc) {
            ptr[0] = firstv;
            ptr[1] = firstv + tmp[1];
            ptr[2] = firstv + tmp[2];
            ptr[3] = firstv + tmp[3];
            ptr[4] = firstv + tmp[4];
            ptr[5] = firstv + tmp[5];
            ptr[6] = firstv + tmp[6];
            ptr[7] = firstv + tmp[7];
        }
        while(ptr < end) *ptr++ = (firstv += inc);
    }
    return start;
}
template<typename Iterator, typename Integer, typename Integer2>
void fastiota(Iterator start, Iterator stop, Integer sv=0) {
    iota(&*start, std::distance(start, stop), sv);
}
template<typename IT, typename IT2=IT, typename IT3=IT2, typename IT4=IT3>
IT *iota(IT *start, IT2 *stop, IT3 firstv=0, IT4 inc=1) {
    return iota(start, stop - start, firstv, inc);
}
} // namespace fastiota
#endif

#endif /* FAST_IOTA_H__ */
