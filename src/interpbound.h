#ifndef INTERPBOUND_H__
#define INTERPBOUND_H__
#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <utility>
#include <type_traits>
#include "sketch/macros.h"

namespace interp {
using std::ptrdiff_t;

template<typename T>
T geomed(T low, T hi) {
    std::conditional_t<(sizeof(T) == 4), float, std::conditional_t<(sizeof(T) == 8), double, long double>> lit, hit, sum;
    std::memcpy(&lit, &low, sizeof(lit));
    std::memcpy(&hit, &hi, sizeof(lit));
    sum = (lit + hit) / 2;
    std::memcpy(&low, &sum, sizeof(lit));
    return low;
}

// if arithmetic is false, interpolate in linearly
// otherwise, interpolate geometrically
template<bool arithmetic=true, typename IT, typename KeyT>
std::pair<IT, bool> search(const IT beg, const IT end, KeyT key) {
    if(beg == end) return {beg, false};
    if(end - beg == 1) {
        return {beg, *beg == key};
    }
    //assert(std::is_sorted(beg, end));
    assert(beg < end);
    if(beg == end) return {end, false};
    using CT = long double;//typename std::conditional<(sizeof(KeyT) > 8), long double, double>::type;
    //DBG_ONLY(std::fprintf(stderr, "[%s] span of %zu elements from %g to %g\n", __PRETTY_FUNCTION__, end - beg, double(*beg), double(*(end - 1)));)
    KeyT max = end[-1];
    KeyT min = *beg;
    //DBG_ONLY(std::fprintf(stderr, "range: %zd. min: %g. max: %g. key: %g\n", end - beg, double(min), double(max), double(key));)
    assert(min <= max);
    if(key > max) return {end, false};
    if(key < min) return {beg, false};
    IT upper = end - 1;
    IT lower = beg;
    //if(key == *lower) return {lower, true};
    //if(key == *upper) return {upper, true};
    IT cid;
    auto mydiff = key - min;
    auto mdiff = max - min;
    CT range = end - beg;
    if constexpr(arithmetic) {
        cid = std::min(lower + std::ptrdiff_t(std::ceil(CT(mydiff) / mdiff * range)), upper);
    } else {
        CT top = mdiff;
        CT bottom = mydiff;
        CT midpoint = std::sqrt(top / bottom) * bottom;
        std::ptrdiff_t median = std::ceil(midpoint / max * range);
        cid = std::min(lower + median, upper);
    }
    //std::fprintf(stderr, "key %zu vs min %zu vs max %zu, yielding %zu/%zu\n", key, min, max, cid - beg, end - beg);
    DBG_ONLY(size_t iternum = 0;)
    for(;;) {
        DBG_ONLY(
            ++iternum;
            std::fprintf(stderr, "Iter %zu, %zd/%zu. key: %g vs current point %g\n", iternum, cid - beg, end - beg, double(key), double(*cid));
        )
        assert(cid <= end);
        assert(cid >= beg || !std::fprintf(stderr, "beg = %zd, end = %zd, pos = %zd\n", size_t(0), end - beg, cid - beg));
        CT value = *cid;
        if(*cid == key) {
            //DBG_ONLY(std::fprintf(stderr, "%zu iters\n", iternum);)
            return {cid, true};
        }
        if(*cid > key) {
            //std::fprintf(stderr, "Lowering upper bound\n");
            upper = cid - 1;
            max = value;
        } else {
            //(*cid < key)
            lower = std::min(cid + 1, end);
            min = value;
        }
        assert(key >= min || key <= max);
        const ptrdiff_t range = upper - lower;
        if constexpr(arithmetic) {
            long double rat = (long double)(key - min) / (max - min);
            //std::fprintf(stderr, "rat: %Lg. num: %Lg. denom: %Lg\n", rat, CT(key - min), CT(max - min));
            ptrdiff_t diff = std::min(std::ptrdiff_t(std::ceil(rat * (range))), range);
            if(upper < lower) {
                std::fprintf(stderr, "Upper is lower by %zd...\n", lower - upper);
                return {cid, false};
            }
            cid = lower + diff;
        } else {
            auto top = key - min, bottom = max - min;
            CT midpoint = geomed(top, bottom);
            ptrdiff_t pos = std::ceil(midpoint / max * range);
            cid = std::max(std::min(lower + pos, upper), lower);
        }
    }
}
template<typename IT, typename KeyT>
IT find(IT beg, IT end, KeyT key) {
    auto [it, v] = search(beg, end, key);
    return v ? it: end;
}

} // namespace interp

#endif
