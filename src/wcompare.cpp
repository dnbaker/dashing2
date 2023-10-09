#include <memory>
#include "sketch/kahan.h"
#include <algorithm>
#include <cassert>
#include "wcompare.h"
#include <numeric>
#include "enums.h"
#include <cstring>
#include "edlib.h"

#include <span>

namespace dashing2 {
//using u128_t = __uint128_t;
template<typename T, size_t BUFSZ=16384> std::vector<T> load_file(std::FILE *fp);

struct PushBackCounter

{
  struct value_type { template<typename T> value_type(const T&) { } };
  void push_back(const value_type &){++count;}
  PushBackCounter(): count(0) {}
  size_t count;
};

template<typename IT, size_t MODE=0>
std::pair<double, double> weighted_compare_mode(const IT *lptr, size_t lhl, const double lhsum, const double *lnptr, const IT *rptr, const double *rnptr, size_t rhl, const double rhsum) noexcept {
    double isz_size = 0, isz_carry = 0.;
    for(size_t lhi = 0, rhi = 0; lhi < lhl && rhi < rhl;) {
        if(lptr[lhi] == rptr[rhi]) {
            const auto increment = std::min(lnptr[lhi++], rnptr[rhi++]);
            if(MODE) {
                sketch::kahan::update(isz_size, isz_carry, increment);
            } else {
                isz_size += increment;
            }
        } else {
            const bool v = lptr[lhi] < rptr[rhi];
            lhi += v; rhi += !v;
        }
    }
    return std::make_pair(isz_size, lhsum + rhsum - isz_size);
}
std::pair<double, double> weighted_compare(const uint64_t *lptr, const double *lnptr, size_t lhl, const double lhsum, const uint64_t *rptr, const double *rnptr, size_t rhl, const double rhsum, bool kahan) noexcept {
    return kahan ? weighted_compare_mode<uint64_t, 1>(lptr, lhl, lhsum, lnptr, rptr, rnptr, rhl, rhsum): weighted_compare_mode<uint64_t, 0>(lptr, lhl, lhsum, lnptr, rptr, rnptr, rhl, rhsum);
}
std::pair<double, double> weighted_compare(const u128_t *lptr, const double *lnptr, size_t lhl, const double lhsum, const u128_t *rptr, const double *rnptr, size_t rhl, const double rhsum, bool kahan) noexcept {
    auto ptr = kahan ? &weighted_compare_mode<u128_t, 1> : &weighted_compare_mode<u128_t, 0>;
    return ptr(lptr, lhl, lhsum, lnptr, rptr, rnptr, rhl, rhsum);
}
size_t hamming_compare(const uint64_t *SK_RESTRICT lptr, size_t lhl, const uint64_t *SK_RESTRICT rptr, size_t rhl) noexcept {
    return std::inner_product(lptr, lptr + std::min(lhl, rhl), rptr, size_t(0), [](auto c, auto x) noexcept {return c += x;}, [](auto lhs, auto rhs) noexcept -> size_t {return lhs == rhs;})
            + (std::max(lhl, rhl) - std::min(lhl, rhl));
}

static constexpr EdlibAlignConfig CONFIG {
    .k = -1,
    .mode = EDLIB_MODE_NW,
    .task = EDLIB_TASK_DISTANCE,
    .additionalEqualities = nullptr,
    .additionalEqualitiesLength = 0,
};

template<typename T>
int32_t edlib_edit_distance(const T&x, const T& y) {
    const char *xptr = reinterpret_cast<const char *>(x.data());
    const char *yptr = reinterpret_cast<const char *>(y.data());
    const int64_t xsize = x.size();
    const int64_t ysize = y.size();
    EdlibAlignResult result = edlibAlign(xptr, xsize, yptr, ysize, CONFIG);
    return result.editDistance;
}

template<typename T>
std::pair<size_t, size_t> mmer_edit_distance_f(std::FILE *lfp, std::FILE *rfp) noexcept {
    // Compute exact edit distance at char level, then divide by the ratio of the size of these registers to characters
    // IE, the edit distance calculated is by chars rather than T
    const auto lhs(load_file<T>(lfp));
    const auto rhs(load_file<T>(rfp));
    const std::span lh_span((const uint8_t *)lhs.data(), (const uint8_t *)(lhs.data() + lhs.size()));
    const std::span rh_span((const uint8_t *)rhs.data(), (const uint8_t *)(rhs.data() + rhs.size()));
    size_t ret = edlib_edit_distance(lh_span, rh_span);
    return {ret / sizeof(T), std::max(lhs.size(), rhs.size())};
}

std::pair<size_t, size_t> mmer_edit_distance(std::FILE *lfp, std::FILE *rfp, bool use128) noexcept {
    auto ptr = use128 ? &mmer_edit_distance_f<u128_t>: &mmer_edit_distance_f<uint64_t>;
    return ptr(lfp, rfp);
}
template<size_t NB>
size_t hamming_compare_f(std::FILE *lfp, std::FILE *rfp) {
    using RT = std::conditional_t<NB == 1, uint8_t, std::conditional_t<NB == 2, uint16_t, std::conditional_t<NB == 4, uint32_t, std::conditional_t<NB == 8, uint64_t, u128_t>>>>;
    size_t ret = 0;
    RT lv, rv;
    for(;;) {
        if(std::feof(lfp)) {
            if(!std::feof(rfp))
                for(;std::fread(&rv, NB, 1, rfp) == 1u;++ret);
            break;
        } else if(std::feof(rfp)) {
            while(std::fread(&lv, NB, 1, lfp) == 1u) ++ret;
            break;
        } else {
            std::fread(&rv, NB, 1, rfp);
            std::fread(&lv, NB, 1, lfp);
            ret += lv == rv;
        }
    }
    return ret;
}
size_t hamming_compare_f64(std::FILE *lfp, std::FILE *rfp) noexcept {return hamming_compare_f<8>(lfp, rfp);}
size_t hamming_compare_f128(std::FILE *lfp, std::FILE *rfp) noexcept {return hamming_compare_f<16>(lfp, rfp);}
size_t set_compare(const uint64_t *lptr, size_t lhl, const uint64_t *rptr, size_t rhl) noexcept {
    PushBackCounter c;
    std::set_intersection(lptr, lptr + lhl, rptr, rptr + rhl, std::back_inserter(c));
    return c.count;
}

double cosine_compare(const uint64_t *lptr, size_t lhl, [[maybe_unused]] const double lhnorm, const uint64_t *rptr, size_t rhl, [[maybe_unused]] const double rhnorm, const double *lnptr, const double *rnptr, bool kahan) noexcept {
    double dotprod = 0, carry = 0.;
    auto increment = [kahan](auto &norm, auto &carry, auto inc) __attribute__((__always_inline__)) {
        if(kahan) sketch::kahan::update(norm, carry, inc); else norm += inc;
    };
    for(size_t lhi = 0, rhi = 0; lhi < lhl && rhi < rhl;) {
        if(lptr[lhi] == rptr[rhi]) increment(dotprod, carry, lnptr[lhi++] * rnptr[rhi++]);
        else {
            const bool v = lptr[lhi] < rptr[rhi];
            lhi += v; rhi += !v;
        }
    }
    return dotprod;
}

template<typename T, size_t BUFSZ>
std::vector<T> load_file(std::FILE *fp) {
    std::unique_ptr<T[]> tmp(new T[BUFSZ]);
    std::vector<T> ret;
    while(!std::feof(fp)) {
        size_t n = std::fread(tmp.get(), sizeof(T), BUFSZ, fp);
        const size_t oldsz = ret.size();
        ret.resize(oldsz + n);
        std::copy(tmp.get(), &tmp[n], &ret[oldsz]);
        if(n != BUFSZ) break;
    }
    return ret;
}
std::pair<double, double>
weighted_compare(std::FILE *lhk, std::FILE *rhk, std::FILE *lhn, std::FILE *rhn, double lhsum, double rhsum, bool use128) noexcept {
    //auto [isz_size, union_size] = weighted_compare(lhk, rhk, lhn, rhn, lhc, rhc);

    u128_t lhv2, rhv2;
    uint64_t lhv, rhv;
    //uint64_t &lhv = *reinterpret_cast<uint64_t *>(&lhv2), &rhv = *reinterpret_cast<uint64_t *>(&rhv2);
    double lhc = 1., rhc = 1.;
    if(std::feof(lhk) || std::feof(rhk)) {
        return {0., lhsum + rhsum};
    }
    // Skip the cardinalities of both
    const size_t incbytes = use128 ? 16: 8;
    auto lh2p = &lhv2; auto lh1p = &lhv;
    auto rh1p = &rhv;  auto rh2p = &rhv2;
    void *const lhp = use128 ? static_cast<void *>(lh2p): static_cast<void *>(lh1p);
    void *const rhp = use128 ? static_cast<void *>(rh2p): static_cast<void *>(rh1p);
#define incl() do {\
        std::fread(lhp, incbytes, 1, lhk); \
        if(lhn) std::fread(&lhc, sizeof(lhc), 1, lhn);\
    } while(0)
#define incr() do {\
        std::fread(rhp, incbytes, 1, rhk);\
        if(rhn) std::fread(&rhc, sizeof(rhc), 1, rhn);\
    } while(0)
#define incb() {incl(); incr();}
    long double isz = 0.;
    incb();
    for(;;) {
        if(use128 ? lhv2 < rhv2 :lhv < rhv) {
            incl();
            // lhv not found, increment lh side
        }
        else if(use128 ? rhv2 < lhv2: rhv < lhv) {
            incr();
        } else {
            isz += std::min(rhc, lhc);
            incb();
        }
        if(std::feof(lhk) || std::feof(rhk)) break;
    }
    return {double(isz), double(lhsum + rhsum - isz)};
}


template<typename T, typename CT>
double cosine_compare_(std::FILE *lhk, std::FILE *rhk, std::FILE *lhn, std::FILE *rhn) noexcept {
    // Returns cosine similariy between the vectors
    //auto [isz_size, union_size] = weighted_compare(lhk, rhk, lhn, rhn, lhc, rhc);
    T lhv, rhv;
    CT lhc = 1., rhc = 1.;
    if(std::feof(lhk) || std::feof(rhk)) return 0.;
    // Skip cardinalities of both files.
    static constexpr size_t incbytes = sizeof(T);
    const auto lhp = &lhv, rhp = &rhv;

    auto isz = 0.L;
    incb();
    for(auto carry = 0.L;;) {
        if(lhv < rhv) incl();
            // lhv not found, increment lh side
        else if(rhv < lhv) incr();
        else {
            sketch::kahan::update(isz, carry, static_cast<long double>(rhc) * static_cast<long double>(lhc));
            incb();
        }
        if(std::feof(lhk) || std::feof(rhk)) break;
    }
#undef incb
#undef incl
#undef incr
    return isz;
}

double cosine_compare(std::FILE *lhk, std::FILE *rhk, std::FILE *lhn, std::FILE *rhn, bool use128keys=false, bool usef32counts=false) noexcept {
    const auto ptr = use128keys ? usef32counts ? &cosine_compare_<u128_t, float>: &cosine_compare_<u128_t, double>: usef32counts ? &cosine_compare_<uint64_t, float>: &cosine_compare_<uint64_t, double>;
    return ptr(lhk, rhk, lhn, rhn);
}

#ifdef WCOMPARE_MAIN

int main(){
}

#endif


} // namespace dashing2
