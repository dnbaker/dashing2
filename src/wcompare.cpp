#include "wcompare.h"
#include <memory>
#include "sketch/kahan.h"
#include <algorithm>
#include <cassert>

namespace dashing2 {
using u128_t = __uint128_t;

struct PushBackCounter
{
  struct value_type { template<typename T> value_type(const T&) { } };
  void push_back(const value_type &){++count;}
  PushBackCounter(): count(0) {}
  size_t count;
};

template<size_t MODE=0>
std::pair<double, double> weighted_compare_mode(const uint64_t *lptr, size_t lhl, const double lhsum, const uint64_t *rptr, size_t rhl, const double rhsum, const double *lnptr, const double *rnptr) {
    double isz_size = 0, isz_carry = 0.;
    auto increment = [](auto &sum, auto &carry, auto inc) __attribute__((__always_inline__)) {
        if(MODE) sketch::kahan::update(sum, carry, inc);
        else    sum += inc;
    };
    for(size_t lhi = 0, rhi = 0; lhi < lhl && rhi < rhl;) {
        if(lptr[lhi] == rptr[rhi]) {
            increment(isz_size, isz_carry, std::min(lnptr[lhi++], rnptr[rhi++]));
        } else {
            const bool v = lptr[lhi] < rptr[rhi];
            lhi += v; rhi += !v;
        }
    }
    return std::make_pair(isz_size, lhsum + rhsum - isz_size);
}
std::pair<double, double> weighted_compare(const uint64_t *lptr, size_t lhl, const double lhsum, const uint64_t *rptr, size_t rhl, const double rhsum, const double *lnptr, const double *rnptr, bool kahan) {
    return kahan ? weighted_compare_mode<1>(lptr, lhl, lhsum, rptr, rhl, rhsum, lnptr, rnptr): weighted_compare_mode<0>(lptr, lhl, lhsum, rptr, rhl, rhsum, lnptr, rnptr);
}
size_t set_compare(const uint64_t *lptr, size_t lhl, const uint64_t *rptr, size_t rhl) {
    PushBackCounter c;
    std::set_intersection(lptr, lptr + lhl, rptr, rptr + rhl, std::back_inserter(c));
    return c.count;
}

double cosine_compare(const uint64_t *lptr, size_t lhl, [[maybe_unused]] const double lhnorm, const uint64_t *rptr, size_t rhl, [[maybe_unused]] const double rhnorm, const double *lnptr, const double *rnptr, bool kahan) {
    double dotprod = 0, carry = 0.;
    auto increment = [kahan](auto &norm, auto &carry, auto inc) __attribute__((__always_inline__)) {
        if(kahan) sketch::kahan::update(norm, carry, inc); else norm += inc;
    };
    // TODO:
    for(size_t lhi = 0, rhi = 0; lhi < lhl && rhi < rhl;) {
        if(lptr[lhi] == rptr[rhi]) {
            increment(dotprod, carry, lnptr[lhi++] * rnptr[rhi++]);
        } else {
            const bool v = lptr[lhi] < rptr[rhi];
            lhi += v; rhi += !v;
        }
    }
    return dotprod;
}
template<typename T>
std::vector<T> load_file(std::FILE *fp) {
    std::unique_ptr<T[]> tmp(new T[16384]);
    std::vector<T> ret;
    while(!std::feof(fp)) {
        size_t n = std::fread(tmp.get(), sizeof(T), 16384, fp);
        const size_t oldsz = ret.size();
        ret.resize(oldsz + n);
        std::copy(tmp.get(), &tmp[n], &ret[oldsz]);
        if(n != 16384)
            break;
    }
    return ret;
}
std::pair<double, double>
weighted_compare(std::FILE *lhk, std::FILE *rhk, std::FILE *lhn, std::FILE *rhn, double lhsum, double rhsum) {
    //auto [isz_size, union_size] = weighted_compare(lhk, rhk, lhn, rhn, lhc, rhc);
    uint64_t lhv, rhv;
    double lhc = 1., rhc = 1.;
    if(std::feof(lhk) || std::feof(rhk)) {
        // If one of the sets is empty, return 0 intersection and union size accordingly
        return {0., lhsum + rhsum};
    }
    auto incl = [&]() {
        assert(lhk);
        std::fread(&lhv, sizeof(lhv), 1, lhk); if(lhn) std::fread(&lhc, sizeof(lhc), 1, lhn);
    };
    auto incr = [&]() {std::fread(&rhv, sizeof(rhv), 1, rhk); if(rhn) std::fread(&rhc, sizeof(rhc), 1, rhn);};
    auto incb = [&]() {incl(); incr();};
    long double isz = 0.;
    incb();
    for(;;) {
        if(lhv < rhv) incl();
            // lhv not found, increment lh side
        else if(rhv < lhv) incr();
        else {
            isz += std::min(rhc, lhc);
            incb();
        }
        if(std::feof(lhk) || std::feof(rhk)) break;
    }
    std::pair<double, double> ret;
    ret.first = isz;
    ret.second = lhsum + rhsum - isz;
    return ret;
}


template<typename T, typename CT>
double cosine_compare_(std::FILE *lhk, std::FILE *rhk, std::FILE *lhn, std::FILE *rhn) {
    // Returns cosine similariy between the vectors
    //auto [isz_size, union_size] = weighted_compare(lhk, rhk, lhn, rhn, lhc, rhc);
    T lhv, rhv;
    CT lhc = 1., rhc = 1.;
    if(std::feof(lhk) || std::feof(rhk)) return 0.;

    auto incl = [&]() {std::fread(&lhv, sizeof(lhv), 1, lhk); if(lhn) std::fread(&lhc, sizeof(lhc), 1, lhn);};
    auto incr = [&]() {std::fread(&rhv, sizeof(rhv), 1, rhk); if(rhn) std::fread(&rhc, sizeof(rhc), 1, rhn);};
    auto incb = [&]() {incl(); incr();};
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
    return isz;
}

double cosine_compare(std::FILE *lhk, std::FILE *rhk, std::FILE *lhn, std::FILE *rhn, bool use128keys=false, bool usef32counts=false) {
    const auto ptr = use128keys ? usef32counts ? &cosine_compare_<u128_t, float>: &cosine_compare_<u128_t, double>: usef32counts ? &cosine_compare_<uint64_t, float>: &cosine_compare_<uint64_t, double>;
    return ptr(lhk, rhk, lhn, rhn);
}


} // namespace dashing2
