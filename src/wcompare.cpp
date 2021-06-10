#include "wcompare.h"
#include "sketch/kahan.h"
#include <algorithm>

namespace dashing2 {

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
    // TODO:
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

} // namespace dashing2
