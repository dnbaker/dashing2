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
std::pair<double, double> weighted_compare_mode(const uint64_t *lptr, size_t lhl, const uint64_t *rptr, size_t rhl, const double *lnptr, const double *rnptr) {
    double union_size = 0, isz_size = 0; // maybe Kahan update?
    double union_carry = 0., isz_carry = 0.;
    auto increment = [](auto &sum, auto &carry, auto inc) {
        if constexpr(MODE)
            sketch::kahan::update(sum, carry, inc);
        else
            sum += inc;
    };
    for(size_t lhi = 0, rhi = 0; lhi < lhl || rhi < rhl;) {
        if(lhi < lhl) {
            if(rhi < rhl) {
                if(lptr[lhi] == rptr[rhi]) {
                    double mnv = std::min(lnptr[lhi], rnptr[rhi]);
                    double mxv = std::max(lnptr[lhi], rnptr[rhi]);
                    increment(isz_size, isz_carry, mnv);
                    increment(union_size, union_carry, mxv);
                    ++lhi; ++rhi;
                } else if(lptr[lhi] < rptr[rhi]) {
                    increment(union_size, union_carry, lnptr[lhi++]);
                } else {
                    increment(union_size, union_carry, rnptr[rhi++]);
                }
            } else increment(union_size, union_carry, lnptr[lhi++]);
        } else if(rhi < rhl) increment(union_size, union_carry, rnptr[rhi++]);
    }
    return std::make_pair(isz_size, union_size);
}
std::pair<double, double> weighted_compare(const uint64_t *lptr, size_t lhl, const uint64_t *rptr, size_t rhl, double *lnptr, double *rnptr, bool kahan) {
    return kahan ? weighted_compare_mode<1>(lptr, lhl, rptr, rhl, lnptr, rnptr): weighted_compare_mode<0>(lptr, lhl, rptr, rhl, lnptr, rnptr);
}
size_t set_compare(const uint64_t *lptr, size_t lhl, const uint64_t *rptr, size_t rhl) {
    PushBackCounter c;
    std::set_intersection(lptr, lptr + lhl, rptr, rptr + rhl, std::back_inserter(c));
    return c.count;
}
}
