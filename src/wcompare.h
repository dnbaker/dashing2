#pragma once
#include <cstdint>
#include <utility>
#include <cstdlib>

namespace dashing2 {
using std::size_t;
using std::uint64_t;
std::pair<double, double> weighted_compare(const uint64_t *lhid, size_t lhn, const uint64_t *rhid, size_t rhn, const double *lhcp, const double *rhcp, int mode=1);
size_t set_compare(const uint64_t *, size_t, const uint64_t *, size_t);
}


