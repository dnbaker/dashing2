#pragma once
#include <cstdint>
#include <utility>
#include <cstdlib>
#include <cstdio>

namespace dashing2 {
using std::size_t;
using std::uint64_t;
// Compute histogram similarity, which normalizes to weighted Jaccard
// Returns {intersection_size, union_size}
std::pair<double, double> weighted_compare(const uint64_t *lptr, const double *lnptr, size_t lhl, const double lhsum, const uint64_t *rptr, const double *rnptr, size_t rhl, const double rhsum, bool kahan=false);
// Compute set similarity, which normalizes to Jaccard
size_t set_compare(const uint64_t *lptr, size_t lhl, const uint64_t *rptr, size_t rhl);
// Compute the dot product between k-mer sets, which can be used to compute cosine similarity/distance.
double cosine_compare(const uint64_t *lptr, size_t lhl, const double lhnorm, const uint64_t *rptr, size_t rhl, const double rhnorm, const double *lnptr, const double *rnptr, bool kahan=false);
std::pair<double, double> weighted_compare(std::FILE *lhk, std::FILE *rhk, std::FILE *lhn, std::FILE *rhn, double lhsum, double rhsum);

}


