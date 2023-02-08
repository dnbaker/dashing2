#pragma once
#include <cstdint>
#include <utility>
#include <cstdlib>
#include <cstdio>
#include "enums.h"

namespace dashing2 {
using std::size_t;
using std::uint64_t;
// Compute histogram similarity, which normalizes to weighted Jaccard
// Returns {intersection_size, union_size}
std::pair<double, double> weighted_compare(const uint64_t *lptr, const double *lnptr, size_t lhl, const double lhsum, const uint64_t *rptr, const double *rnptr, size_t rhl, const double rhsum, bool kahan=false) noexcept;
std::pair<double, double> weighted_compare(const u128_t *lptr, const double *lnptr, size_t lhl, const double lhsum, const u128_t *rptr, const double *rnptr, size_t rhl, const double rhsum, bool kahan=true) noexcept;
// Compute set similarity, which normalizes to Jaccard
size_t set_compare(const uint64_t *lptr, size_t lhl, const uint64_t *rptr, size_t rhl) noexcept;
// Compute the dot product between k-mer sets, which can be used to compute cosine similarity/distance.
double cosine_compare(const uint64_t *lptr, size_t lhl, const double lhnorm, const uint64_t *rptr, size_t rhl, const double rhnorm, const double *lnptr, const double *rnptr, bool kahan=false) noexcept;
size_t hamming_compare(const uint64_t *SK_RESTRICT lptr, size_t lhl, const uint64_t *SK_RESTRICT rptr, size_t rhl) noexcept;

std::pair<double, double> weighted_compare(std::FILE *lhk, std::FILE *rhk, std::FILE *lhn, std::FILE *rhn, double lhsum, double rhsum, bool use128) noexcept;

size_t hamming_compare_f64(std::FILE *lfp, std::FILE *rfp) noexcept;
size_t hamming_compare_f128(std::FILE *lfp, std::FILE *rfp)noexcept;
std::pair<size_t, size_t> mmer_edit_distance(std::FILE *lfp, std::FILE *rfp, bool use128=true) noexcept;
// Computes edit distance between two file pointer's data
// If true, uses __uint128_t
// otherwise, uint64_t

}


