#pragma once
#include "flat_hash_map/flat_hash_map.hpp"
#include "./hash.h"

namespace dashing2 {

void merge(ska::flat_hash_map<uint64_t, int32_t> &lhs, const ska::flat_hash_map<uint64_t, int32_t> &rhs);
void merge(ska::flat_hash_map<uint64_t, double> &lhs, const ska::flat_hash_map<uint64_t, double> &rhs);
void merge(ska::flat_hash_map<u128_t, int32_t, FHasher> &lhs, const ska::flat_hash_map<u128_t, int32_t, FHasher> &rhs);
void merge(ska::flat_hash_map<u128_t, double, FHasher> &lhs, const ska::flat_hash_map<u128_t, double, FHasher> &rhs);

}
