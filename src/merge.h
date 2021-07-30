#pragma once
#include "./hash.h"
#include "counter.h"

namespace dashing2 {

void merge(flat_hash_map<uint64_t, int32_t> &lhs, const flat_hash_map<uint64_t, int32_t> &rhs);
void merge(flat_hash_map<uint64_t, double> &lhs, const flat_hash_map<uint64_t, double> &rhs);
void merge(flat_hash_map<u128_t, int32_t, FHasher> &lhs, const flat_hash_map<u128_t, int32_t, FHasher> &rhs);
void merge(flat_hash_map<u128_t, double, FHasher> &lhs, const flat_hash_map<u128_t, double, FHasher> &rhs);

}
