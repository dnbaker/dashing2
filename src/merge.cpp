#include "merge.h"

namespace dashing2 {

template<typename KeyT, typename VT, typename Hash>
void merge_template(ska::flat_hash_map<KeyT, VT, Hash> &lhs, const ska::flat_hash_map<KeyT, VT, Hash> &rhs) {
    typename ska::flat_hash_map<KeyT, VT, Hash>::iterator it;
    for(const auto &pair: rhs)
        if((it = lhs.find(pair.first)) != lhs.end())
            it->second += pair.second;
}


void merge(ska::flat_hash_map<uint64_t, int32_t> &lhs, const ska::flat_hash_map<uint64_t, int32_t> &rhs) {
    merge_template(lhs, rhs);
}
void merge(ska::flat_hash_map<uint64_t, double> &lhs, const ska::flat_hash_map<uint64_t, double> &rhs) {
    merge_template(lhs, rhs);
}

void merge(ska::flat_hash_map<u128_t, int32_t, FHasher> &lhs, const ska::flat_hash_map<u128_t, int32_t, FHasher> &rhs) {
    merge_template(lhs, rhs);
}
void merge(ska::flat_hash_map<u128_t, double, FHasher> &lhs, const ska::flat_hash_map<u128_t, double, FHasher> &rhs) {
    merge_template(lhs, rhs);
}

}
