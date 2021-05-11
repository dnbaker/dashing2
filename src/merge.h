#pragma once
#include "flat_hash_map/flat_hash_map.hpp"

namespace dashing2 {
template<typename KeyT, typename VT, typename Hash>
void merge(ska::flat_hash_map<KeyT, VT, Hash> &lhs, const ska::flat_hash_map<KeyT, VT, Hash> &rhs) {
    typename ska::flat_hash_map<KeyT, VT, Hash>::iterator it;
    for(const auto &pair: rhs)
        if((it = lhs.find(pair.first)) != lhs.end())
            it->second += pair.second;
}

}
