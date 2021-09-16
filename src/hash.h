#pragma once
#include "sketch/hash.h"
#include "enums.h"
#include "robin_hood.h"
namespace dashing2 {

namespace hash = sketch::hash;
struct FHasher {
    using FastRevHash = hash::CEHasher;
    FastRevHash rhasher_;
    FHasher() {}
    template<typename T>
    INLINE decltype(auto) hash(T x) const {
        return this->operator()(x);
    }
    template<typename T>
    static INLINE decltype(auto) hashi(T x) {
        return FHasher()(x);
    }
    INLINE uint32_t operator()(uint32_t x) const {
        return rhasher_(x);
    }
    INLINE uint64_t operator()(uint64_t x) const {
        return rhasher_(x);
    }
    INLINE uint64_t operator()(u128_t x) const {
        return rhasher_(uint64_t(x>>64)) ^ rhasher_(uint64_t(x));
    }
};

template<typename Key, typename V,
         typename Hash=std::conditional_t<std::is_same_v<u128_t, Key>, FHasher, robin_hood::hash<Key>>>
using flat_hash_map = robin_hood::unordered_flat_map<Key, V, Hash>;

template<typename Key,
         typename Hash=std::conditional_t<
                std::is_same_v<u128_t, Key>, FHasher, robin_hood::hash<Key>
         >
        >
using flat_hash_set = robin_hood::unordered_flat_set<Key>;

}
