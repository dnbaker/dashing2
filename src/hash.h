#pragma once
#include "sketch/hash.h"
#include "enums.h"
namespace dashing2 {

namespace hash = sketch::hash;
struct FHasher {
    using FastRevHash = hash::CEHasher;
    FastRevHash rhasher_;
    FHasher() {}
    INLINE uint64_t operator()(u128_t x) const {
        return rhasher_(uint64_t(x>>64)) ^ rhasher_(uint64_t(x));
    }
};

}
