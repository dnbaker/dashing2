#pragma once
#include "sketch/hash.h"
#include "enums.h"
namespace dashing2 {

namespace hash = sketch::hash;
struct FHasher {
    using FastRevHash = hash::FusedReversible3<hash::InvMul, hash::RotL33, hash::MultiplyAddXoRot<16>>;
    FastRevHash rhasher_;
    FHasher() {}
    INLINE uint64_t operator()(u128_t x) const {
        return rhasher_(uint64_t(x>>64)) ^ rhasher_(uint64_t(x));
    }
};

}
