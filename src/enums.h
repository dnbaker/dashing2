#pragma once

namespace dashing2 {


enum DataType {
    FASTX,
    BIGWIG,
    BED
};


enum SketchSpace {
    SPACE_SET,      // MinHash/SetSketch/HLL
    SPACE_MULTISET, // Weighted MinHash -- e.g., Bag or Tree
    SPACE_PSET,      // ProbMinHash
    SPACE_EDIT_DISTANCE // edit distance -- implies OMH and forces us to use strings as input instead of bags of k-mers
};

enum CountingType {
    EXACT_COUNTING,
    COUNTSKETCH_COUNTING
};

#ifndef SKETCH_FLOAT_TYPE
#define SKETCH_FLOAT_TYPE double
#endif
using RegT = SKETCH_FLOAT_TYPE;
#undef SKETCH_FLOAT_TYPE

using u128_t = __uint128_t;

}
