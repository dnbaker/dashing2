#pragma once
#include <string>

namespace dashing2 {


enum DataType {
    FASTX,
    BIGWIG,
    BED,
    LEAFCUTTER
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
    // Add or substract each item to its assigned bucket using a random variable seeded by the key
    // Usually estimates higher-hitter featuers well
    // Pseudocode:
    // data[hash(key) % len(data)] += hash(key) & 1 ? inc: -inc;
    // Then use std::abs(data[i]) for each element after counting
    // This reduces the sample space at some inexactness, but the biggest elements will remain the biggest
};

#ifndef SKETCH_FLOAT_TYPE
#define SKETCH_FLOAT_TYPE double
#endif
using RegT = SKETCH_FLOAT_TYPE;
#undef SKETCH_FLOAT_TYPE

using u128_t = __uint128_t;

enum KmerSketchResultType {
    ONE_PERM = 0,       // Faster (3-4x) than Full, comparable accuracy for both cardinality and set similarities
    // This is a stochastically-averaged generalized HyperLogLog
    // Constant-time updates,
    FULL_SETSKETCH = 1, // Not stochastically-averaged; potentially better LSH properties
    // This is a generalized HyperLogLog
    FULL_MMER_SET = 2,
    /*
     * Convert the genome into a k-mer set; uses a hash table
    */
    FULL_MMER_SEQUENCE = 3,
    /*
    Convert the genome into a list of k-mers/minimizers; could be used for minimizer index generation
    */
    FULL_MMER_COUNTDICT = 4
    /*
        Convert into a k-mer: count dictionary.
    */
};
static inline std::string to_string(KmerSketchResultType t) {
    if(t == ONE_PERM) {return "OnePermutationSetSketch";}
    if(t == FULL_SETSKETCH) return "FullSetSketch";
    if(t == FULL_MMER_SET) return "FullMmerSet";
    if(t == FULL_MMER_SEQUENCE) return "FullMmerSequence";
    return "FullMmerCountdict";
}


}
