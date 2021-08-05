#pragma once
#include <string>

namespace dashing2 {

#ifndef SKETCH_FLOAT_TYPE
#define SKETCH_FLOAT_TYPE double
#endif
using RegT = SKETCH_FLOAT_TYPE;
#undef SKETCH_FLOAT_TYPE

using u128_t = __uint128_t;



// Determines the file-type being consumed
enum DataType {
    FASTX,
    BIGWIG, //--bigwig
    BED,    //--bed
    LEAFCUTTER //--leafcutter
};


enum SketchSpace {
    SPACE_SET,      // MinHash/SetSketch/HLL
    SPACE_MULTISET, // Weighted MinHash -- e.g., Bag or Tree
    SPACE_PSET,      // ProbMinHash
    SPACE_EDIT_DISTANCE, // edit distance -- implies OMH and forces us to use strings as input instead of bags of k-mers
    EDIT_DISTANCE = SPACE_EDIT_DISTANCE,
    MULTISET = SPACE_MULTISET
};

enum CountingType {
    EXACT_COUNTING,
    COUNTSKETCH_COUNTING,
    // Add or substract each item to its assigned bucket using a random variable seeded by the key
    // Usually estimates higher-hitter featuers well
    // Pseudocode:
    // data[hash(key) % len(data)] += hash(key) & 1 ? inc: -inc;
    // Then use std::abs(data[i]) for each element after counting
    // This reduces the sample space at some inexactness, but the biggest elements will remain the biggest

    //Not implemented, but we may expand this:
    COUNTMIN_COUNTING,
    CQF_COUNTING // Not implemen
};

#define THROW_EXCEPTION(...) do {\
        auto exception__ = __VA_ARGS__;\
        std::cerr << "Exception " << exception__.what() << " from " << std::this_thread::get_id() << '\n';\
        throw exception__;\
    } while(0)


enum KmerSketchResultType {
    ONE_PERM = 0,       // Faster (3-4x) than Full, comparable accuracy for both cardinality and set similarities
    // This is a stochastically-averaged generalized HyperLogLog
    // Constant-time updates,
    FULL_SETSKETCH, // Not stochastically-averaged; potentially better LSH properties
    // This is a generalized HyperLogLog
    FULL_MMER_SET,
    /*
     * Convert the genome into a k-mer set; uses a hash table, which is then converted into a sorted hash list
    */
    FULL_MMER_COUNTDICT,
    /*
        Convert into a k-mer: count dictionary.
    */
    FULL_MMER_SEQUENCE,
    /*
    Convert the genome into a list of k-mers/minimizers; could be used for minimizer index generation
    */
};

enum OutputKind {
    SYMMETRIC_ALL_PAIRS,
    ASYMMETRIC_ALL_PAIRS,
    KNN_GRAPH, // Fixed top-k neighbors
    NN_GRAPH_THRESHOLD, // Variable number of similarities, as given by threshold
    PANEL,
    DEDUP
};

enum OutputFormat {
    MACHINE_READABLE,
    HUMAN_READABLE,
    BINARY = MACHINE_READABLE
};

std::string to_string(KmerSketchResultType t);
std::string to_string(SketchSpace ss);
std::string to_string(DataType dt);
std::string to_string(CountingType ct);
std::string to_string(OutputKind ok);
std::string to_string(OutputFormat of);
std::string trim_folder(const std::string &s);
struct Dashing2Options;

std::string to_suffix(const Dashing2Options &opts);

}
