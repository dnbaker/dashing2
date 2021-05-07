#ifndef D2_H_H___
#define D2_H_H___
#include <memory>
#include <vector>
#include "bonsai/encoder.h"
#include "bigWig.h"
#include "flat_hash_map/flat_hash_map.hpp"
#include "xxHash/xxh3.h"
#include "sketch/setsketch.h"
#include "sketch/bmh.h"

#ifndef SKETCH_FLOAT_TYPE
#define SKETCH_FLOAT_TYPE double
#endif
using RegT = SKETCH_FLOAT_TYPE;

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
using u128_t = __uint128_t;

using namespace sketch;
namespace hash = sketch::hash;
using FastRevHash = hash::FusedReversible3<hash::InvMul, hash::RotL33, hash::MultiplyAddXoRot<16>>;

struct FHasher {
    FastRevHash rhasher_;
    FHasher() {}
    INLINE uint64_t operator()(u128_t x) const {
        return rhasher_(uint64_t(x>>64)) ^ rhasher_(uint64_t(x));
    }
};

struct Counter {
    CountingType ct_;
    ska::flat_hash_map<uint64_t, uint32_t> c64_;
    ska::flat_hash_map<u128_t, uint32_t, FHasher> c128_;
    std::vector<int32_t> count_sketch_;
    schism::Schismatic<uint64_t> s64_;
    Counter(size_t cssize=0): ct_(cssize ? COUNTSKETCH_COUNTING: EXACT_COUNTING), count_sketch_(cssize), s64_(cssize + !cssize) {}
    static constexpr uint64_t BM64 = 0x8000000000000000ull;
    void add(u128_t x) {
        if(ct_ == EXACT_COUNTING) {
            auto it = c128_.find(x);
            if(it == c128_.end()) c128_.emplace(x, 1);
            else ++it->second;
        } else {
            const auto hv = sketch::hash::WangHash::hash(uint64_t(x)) ^ sketch::hash::WangHash::hash(x >> 64);
            count_sketch_[s64_.mod(hv)] += ((hv & BM64) ? 1: -1);
        }
    }
    void add(uint64_t x) {
        if(ct_ == EXACT_COUNTING) {
            auto it = c64_.find(x);
            if(it == c64_.end()) c64_.emplace(x, 1);
            else ++it->second;
        } else {
            const auto hv = sketch::hash::WangHash::hash(x);
            count_sketch_[s64_.mod(hv)] += ((hv & BM64) ? 1: -1);
        }
    }
};
struct ParseOptions {

    // K-mer options
    int k_;
    bns::Spacer sp_;
    bns::Encoder<> enc_;
    bns::RollingHasher<uint64_t> rh_;
    bns::RollingHasher<u128_t> rh128_;
    bns::RollingHashingType rht_;
    bool parse_by_seq_ = false;
    bool trim_chr_ = true;
    size_t sketchsize_ = 2048;
    bool one_perm_ = true;

    // Whether to sketch multiset, set, or discrete probability distributions

    SketchSpace sspace_;
    CountingType count_;
    DataType dtype_;
    ParseOptions(int k, int w=-1, std::string spacing="", bns::RollingHashingType rht=bns::DNA, SketchSpace space=SPACE_SET, CountingType count=EXACT_COUNTING, DataType dtype=FASTX):
        k_(k), sp_(k, w > 0 ? w: k, spacing.data()), enc_(sp_), rh_(k), rh128_(k), rht_(rht), sspace_(space), count_(count), dtype_(dtype) {}
    ParseOptions &parse_by_seq() {parse_by_seq_ = true; return *this;}
    ParseOptions &parse_bigwig() {dtype_ = BIGWIG; return *this;}
    ParseOptions &parse_bed() {dtype_ = BED; return *this;}
    ParseOptions &parse_protein() {rh_.enctype_ = rh128_.enctype_ = rht_ = bns::PROTEIN; return *this;}
    ParseOptions &sketchsize(size_t sz) {sketchsize_ = sz; return *this;}
};

using FullSetSketch = sketch::CSetSketch<RegT>;
using OPSetSketch = sketch::setsketch::OPCSetSketch<RegT>;
using BagMinHash = sketch::BagMinHash2<RegT>;

enum OutputKind {
    SYMMETRIC_ALL_PAIRS,
    ASYMMETRIC_ALL_PAIRS,
    KNN_GRAPH,
    NN_GRAPH_THRESHOLD // Variable number of similarities, as given by threshold
};

struct Collection {
    ParseOptions opts_;
    std::vector<std::string> names;
    std::vector<std::vector<uint64_t>> sketch_data;
    std::vector<std::vector<uint64_t>> full_kmers; // Optional
    // For BigWig and BED files, full_kmers will be empty,
    std::vector<std::string> full_strings; // For edit-distance clustering only.
    OutputKind kind_;
    double threshold_; // if threshold_ < 0, use -fraction to keep only items above *fraction* similarity
                       // if threshold_ <= -1, then emit all comparisons
};


std::vector<RegT> bed2sketch(std::string path, const ParseOptions &opts);

#endif
