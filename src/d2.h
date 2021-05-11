#ifndef D2_H_H___
#define D2_H_H___
#ifdef _OPENMP
#include <omp.h>
#endif
#include <memory>
#include <vector>
#include "bonsai/encoder.h"
#include "bigWig.h"
#include "flat_hash_map/flat_hash_map.hpp"
#include "xxHash/xxh3.h"
#include "sketch/setsketch.h"
#include "sketch/bmh.h"
#include "enums.h"
#include "counter.h"


namespace dashing2 {
using namespace sketch;

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
    size_t countsketchsize_ = 80'000'000ull; // Default to 80-million count-sketch
    bool one_perm_ = true;
    bool by_chrom_ = false;
    bool bed_parse_normalize_intervals_ = false;
    size_t cssize_ = 0;
    bool save_kmers_ = false;
    bool save_kmercounts_ = false;

    // Whether to sketch multiset, set, or discrete probability distributions
 
    SketchSpace sspace_;
    CountingType count_;
    DataType dtype_;
    bool use128_ = false;
    size_t nthreads_;
    ParseOptions(int k, int w=-1, std::string spacing="", bns::RollingHashingType rht=bns::DNA, SketchSpace space=SPACE_SET, CountingType count=EXACT_COUNTING, DataType dtype=FASTX, size_t nthreads=0, bool use128=false):
        k_(k), sp_(k, w > 0 ? w: k, spacing.data()), enc_(sp_), rh_(k), rh128_(k), rht_(rht), sspace_(space), count_(count), dtype_(dtype), use128_(use128) {
        if(nthreads <= 0) {
            if(char *s = std::getenv("OMP_NUM_THREADS"))
                nthreads = std::max(std::atoi(s), 1);
        } else nthreads = 1;
        nthreads_ = nthreads;
    }
    ParseOptions &parse_by_seq() {parse_by_seq_ = true; return *this;}
    ParseOptions &parse_bigwig() {dtype_ = BIGWIG; return *this;}
    ParseOptions &parse_bed() {dtype_ = BED; return *this;}
    ParseOptions &parse_protein() {rh_.enctype_ = rh128_.enctype_ = rht_ = bns::PROTEIN; return *this;}
    ParseOptions &sketchsize(size_t sz) {sketchsize_ = sz; return *this;}
    ParseOptions &use128(size_t sz) {use128_ = sz; return *this;}
    ParseOptions &cssize(size_t sz) {cssize_ = sz; return *this;}
    ParseOptions &nthreads(size_t nthreads) {nthreads_ = nthreads; OMP_ONLY(omp_set_num_threads(nthreads);) return *this;}
    size_t nthreads() const {return nthreads_;}
    size_t sketchsize() const {return sketchsize_;}
    size_t use128() const {return use128_;}
    size_t cssize() const {return cssize_;}
};




using FullSetSketch = sketch::CSetSketch<RegT>;
using OPSetSketch = sketch::setsketch::OPCSetSketch<RegT>;
using BagMinHash = sketch::BagMinHash2<RegT>;
using ProbMinHash = sketch::pmh2_t<RegT>;
using OrderMinHash = sketch::omh::OMHasher<RegT>;

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
}
//std::vector<RegT> reduce(ska::flat_hash_map<std::string, std::vector<RegT>> &map);

#endif
