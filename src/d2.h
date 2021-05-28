#ifndef D2_H_H___
#define D2_H_H___
#ifdef _OPENMP
#include <omp.h>
#endif
#include "enums.h"
#include <memory>
#include <vector>
#include "bonsai/encoder.h"
#include "bigWig.h"
#include "flat_hash_map/flat_hash_map.hpp"
#include "xxHash/xxh3.h"
#include "sketch/setsketch.h"
#include "sketch/bmh.h"
#include "counter.h"
#include "filterset.h"


namespace dashing2 {
using namespace sketch;

// To allow for 64-bit set identifiers, compile with -DLSHIDTYPE=uint64_t
#ifndef LSHIDTYPE
#define LSHIDTYPE uint32_t
#endif
using LSHIDType = LSHIDTYPE;
#undef LSHIDTYPE

#ifndef DASHING2_INDEX_FLOAT_TYPE
#define DASHING2_INDEX_FLOAT_TYPE float
#endif
using LSHDistType = DASHING2_INDEX_FLOAT_TYPE;
#undef DASHING2_INDEX_FLOAT_TYPE


struct IntervalSketchResult {
    using Map = ska::flat_hash_map<std::string, std::vector<RegT>>;
    std::unique_ptr<Map> chrmap_;
    std::unique_ptr<std::vector<RegT>> global_;
};

template<typename F>
void for_each_substr(const F &func, const std::string &s, const int sep=' ') {
    const char *p;
    if((p = std::strchr(s.data(), sep)) == nullptr) {
        func(s);
        return;
    }
    const char *p2 = s.data();
    std::string tmp(p2, p);
    for(;;) {
        func(tmp);
        std::swap(p2, ++p);
        if((p = std::strchr(p2, sep)) == nullptr) {
            tmp = p2;
            func(tmp);
            break;
        }
        tmp = std::string(p2, p);
        if(std::all_of(tmp.begin(), tmp.end(), [](auto x) {return std::isspace(x);})) break;
    }
}

struct Dashing2Options {

    // K-mer options
    int k_, w_;
    bns::Spacer sp_;
    bns::Encoder<> enc_;
    bns::RollingHasher<uint64_t> rh_;
    bns::RollingHasher<u128_t> rh128_;
    bns::RollingHashingType rht_;
    bool parse_by_seq_ = false;
    bool trim_chr_ = true;
    size_t sketchsize_ = 2048;
    double count_threshold_ = 0.;
    bool one_perm_ = true;
    KmerSketchResultType kmer_result_ = ONE_PERM;
    bool by_chrom_ = false;
    bool bed_parse_normalize_intervals_ = false;
    size_t cssize_ = 0;
    bool save_kmers_ = false;
    bool save_kmercounts_ = false;
    bool homopolymer_compress_minimizers_ = false;
    bool trim_folder_paths_ = false; // When writing output files, write to cwd instead of the directory the files came from
    bool cache_sketches_ = false;
    bool build_mmer_matrix_ = false;
    bool build_count_matrix_ = false;
    bool build_sig_matrix_ = true;
    std::string outprefix_;

    // Whether to sketch multiset, set, or discrete probability distributions

    SketchSpace sspace_;
    DataType dtype_;
    bool use128_ = false;
    unsigned nthreads_;

    std::unique_ptr<FilterSet> fs_;
    
    Dashing2Options(int k, int w=-1, bns::RollingHashingType rht=bns::DNA, SketchSpace space=SPACE_SET, DataType dtype=FASTX, size_t nt=0, bool use128=false, std::string spacing=""):
        k_(k), w_(w), sp_(k, w > 0 ? w: k, spacing.data()), enc_(sp_), rh_(k), rh128_(k), rht_(rht), sspace_(space), dtype_(dtype), use128_(use128) {
        std::fprintf(stderr, "Dashing2 made with k = %d, w = %d, space = %s, datatype = %s\n", k, w, to_string(sspace_).data(), to_string(dtype_).data());
        if(nt <= 0) {
            if(char *s = std::getenv("OMP_NUM_THREADS"))
                nt = std::max(std::atoi(s), 1);
            else nt = 1;
        }
        nthreads(nt);
    }
    bool trim_folder_paths() const {
        return trim_folder_paths_ || outprefix_.size();
    }
    auto w() const {return w_;}
    void w(int neww) {w_ = neww; sp_.resize(k_, w_); rh128_.window(neww); rh_.window(neww);}
    Dashing2Options &parse_by_seq(bool v) {parse_by_seq_ = v; return *this;}
    Dashing2Options &parse_bigwig() {dtype_ = BIGWIG; return *this;}
    Dashing2Options &parse_bed() {dtype_ = BED; return *this;}
    Dashing2Options &parse_protein() {rh_.enctype_ = rh128_.enctype_ = rht_ = bns::PROTEIN; return *this;}
    Dashing2Options &sketchsize(size_t sz) {sketchsize_ = sz; return *this;}
    Dashing2Options &use128(size_t sz) {use128_ = sz; return *this;}
    Dashing2Options &cssize(size_t sz) {cssize_ = sz; return *this;}
    Dashing2Options &nthreads(unsigned nthreads) {nthreads_ = nthreads; OMP_ONLY(omp_set_num_threads(nthreads);) return *this;}
    Dashing2Options &kmer_result(KmerSketchResultType kmer_result) {kmer_result_ = kmer_result; return *this;}
    Dashing2Options &cache_sketches(bool v) {cache_sketches_ = v;return *this;}
    Dashing2Options &save_kmers(bool v) {save_kmers_ = v;return *this;}
    Dashing2Options &save_kmercounts(bool v) {save_kmercounts_ = v;return *this;}
    Dashing2Options &outprefix(const std::string &v) {outprefix_ = v;return *this;}
    bool parse_by_seq() {return parse_by_seq_;}
    std::string &outprefix() {return outprefix_;}
    const std::string &outprefix() const {return outprefix_;}
    bool cache_sketches() const {return cache_sketches_;}
    unsigned nthreads() const {return nthreads_;}
    unsigned kmer_result() const {return kmer_result_;}
    size_t sketchsize() const {return sketchsize_;}
    bool use128() const {return use128_;}
    size_t cssize() const {return cssize_;}
    CountingType ct() const {return cssize_ > 0 ? COUNTSKETCH_COUNTING: EXACT_COUNTING;}
    CountingType count() const {return ct();}
};




using FullSetSketch = sketch::CSetSketch<RegT>;
using OPSetSketch = sketch::setsketch::OPCSetSketch<RegT>;
using BagMinHash = sketch::BagMinHash2<RegT>;
using ProbMinHash = sketch::pmh2_t<RegT>;
using OrderMinHash = sketch::omh::OMHasher<RegT>;

} // namespace dashing2
//std::vector<RegT> reduce(ska::flat_hash_map<std::string, std::vector<RegT>> &map);

#endif
