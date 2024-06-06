#pragma once
#ifndef DASHING2_FASTX_SKETCH_H__
#define DASHING2_FASTX_SKETCH_H__
#include "d2.h"
#include "mmvec.h"
#include "tmpseqs.h"
#include "cmp_main.h"
#include <variant>

namespace dashing2 {
using std::to_string;
template<typename T>
static inline std::string to_string(const T *ptr) {
    std::ostringstream oss;
    oss << static_cast<const void *>(ptr);
    return oss.str();
}

static bool seqs_in_memory = false;

struct Dashing2DistOptions;

struct SketchingResult {
    SketchingResult(): sequences_(seqs_in_memory) {}
    SketchingResult(SketchingResult &&o) = default;
    SketchingResult(const SketchingResult &o) = delete;
    SketchingResult(SketchingResult &o) = delete;
    SketchingResult &operator=(SketchingResult &&o) = default;
    SketchingResult &operator=(const SketchingResult &o) = delete;

    std::vector<std::string> names_; // List of files, potentially multiple per line
    std::vector<std::string> destination_files_; // Contains sketches/kmer-sets,kmer-sequences etc.
    std::vector<std::string> kmerfiles_;         // Contains file-paths for k-mers, if saved.
    std::vector<std::string> kmercountfiles_;    // Contains k-mer counts, if saved
    // kmerfiles and kmercountfiles are unset for bed files
    std::vector<uint32_t> nperfile_; // This is either empty (in which case each filename/row has its own sketch)
                                     // Or, this contains a list indicating the number of sketches created for each file/line
    std::vector<double> cardinalities_;
    tmpseq::MemoryOrRAMSequences sequences_;
    // This is only filled if sspace is SPACE_EDIT_DISTANCE and
    // we are using LSH only for pre-filtering but performing exact distance calculations via edit distance
    mm::vector<RegT> signatures_; //holds register values of sketches
    // TODO: mmap these matrices to reduce peak memory footprint
    mm::vector<uint64_t> kmers_;
    std::vector<float> kmercounts_; // Contains counts for k-mers, if desired
    // This contains the k-mers corresponding to signatures, if asked for 128-bit k-mers, these are stored in chunks of 2 64-bit integers.
    size_t nq = 0;
    size_t total_seqs() const {
        // Sum of nperfile if nonempty
        // otherwise, just one sequence/bag of k-mers per "name"
        return nperfile_.size() ? std::accumulate(nperfile_.begin(), nperfile_.end(), size_t(0)): names_.size();
    }
    std::string str() const;
    static SketchingResult merge(SketchingResult *start, size_t n, const std::vector<std::string> &);
    void print();
    size_t nqueries() const {return nq;}
    void nqueries(size_t nqnew) {nq = nqnew;}
};
using FastxSketchingResult = SketchingResult;
void seq_resize(std::vector<std::string>& seqs, const size_t num_seqs);
void seq_resize(tmpseq::Seqs&seqs, const size_t num_seqs) noexcept;
void seq_resize(tmpseq::MemoryOrRAMSequences& , const size_t num_seqs);
int32_t num_threads();

//Added declaration of load_copy from fastxsketch.cpp so that I can call it in my own project to load sketches from files into SketchingResult objects
template<typename T, size_t chunk_size = 65536> //this makes load copy callable with all points of types for T, e.g. doubles, ints etc. 
size_t load_copy(const std::string &path, T *ptr, double *cardinality, const size_t ss);


FastxSketchingResult &fastx2sketch(FastxSketchingResult &res, Dashing2Options &opts, const std::vector<std::string> &paths, std::string path);
FastxSketchingResult &fastx2sketch_byseq(FastxSketchingResult &res, Dashing2DistOptions &opts, const std::string &path, kseq_t *kseqs, std::string outpath, bool parallel=false, const size_t seqs_per_batch = 8192);
std::string makedest(Dashing2Options &opts, const std::string &path, bool iskmer=false);

namespace variation {
using ByteSetS = sketch::setsketch::CFByteSetS;
using NibbleSetS = sketch::setsketch::CFNibbleSetS;
using ShortSetS = sketch::setsketch::CFShortSetS;
using UintSetS = sketch::setsketch::CFUintSetS;
using VSetSketch = std::variant<NibbleSetS, ByteSetS, ShortSetS, UintSetS>;

INLINE const RegT *getdata(VSetSketch &o) {
    const RegT *ret;
    std::visit([&ret](auto &x) {ret = (const RegT *)x.data();}, o);
    return ret;
}
INLINE double getcard(VSetSketch &o) {
    double ret;
    std::visit([&ret](auto &x) {ret = x.getcard();}, o);
    return ret;
}
} // namespace variation
using variation::VSetSketch;
} // namespace dashing2

#endif
