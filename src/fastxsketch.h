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
    std::vector<double> cardinalities_; //apparently holds the cardinalties of a particular input sequence
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
template<typename T, size_t chunk_size = 65536> 
size_t load_copy(const std::string &path, T *ptr, double *cardinality, const size_t ss) {
    T *const origptr = ptr;
    if(path.size() > 3 && std::equal(path.data() + path.size() - 3, &path[path.size()], ".gz")) { //case for .gz files
        gzFile fp = gzopen(path.data(), "rb");
        if(!fp) return 0; //THROW_EXCEPTION(std::runtime_error(std::string("Failed to open file at ") + path));
        gzread(fp, cardinality, sizeof(*cardinality)); //read cardinality into cardinality variable
        for(int nr; 
            !gzeof(fp) && (nr = gzread(fp, ptr, sizeof(T) * chunk_size)) == sizeof(T) * chunk_size;
            ptr += nr / sizeof(T)); //read the registers into 
        gzclose(fp);
        return ptr - origptr;
    } else if(path.size() > 3 && std::equal(path.data() + path.size() - 3, &path[path.size()], ".xz")) { //case for .xz files
        auto cmd = std::string("xz -dc ") + path;
        std::FILE *fp = ::popen(cmd.data(), "r");
        if(fp == 0) return 0;
        std::fread(cardinality, sizeof(*cardinality), 1, fp);
        for(auto up = (uint8_t *)ptr;!std::feof(fp) && std::fread(up, sizeof(T), chunk_size, fp) == chunk_size; up += chunk_size * sizeof(T));
        ::pclose(fp);
        return ptr - origptr;
    }
    std::FILE *fp = bfopen(path.data(), "rb");
    if(!fp) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open ") + path));
    std::fread(cardinality, sizeof(*cardinality), 1, fp);
    const int fd = ::fileno(fp);
    size_t sz = 0;
    if(!::isatty(fd)) {
        struct stat st;
        if(::fstat(fd, &st)) THROW_EXCEPTION(std::runtime_error(std::string("Failed to fstat") + path));
        if(!st.st_size) {
            std::fprintf(stderr, "Warning: Empty file found at %s\n", path.data());
            return 0;
        }
        size_t expected_bytes = st.st_size - 8;
        const size_t expected_sketch_nb = ss * sizeof(T);
        if(expected_bytes != expected_sketch_nb) {
            std::fprintf(stderr, "Expected %zu bytes of sketch, found %zu\n", expected_sketch_nb, expected_bytes);
        }
        size_t nb = std::fread(ptr, 1, expected_bytes, fp);
        if(nb != expected_bytes) {
            std::fprintf(stderr, "Read %zu bytes instead of %zu for file %s\n", nb, expected_bytes, path.data());
            perror("Error reading.");
            THROW_EXCEPTION(std::runtime_error("Error in reading from file"));
        }
        sz = expected_bytes / sizeof(T);
    } else {
        auto up = (uint8_t *)ptr;
        for(;!std::feof(fp) && std::fread(up, sizeof(T), chunk_size, fp) == chunk_size; up += chunk_size * sizeof(T));
        sz = (up - (uint8_t *)ptr) / sizeof(T);
    }
    std::fclose(fp);
    return sz;
}

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
