#ifndef D2_CMP_H__
#define D2_CMP_H__
#include "d2.h"
#include "fastxsketch.h"

namespace dashing2 {

enum Measure {
    SIMILARITY,  // Jaccard or set similarity. Default behavior
    CONTAINMENT, // (A & B) / A
    SYMMETRIC_CONTAINMENT, // (A & B) / min(|A|, |B|)
    POISSON_LLR, // -log-transformed similarity
    INTERSECTION, // |A & B|
    UNION_SIZE,   // |A | B|
    M_EDIT_DISTANCE, // Whether in minimizer-sequence space or sequence edit distance space.
    MASH_DISTANCE = POISSON_LLR
};
static inline std::string to_string(Measure m) {
    if(m == SIMILARITY) return "SIMILARITY";
    if(m == CONTAINMENT) return "CONTAINMENT";
    if(m == SYMMETRIC_CONTAINMENT) return "SYMMETRIC_CONTAINMENT";
    if(m == POISSON_LLR) return "POISSON_LLR";
    if(m == INTERSECTION) return "INTERSECTION";
    if(m == UNION_SIZE) return "UNION_SIZE";
    if(m == M_EDIT_DISTANCE) return "EDIT_DISTANCE";
    return "UNKNOWN";
}

struct SketchingResult;

template<Measure measure>
struct is_symmetric: public std::true_type {};
template<> struct is_symmetric<CONTAINMENT>: public std::false_type {};
#if __cplusplus >= 201703L
template<Measure measure>
static constexpr bool is_symmetric_v = is_symmetric<measure>::value;
#endif
static constexpr inline bool symmetric(Measure msr) {
    switch(msr) {
        case CONTAINMENT: return false;
        default: return true;
    }
}
static constexpr inline bool distance(Measure msr) {
    switch(msr) {
        case UNION_SIZE: case INTERSECTION: case SIMILARITY: case CONTAINMENT: return false;
        default: return true;
    }
}


struct Dashing2DistOptions: public Dashing2Options {
    OutputKind output_kind_;
    OutputFormat output_format_;
    // If <= 0, uses setsketch to compress before distances
    // If 1 (> 0), generates b-bit signatures and truncates
    int num_neighbors_ = -1; // Only emits top-"nn" neighbors
    double min_similarity_ = -1.; // Only emit similarities which are above min_similarity_ if nonnegative
    mutable void *compressed_ptr_ = nullptr;
    mutable Measure measure_ = SIMILARITY;
    std::string outfile_path_;
    mutable bool exact_kmer_dist_ = false;
    bool refine_exact_ = false;
    size_t cmp_batch_size_ = 16;
    unsigned int nLSH = 2;
    Dashing2DistOptions(Dashing2Options &opts, OutputKind outres, OutputFormat of, double nbytes_for_fastdists=-1, int truncate_method=0, int nneighbors=-1, double minsim=-1., std::string outpath="/workspaces/BMI_Project/lib_integration_tests/dashing2/distout.txt", bool exact_kmer_dist=false, bool refine_exact=false, int nlshsubs=3): //modified initialization of outpath to not be empty -> this should just be temporary
        Dashing2Options(opts), output_kind_(outres), output_format_(of), outfile_path_(outpath), exact_kmer_dist_(exact_kmer_dist), refine_exact_(refine_exact), nLSH(nlshsubs)
    {
        if(verbosity) {
            std::fprintf(stderr, "[%s] output format is %s\n", __PRETTY_FUNCTION__, ::dashing2::to_string(output_format_).data());
        }
        set_sketch_compressed();
        if(nbytes_for_fastdists < 0) nbytes_for_fastdists = sizeof(RegT);
        if(std::fmod(nbytes_for_fastdists, 1.)) {
            if(nbytes_for_fastdists != 0.5) THROW_EXCEPTION(std::runtime_error("Can only do 1, 2, 4, 8, or 0.5 bytes per register"));
        } else if(int(nbytes_for_fastdists) & (int(nbytes_for_fastdists - 1)))
            THROW_EXCEPTION(std::runtime_error("Can't compress sketches to non-power of 2 register size. Should be < sizeof(RegT)"));
        fd_level_ = nbytes_for_fastdists;
        truncation_method_ = truncate_method;
        num_neighbors_ = nneighbors;
        min_similarity_ = minsim;
        if(this->kmer_result_ >= FULL_MMER_SET)
            exact_kmer_dist_ = true;
        if(outfile_path_.empty() || outfile_path_ == "-") outfile_path_ = "/dev/stdout";
        if(nLSH < 1) nLSH = 1;
        if(this->sketch_compressed_set) {
            if(size_t rem = opts.sketchsize_ % sizeof(RegT) / opts.fd_level_; rem != 0) {
                opts.sketchsize_ += sizeof(RegT) / opts.fd_level_ - rem;
                std::fprintf(stderr, "Sketchsize is not %zu-bit register multiple; padding the number of registers to fit. New number of registers: %zu\n", sizeof(RegT) * 8, opts.sketchsize_);
            }
            const size_t mul = 8 / nbytes_for_fastdists;
            if(const size_t rem = opts.sketchsize_ % size_t(8 / nbytes_for_fastdists); rem) {
                std::fprintf(stderr, "When sketching compressed, always pad to 64-bit sets.\n");
                opts.sketchsize_ += mul - rem;
            }
        }
        validate();
    }
    bool truncate_mode() const {return truncation_method_;}
    void validate() const {
        Dashing2Options::validate();
        if(num_neighbors_ > 0 && min_similarity_ > 0.) {
            THROW_EXCEPTION(std::invalid_argument("invalid: nn > 0 and minsim > 0. Pick either top-k or minimum similarity. (Can't do both.)"));
        }
        if((sspace_ == SPACE_PSET || sspace_ == SPACE_EDIT_DISTANCE) && measure_ == INTERSECTION) {
            std::fprintf(stderr, "Can't estimate intersection sizes for ProbMinHash due to the implicit normalization. Reverting to similarity\n");
            measure_ = SIMILARITY;
        }
        if((sspace_ == SPACE_PSET || sspace_ == SPACE_EDIT_DISTANCE) && measure_ == UNION_SIZE) {
            std::fprintf(stderr, "Can't estimate union sizes for ProbMinHash due to the implicit normalization. Reverting to similarity\n");
            measure_ = SIMILARITY;
        }
        if((sspace_ == SPACE_EDIT_DISTANCE) && measure_ != SIMILARITY && measure_ != M_EDIT_DISTANCE) {
            std::fprintf(stderr, "Edit distance LSH, but measure is not similarity. Switching to edit distance so that it is well defined\n");
             measure_ = M_EDIT_DISTANCE;
        }
        if(sketch_compressed_set) {
            if(this->kmer_result_ != FULL_SETSKETCH)
                THROW_EXCEPTION(std::invalid_argument("Sketch compressed is only available for FullSetSketch."));
            if(this->compressed_b_ < 1.L) THROW_EXCEPTION(std::invalid_argument("base must be >= 1."));
            if(this->compressed_a_ <= 0.L) THROW_EXCEPTION(std::invalid_argument("offset a must be > 0."));
        }
#ifdef __aarch64__
        if((truncation_method_ <= 0) && (fd_level_ < 1.)) {
            std::fprintf(stderr, "Warning: on apple silicon, accuracy may be reduced by the lack of 80-bit floating-point support. Consider using --fastcmp/--regsize 1 or greater.\n");
        }
#endif
    }
};
void cmp_core(const Dashing2DistOptions &ddo, SketchingResult &res);
LSHDistType compare(const Dashing2DistOptions &opts, const SketchingResult &result, size_t i, size_t j);
void emit_rectangular(const Dashing2DistOptions &opts, const SketchingResult &result);
size_t default_batchsize(size_t &batch_size, const Dashing2DistOptions &opts);


}

#endif
