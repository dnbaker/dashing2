#ifndef D2_CMP_H__
#define D2_CMP_H__
#include "d2.h"
#include "fastxsketch.h"

namespace dashing2 {

#if 0
enum OutputKind {
    SYMMETRIC_ALL_PAIRS,
    ASYMMETRIC_ALL_PAIRS,
    KNN_GRAPH, // Fixed top-k neighbors
    NN_GRAPH_THRESHOLD, // Variable number of similarities, as given by threshold
    DEDUP
}
enum OutputFormat {
    MACHINE_READABLE,
    HUMAN_READABLE,
    BINARY = MACHINE_READABLE
};
#endif
enum Measure {
    SIMILARITY,  // Jaccard or set similarity. Default behavior
    CONTAINMENT, // (A & B) / A
    SYMMETRIC_CONTAINMENT, // (A & B) / min(|A|, |B|)
    POISSON_LLR, // -log-transformed similarity
    INTERSECTION, // |A & B|
    M_EDIT_DISTANCE = POISSON_LLR, // Whether in minimizer-sequence space or sequence edit distance space.
    MASH_DISTANCE = POISSON_LLR
};
static inline std::string to_string(Measure m) {
    if(m == SIMILARITY) return "SIMILARITY";
    if(m == CONTAINMENT) return "CONTAINMENT";
    if(m == SYMMETRIC_CONTAINMENT) return "SYMMETRIC_CONTAINMENT";
    if(m == POISSON_LLR) return "POISSON_LLR";
    if(m == INTERSECTION) return "INTERSECTION";
    return "UNKNOWN";
}

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
        case INTERSECTION: case SIMILARITY: case CONTAINMENT: return false;
        default: return true;
    }
}


struct Dashing2DistOptions: public Dashing2Options {
    OutputKind output_kind_;
    OutputFormat output_format_;
    double fd_level_; // Determines the number of bytes to which the minhash sketches are compressed
    int truncation_method_ = 0;
    // If <= 0, uses setsketch to compress before distances
    // If 1 (> 0), generates b-bit signatures and truncates
    int num_neighbors_ = -1; // Only emits top-"nn" neighbors
    double min_similarity_ = -1.; // Only emit similarities which are above min_similarity_ if nonnegative
    void *compressed_ptr_ = nullptr;
    double compressed_b_ = -1.;
    double compressed_a_ = -1.;
    mutable Measure measure_ = SIMILARITY;
    std::string outfile_path_;
    bool exact_kmer_dist_ = false;
    bool refine_exact_ = false;
    Dashing2DistOptions(Dashing2Options &opts, OutputKind outres, OutputFormat of, double nbytes_for_fastdists=-1, int truncate_method=0, int nneighbors=-1, double minsim=-1., std::string outpath="", bool exact_kmer_dist=false, bool refine_exact=false): Dashing2Options(std::move(opts)), output_kind_(outres), output_format_(of), outfile_path_(outpath), exact_kmer_dist_(exact_kmer_dist), refine_exact_(refine_exact)
    {
        if(nbytes_for_fastdists < 0) nbytes_for_fastdists = sizeof(RegT);
        if(std::fmod(nbytes_for_fastdists, 1.)) {
            if(nbytes_for_fastdists != 0.5) throw std::runtime_error("Can only do 1, 2, 4, 8, or 0.5 bytes per register");
        } else if(int(nbytes_for_fastdists) & (int(nbytes_for_fastdists - 1)))
            throw std::runtime_error("Can't compress sketches to non-power of 2 register size. Should be < sizeof(RegT)");
        fd_level_ = nbytes_for_fastdists;
        truncation_method_ = truncate_method;
        num_neighbors_ = nneighbors;
        min_similarity_ = minsim;
        if(this->kmer_result_ >= FULL_MMER_SET) {
            exact_kmer_dist_ = true;
        }
        validate();
    }
    ~Dashing2DistOptions() {
        std::free(compressed_ptr_);
    }
    void validate() const {
        Dashing2Options::validate();
        if(num_neighbors_ > 0 && min_similarity_ > 0.) {
            throw std::invalid_argument("invalid: nn > 0 and minsim > 0. Pick either top-k or minimum similarity. (Can't do both.)");
        }
        if((sspace_ == SPACE_PSET || sspace_ == SPACE_EDIT_DISTANCE) && measure_ == INTERSECTION) {
            std::fprintf(stderr, "Can't estimate intersection sizes for ProbMinHash due to the implicit normalization. Reverting to similarity\n");
            measure_ = SIMILARITY;
        }
        if((sspace_ == SPACE_EDIT_DISTANCE) && measure_ != SIMILARITY && measure_ != M_EDIT_DISTANCE) {
            std::fprintf(stderr, "Edit distance LSH, but measure is not similarity. Switching to edit distance so that it is well defined\n");
             measure_ = M_EDIT_DISTANCE;
        }
    }
};
void cmp_core(Dashing2DistOptions &ddo, SketchingResult &res);
LSHDistType compare(Dashing2DistOptions &opts, const SketchingResult &result, size_t i, size_t j);

}

#endif
