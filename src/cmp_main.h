#ifndef D2_CMP_H__
#define D2_CMP_H__
#include "d2.h"

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
struct Dashing2DistOptions: public Dashing2Options {
    OutputKind output_kind_;
    OutputFormat output_format_;
    int fd_level_; // Determines the number of bytes to which the minhash sketches are compressed
    int truncation_method_ = 0;
    // If 0, uses setsketch to compress before distances
    // If 1, generates b-bit signatures and truncates
    int num_neighbors_ = -1; // Only emits top-"nn" neighbors
    double min_similarity_ = -1.; // Only emit similarities which are above min_similarity_ if nonnegative
    Dashing2DistOptions(Dashing2Options &opts, OutputKind outres, OutputFormat of, size_t nbytes_for_fastdists=-1, bool truncate_bbits=false, int nneighbors=-1., double minsim=-1.): Dashing2Options(opts), output_kind_(outres), output_format_(of) {
        if(nbytes_for_fastdists == -1) nbytes_for_fastdists == sizeof(RegT);
        switch(nbytes_for_fastdists) {
            default: throw std::runtime_error("Can't compress sketches to non-power of 2 register size. Should be < sizeof(RegT)"); break;
            case 1: case 2: case 4: case 8: case 16: break;
        }
        fd_level_ = nbytes_for_fastdists;
        if(truncate_bbits) truncation_method_ = 1;
        num_neighbors_ = nneighbors;
        min_similarity_ = minsim;
    }
};

}

#endif
