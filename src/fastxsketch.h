#pragma once
#include "d2.h"

namespace dashing2 {
using std::to_string;
template<typename T>
static inline std::string to_string(const T *ptr) {
    char buf_[64];
    return std::string(buf_, std::sprintf(buf_, "%p", static_cast<const void *>(ptr)));
}
struct SketchingResult {
    std::vector<std::string> names_; // List of files, potentially multiple per line
    std::vector<std::string> destination_files_; // Contains sketches/kmer-sets,kmer-sequences etc.
    std::vector<std::string> kmerfiles_;         // Contains file-paths for k-mers, if saved.
    std::vector<std::string> kmercountfiles_;    // Contains k-mer counts, if saved
    // kmerfiles and kmercountfiles are unset for bed files
    std::vector<uint32_t> nperfile_; // This is either empty (in which case each filename/row has its own sketch)
                                     // Or, this contains a list indicating the number of sketches created for each file/line
    std::vector<RegT> signatures_; // Signatures, packed into a single array
    std::vector<uint64_t> kmers_;  // This contains the k-mers corresponding to signatures, if asked for
    std::vector<uint32_t> kmercounts_; // Contains counts for k-mers, if desired
    const ParseOptions *options_ = nullptr;
    size_t total_seqs() const {
        // Sum of nperfile if nonempty
        // otherwise, just one sequence/bag of k-mers per "name"
        return nperfile_.size() ? std::accumulate(nperfile_.begin(), nperfile_.end(), size_t(0)): names_.size();
    }
    std::string str() const;
    void print();
};
using FastxSketchingResult = SketchingResult;

FastxSketchingResult fastx2sketch(ParseOptions &opts, std::vector<std::string> &paths);
}
