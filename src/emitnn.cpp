#include "emitnn.h"
namespace dashing2 {

// Choose to emit KNN graph in CSR-format
// source = row, dest = column
// Format:
// 16 bytes: uint64_t nids, uint64_t nnz
// (nids + 1) * 8 bytes: indptr in uint64_t
// nnz * sizeof(LSHIDType): indices in LSHIDType (default uint32_t)
// nnz * sizeof(LSHDistType): data in LSHDistType (default float)
void emit_neighbors(std::vector<pqueue> &lists, Dashing2DistOptions &opts, const SketchingResult &result) {
    auto emitstart = std::chrono::high_resolution_clock::now();
    std::string &outname = opts.outfile_path_;
    std::FILE *ofp = stdout;
    if(outname.size() && (ofp = std::fopen(outname.data(), "wb")) == nullptr)
        throw std::runtime_error(std::string("Failed to open file ") + outname + " for writing");
    if(opts.output_format_ == HUMAN_READABLE) {
        std::fprintf(ofp, "#Collection\tNeighbor lists -- name:distance, separated by tabs\n");
        for(size_t i = 0; i < lists.size(); ++i) {
            auto &l = lists[i];
            std::fprintf(ofp, "%s", result.names_[i].data());
            for(size_t j = 0; j < l.size(); ++j) {
                const auto [msr, rhid] = l[j];
                std::fputc('\t', ofp);
                std::fwrite(result.names_[rhid].data(), 1, result.names_[rhid].size(), ofp);
                std::fprintf(ofp, ":%0.8g", msr);
            }
            std::fputc('\n', ofp);
        }
    } else {
        std::vector<uint64_t> indptr(lists.size() + 1);
        for(size_t i = 0; i < lists.size(); ++i) {
            auto &l(lists[i]);
            auto lsz = l.size();
            indptr[i + 1] = indptr[i] + lsz;
        }
        const size_t nnz = indptr.back();
        uint64_t dims [] {lists.size(), nnz};
        std::vector<LSHIDType> indices(nnz);
        checked_fwrite(ofp, dims, sizeof(dims));
        checked_fwrite(ofp, indptr.data(), (indptr.size() * sizeof(uint64_t)));
        for(size_t i = 0; i < lists.size(); ++i)
            for(size_t j = 0; j < lists[i].size(); ++j)
                checked_fwrite(ofp, &lists[i][j].second, sizeof(LSHIDType));
        for(size_t i = 0; i < lists.size(); ++i)
            for(size_t j = 0; j < lists[i].size(); ++j)
                checked_fwrite(ofp, &lists[i][j].first, sizeof(LSHDistType));
    }
    if(ofp != stdout) std::fclose(ofp);
    auto emitstop = std::chrono::high_resolution_clock::now();
    std::fprintf(stderr, "Emitting list took %Lgs.\n", std::chrono::duration<long double, std::ratio<1, 1>>(emitstop - emitstart).count());
}


} // namespace dashing2
