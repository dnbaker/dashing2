#include "refine.h"
namespace dashing2 {

void refine_results(std::vector<pqueue> &lists, Dashing2DistOptions &opts, const SketchingResult &result) {
    //LSHDistType compare(Dashing2DistOptions &opts, const SketchingResult &result, size_t i, size_t j);
    const LSHDistType mult = distance(opts.measure_) ? 1.: -1.;
    // 1. Perform full distance computations over the LSH-selected candidates
    OMP_PFOR_DYN
    for(size_t i = 0; i < lists.size(); ++i) {
        const size_t lhid = i;
        auto &l = lists[i];
        // Remove self from the list
        // We'll add self-connections later
        if(auto it = std::find_if(l.begin(), l.end(), [lhid](auto x) {return x.second == lhid;}); it != l.end()) {
            std::fprintf(stderr, "Found self unexpectedly in list: %zu at %zd\n", i, it - l.begin());
        }
        //l.erase(std::remove_if(l.begin(), l.end(), [](auto x) {return x.second == i;}));
        auto beg = l.begin(), e = l.end();
        const size_t lsz = l.size();
        std::fprintf(stderr, "Processing seqset %zu/%s\n", i, result.names_[i].data());
        if(opts.num_neighbors_ > 0 && size_t(opts.num_neighbors_) < lsz) {
            // -- as above, for the KNN-format
            std::fprintf(stderr, "lsz: %zu. nn: %u\n", lsz, opts.num_neighbors_);
            std::fprintf(stderr, "Filtering for top-%d neighbors\n", opts.num_neighbors_);
            for(size_t j = 0; j < lsz; ++j) {
                auto &[dist, id] = l[j];
                dist = mult * compare(opts, result, lhid, id);
            }
            std::sort(beg, e);
            if(size_t(opts.num_neighbors_) < l.size() - 1)
                l.resize(opts.num_neighbors_);
        } else if(opts.min_similarity_ > 0.) {
            // -- as above, for the NN_GRAPH_THRESHOLD format
            std::fprintf(stderr, "Filtering for minimum similarity\n");
            // This stopping after 5 consecutive beyond-threshold points is purely heuristic
            // and may change.
            size_t failures = 0;
            for(size_t j = 0; j < lsz; ++j) {
                auto &[dist, id] = l[j];
                auto v = compare(opts, result, lhid, id);
                const bool pass = distance(opts.measure_) ? v < opts.min_similarity_: v >= opts.min_similarity_;
                if(!pass) {
                    dist = std::numeric_limits<LSHDistType>::max();
                    if(++failures == 5) {
                        l.resize(j);
                        break;
                    }
                } else {
                    dist = v * mult;
                    failures = 0;
                }
            }
            std::sort(l.begin(), l.end());
            l.erase(std::find_if(l.begin(), l.end(), [dist=distance(opts.measure_),ms=opts.min_similarity_](PairT x) {
                return dist ? x.first < ms: -x.first >= ms;
            }), l.end());
            std::fprintf(stderr, "Smallest similarity: %g. threshold: %g\n", mult * l.back().first, opts.min_similarity_);
        } else {
            std::transform(beg, e, beg, [&](PairT x) -> PairT {return {mult * compare(opts, result, lhid, x.second), x.second};});
            std::sort(beg, e);
        }
        // Now that we've selected the top-k/bottom-k (similarity/distance), multiply
        if(!distance(opts.measure_)) {
            std::transform(beg, e, beg, [&](PairT x) {return PairT{-x.first, x.second};});
            //std::reverse(beg, e);
        }
    }
}
} // namespace dashing2
