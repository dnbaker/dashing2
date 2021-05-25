#include "refine.h"
namespace dashing2 {

void refine_results(std::vector<std::vector<PairT>> &lists, Dashing2DistOptions &opts, const SketchingResult &result) {
    //LSHDistType compare(Dashing2DistOptions &opts, const SketchingResult &result, size_t i, size_t j);
    const LSHDistType mult = distance(opts.measure_) ? 1.: -1.;
    // 1. Perform full distance computations over the LSH-selected candidates
    OMP_PFOR_DYN
    for(size_t i = 0; i < lists.size(); ++i) {
        const size_t lhid = i;
        auto beg = lists[i].begin(), e = lists[i].end();
        auto &l = lists[i];
        const size_t lsz = lists[i].size();
        // Skip self for nn graph, we'll add self-connections later
        if(opts.num_neighbors_ > 0 && size_t(opts.num_neighbors_) < lsz) {
            std::fprintf(stderr, "Filtering for top-%d neighbors\n", opts.num_neighbors_);
            for(size_t j = 0; j < lsz; ++j) {
                // -- as above, for the KNN-format
                l[j].first = (l[j].second == lhid)
                ? std::numeric_limits<LSHDistType>::max()
                : mult * compare(opts, result, lhid, l[j].second);
            }
            std::sort(beg, e);
            l.resize(opts.num_neighbors_);
        } else if(opts.min_similarity_ > 0.) {
            std::fprintf(stderr, "Filtering for minimum similarity\n");
            // -- as above, for the NN_GRAPH_THRESHOLD format
            size_t failures = 0;
            for(size_t j = 0; j < lsz; ++j) {
                auto &[dist, id] = l[j];
                if(id == i) {
                    dist = std::numeric_limits<LSHDistType>::max();
                    continue;
                }
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
        if(distance(opts.measure_)) {
            std::transform(beg, e, beg, [&](PairT x) {return PairT{-x.first, x.second};});
        }
    }
}
} // namespace dashing2
