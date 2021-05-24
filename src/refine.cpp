#include "refine.h"
namespace dashing2 {

void refine_results(std::vector<std::vector<PairT>> &lists, Dashing2DistOptions &opts, const SketchingResult &result) {
    LSHDistType compare(Dashing2DistOptions &opts, const SketchingResult &result, size_t i, size_t j);
    const LSHDistType mult = distance(opts.measure_) ? 1.: -1.;
    // 1. Perform full distance computations over the LSH-selected candidates
    OMP_PFOR_DYN
    for(size_t i = 0; i < lists.size(); ++i) {
        const size_t lhid = i;
        auto beg = lists[i].begin(), e = lists[i].end();
        auto &l = lists[i];
        const size_t lsz = lists[i].size();
        if(opts.num_neighbors_ > 0 && size_t(opts.num_neighbors_) < lsz) {
            for(size_t j = 0; j < lsz; ++j) {
                l[j].first = mult * compare(opts, result, lhid, x.second);
            }
            std::sort(beg, e);
            l.resize(opts.num_neighbors_);
        } else if(opts.min_similarity_ > 0.) {
            std::fprintf(stderr, "Filtering for minimum similarity\n");
            size_t failures = 0;
            for(size_t j = 0; j < lsz; ++j) {
                auto v = compare(opts, result, lhid, l[j].second);
                const bool pass = distance(opts.measure_) ? v < opts.min_similarity_: v >= opts.min_similarity_;
                if(!pass) {
                    l[j] = {std::numeric_limits<LSHDistType>::max(), LSHIDType(0)};
                    if(++failures == 5) {
                            l.resize(j);
                            break;
                        }
                    }
                } else {
                    l[j].second = v * mult;
                    failures = 0;
                }
            }
            std::sort(l.begin(), l.end());
            l.erase(std::find_if(l.begin(), l.end(), [dist=distance(opts.measure_),ms=opts.min_similarity_](auto x) {
                return dist ? x.first < ms: -x.first >= ms;
            }), l.end());
            std::fprintf(stderr, "Smallest similarity: %g. threshold: %g\n", mult * l.back().first, opts.min_similarity_);
        } else {
            std::transform(beg, e, beg, [&](auto x) -> PairT {return {mult * compare(opts, result, lhid, x.second), x.second};});
            std::sort(beg, e);
        }
    }
}
} // namespace dashing2
