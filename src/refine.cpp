#include "refine.h"
namespace dashing2 {

static constexpr LSHDistType MDIST = std::numeric_limits<LSHDistType>::max();

void refine_results(std::vector<pqueue> &lists, const Dashing2DistOptions &opts, const SketchingResult &result) {
    //LSHDistType compare(Dashing2DistOptions &opts, const SketchingResult &result, size_t i, size_t j);
    const LSHDistType mult = distance(opts.measure_) ? 1.: -1.;
    // 1. Perform full distance computations over the LSH-selected candidates
    if(opts.refine_exact_ && !opts.exact_kmer_dist_) {
        if(opts.kmer_result_ <= FULL_SETSKETCH && opts.compressed_ptr_) {
            opts.compressed_ptr_ = 0;
        } else {
            opts.exact_kmer_dist_ = true;
        }
    }

    auto refinstart = std::chrono::high_resolution_clock::now();
    OMP_PFOR_DYN
    for(size_t i = 0; i < lists.size(); ++i) {
        const size_t lhid = i;
        auto &l = lists[i];
        // Selves are not in the list; We'll add self-connections later
        auto beg = l.begin(), e = l.end();
        const size_t lsz = l.size();
        DBG_ONLY(std::fprintf(stderr, "Processing seqset %zu/%s\n", i, result.names_[i].data());)
        if(opts.num_neighbors_ > 0 && size_t(opts.num_neighbors_) < lsz) {
            // -- as above, for the KNN-format
#ifndef NDEBUG
            std::fprintf(stderr, "lsz: %zu. nn: %u\n", lsz, opts.num_neighbors_);
            std::fprintf(stderr, "Filtering for top-%d neighbors\n", opts.num_neighbors_);
#endif
            for(size_t j = 0; j < lsz; ++j) {
                auto &[dist, id] = l[j];
                dist = mult * compare(opts, result, lhid, id);
            }
            std::sort(beg, e);
            if(size_t(opts.num_neighbors_) < l.size() - 1)
                l.resize(opts.num_neighbors_);
        } else if(opts.min_similarity_ > 0.) {
            static constexpr size_t EARLY_FAILURE_EXIT_THRESHOLD = 25u;
            // -- as above, for the NN_GRAPH_THRESHOLD format
            // This stopping after `EARLY_FAILURE_EXIT_THRESHOLD` consecutive beyond-threshold points is purely heuristic
            // and may change.
            size_t failures = 0;
            for(size_t j = 0; j < lsz; ++j) {
                auto &[dist, id] = l[j];
                auto v = compare(opts, result, lhid, id);
                const bool pass = distance(opts.measure_) ? v < opts.min_similarity_: v >= opts.min_similarity_;
                if(!pass) {
                    dist = MDIST;
                    if(++failures == EARLY_FAILURE_EXIT_THRESHOLD) {
                        l.resize(j);
                        break;
                    }
                } else {
                    dist = v * mult;
                    failures = 0;
                }
            }
            l.erase(std::remove_if(l.begin(), l.end(), [dist=distance(opts.measure_),ms=opts.min_similarity_](const PairT x) -> bool {
                return x.first == MDIST || (dist ? x.first > ms: -x.first < ms);
            }), l.end());
            std::sort(l.begin(), l.end());
#ifndef NDEBUG
            if(!l.empty()) {
                std::fprintf(stderr, "Smallest comparison: %g. threshold: %g\n", l.back().first, opts.min_similarity_);
            }
#endif
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
    auto refinstop = std::chrono::high_resolution_clock::now();

    std::fprintf(stderr, "List refinement took %Lgs.\n", std::chrono::duration<long double, std::ratio<1, 1>>(refinstop - refinstart).count());
}
} // namespace dashing2
