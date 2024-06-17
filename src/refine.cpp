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
        if(opts.num_neighbors_ > 0) {
            // -- as above, for the KNN-format
            for(size_t j = 0; j < lsz; ++j) {
                auto &[dist, id] = l[j];
                dist = mult * compare(opts, result, lhid, id, nullptr);
            }
            std::sort(beg, e);
            if(!distance(opts.measure_)) {
                //Trimming all neighbors with 0 similarity.
                l.erase(std::find_if(beg, e, [](const auto &x) {return x.first == 0.;}), e);
                beg = l.begin(), e = l.end();
            }
            if(size_t(opts.num_neighbors_) < l.size()) {
                l.erase(std::find_if(beg + opts.num_neighbors_, e, [bs=l[opts.num_neighbors_ - 1].first](const auto &x) {return x.first > bs;}),
                        e);
            }

        } else if(opts.min_similarity_ > 0.) {
            static constexpr size_t EARLY_FAILURE_EXIT_THRESHOLD = 20u;
            // -- as above, for the NN_GRAPH_THRESHOLD format
            // This stopping after `EARLY_FAILURE_EXIT_THRESHOLD` consecutive beyond-threshold points is purely heuristic
            // and may change.
            size_t failures = 0;
            for(size_t j = 0; j < lsz; ++j) {
                auto &[dist, id] = l[j];
                auto v = compare(opts, result, lhid, id, nullptr);
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
        } else {
            std::transform(beg, e, beg, [&](PairT x) -> PairT {return {mult * compare(opts, result, lhid, x.second, nullptr), x.second};});
            std::sort(beg, e);
        }
        // Now that we've selected the top-k/bottom-k (similarity/distance), multiply
        if(!distance(opts.measure_)) {
            std::transform(beg, e, beg, [&](PairT x) {return PairT{-x.first, x.second};});
        }
    }
    auto refinstop = std::chrono::high_resolution_clock::now();

    std::fprintf(stderr, "List refinement took %Lgs.\n", std::chrono::duration<long double, std::ratio<1, 1>>(refinstop - refinstart).count());
}
} // namespace dashing2
