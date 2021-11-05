#include "cmp_main.h"
#include "minispan.h"
#include <vector>
#include <mutex>
#include "index_build.h"
#include "dedup_core.h"

namespace dashing2 {


// std::less -> vector.back() == heap.top(), where it is the smallest of the top
// if std::less<void>{}(vector.back(), newitem), then pop heap, add item and push heap
// We store similarities as -x in our nearest neighbor lists so that
// we only compile one heap update step -- for smallest distances
void pqueue::erase(typename std::priority_queue<PairT>::container_type::iterator it, typename std::priority_queue<PairT>::container_type::iterator oit) {
    this->c.erase(it, oit);
}


void update(pqueue &x, flat_hash_set<LSHIDType> &xset, const PairT &item, const int topk, size_t k, std::mutex &mut) {
    const auto [dist, id] = item;
    if(xset.find(id) != xset.end()) {
        DBG_ONLY(std::fprintf(stderr, "id %u is already present\n", int(id));)
        return;
    }
    if(topk <= 0 || x.size() < k) {
        std::lock_guard<std::mutex> lock(mut);
        xset.insert(id);
        x.push(item);
        return;
        //assert(*std::min_element(x.begin(), x.end()) == x.front());
    }
    //DBG_ONLY(std::fprintf(stderr, "Updating item %g/%u\n", item.first, item.second);)
    if(item.first <= x.top().first) {
        DBG_ONLY(std::fprintf(stderr, "New top before update: %g/Size %zu, with new item %g/%u added\n", x.front().first, x.size(), item.first, item.second);)
        std::lock_guard<std::mutex> lock(mut);
        auto old = x.top();
        // If the rank is the same as k - 1, save both
        // otherwise, discard the old one, since it's not good enough
        if(old.first != item.first) {
            xset.erase(old.second);
            x.pop();
        }
        x.push(item);
    } DBG_ONLY(else std::fprintf(stderr, "Count %g was not sufficient to be included. Current top: %g\n", item.first, x.top().first);)
}

#define ALL_CASE_NS\
               CASE_N(8, uint64_t);\
               CASE_N(4, uint32_t);\
               CASE_N(2, uint16_t);\
               CASE_N(1, uint8_t);\
                default: __builtin_unreachable();

std::vector<pqueue> build_index(SetSketchIndex<LSHIDType, LSHIDType> &idx, const Dashing2DistOptions &opts, const SketchingResult &result) {
    const bool index_compressed = opts.sketch_compressed_set;
    // Builds the LSH index and populates nearest-neighbor lists in parallel
    const size_t ns = result.names_.size();
    const int topk = opts.min_similarity_ > 0. ? -1: opts.num_neighbors_ > 0 ? 1: 0;
    static constexpr const LSHDistType INFLATE_FACTOR = 3.5;
    // Make the similarities negative so that the smallest items are the ones with the highest similarities
    size_t ntoquery = opts.num_neighbors_ <= 0 ? (maxcand_global <= 0 ? ns - 1: size_t(maxcand_global))
                                               : std::min(ns - 1, size_t(opts.num_neighbors_ * INFLATE_FACTOR));
    std::vector<pqueue> neighbor_lists(ns);
    if(opts.output_kind_ == KNN_GRAPH && opts.num_neighbors_ > 0)
        for(auto &n: neighbor_lists)
            n.reserve(opts.num_neighbors_);
    std::vector<flat_hash_set<LSHIDType>> neighbor_sets(ns);
    std::unique_ptr<std::mutex[]> mutexes(new std::mutex[ns]);
    auto idxstart = std::chrono::high_resolution_clock::now();
    // Build the index
    const bool indexing_compressed = index_compressed && opts.fd_level_ >= 1. && opts.fd_level_ < sizeof(RegT) && opts.kmer_result_ < FULL_MMER_SET;
    idx.size(ns);
    OMP_PFOR
    for(size_t i  = 0; i < ns; ++i) {
        if(indexing_compressed) {
            switch(int(opts.fd_level_)) {
#define CASE_N(digit, TYPE) case digit: idx.update(minispan<TYPE>((TYPE *)opts.compressed_ptr_ + opts.sketchsize_ * i, opts.sketchsize_), i); break
            ALL_CASE_NS
#undef CASE_N
            }
        } else {
            idx.update(minispan<RegT>(&result.signatures_[opts.sketchsize_ * i], opts.sketchsize_), i);
        }
    }
    auto idxstop = std::chrono::high_resolution_clock::now();
    // Build neighbor lists
    // Currently parallelizing the outer loop,
    // but the inner might be worth trying
    OMP_PFOR_DYN
    for(size_t id = 0; id < ns; ++id) {
        //std::fprintf(stderr, "%zu\t%zu\n", id, ns);
        std::tuple<std::vector<LSHIDType>, std::vector<uint32_t>, std::vector<uint32_t>> query_res;
        if(indexing_compressed && opts.fd_level_ >= 1. && opts.fd_level_ < sizeof(RegT) && opts.kmer_result_ < FULL_MMER_SET) {
            switch(int(opts.fd_level_)) {
#define CASE_N(i, TYPE) \
        case i: {query_res = idx.query_candidates(\
            minispan<TYPE>((TYPE *)opts.compressed_ptr_ + opts.sketchsize_ * id, opts.sketchsize_),\
            ntoquery);\
        } break
               ALL_CASE_NS
#undef CASE_N
            }
        } else
            query_res = idx.query_candidates(minispan<RegT>(&result.signatures_[opts.sketchsize_ * id], opts.sketchsize_), ntoquery);
        auto &[ids, counts, npr] = query_res;
        const size_t idn = ids.size();
        for(size_t j = 0; j < idn; ++j) {
            const LSHIDType oid = ids[j];
            if(id == oid) continue; // Don't track one's self
            const auto cd(-LSHDistType(counts[j]));
            update(neighbor_lists[oid], neighbor_sets[oid], PairT{cd, id}, topk, ntoquery, mutexes[oid]);
            update(neighbor_lists[id], neighbor_sets[id], PairT{cd, oid}, topk, ntoquery, mutexes[id]);
        }
    }
    OMP_PFOR_DYN
    for(size_t i = 0; i < ns; ++i)
        neighbor_lists[i].sort();
    auto knnstop = std::chrono::high_resolution_clock::now();

VERBOSE_ONLY(
    for(size_t id = 0; id < ns; ++id) {
        std::fprintf(stderr, "ID %s/%zu/%zu has a list of %zu long\n", result.names_[id].data(), id + 1, ns, neighbor_lists[id].size());
        for(size_t i = 0; i < neighbor_lists[id].size(); ++i) {
            size_t my_id = neighbor_lists[id][i].second;
            if(my_id > ns) {
                std::fprintf(stderr, "Impossible id (%d)? d = %g, id = %zu\n", int(i), neighbor_lists[id][i].first, my_id);
                std::exit(1);
            }
            std::fprintf(stderr, "First sequence list: %s/%zu\n", result.names_.at(neighbor_lists[id][i].second).data(), my_id);
        }
    }
)
    std::fprintf(stderr, "Index building took %Lgs. KNN Generation took %Lgs\n", std::chrono::duration<long double, std::ratio<1, 1>>(idxstop - idxstart).count(), std::chrono::duration<long double, std::ratio<1, 1>>(knnstop - idxstop).count());
    return neighbor_lists;
}
std::vector<pqueue> build_exact_graph(SetSketchIndex<LSHIDType, LSHIDType> &, const Dashing2DistOptions &opts, const SketchingResult &result, const bool) {
    // Builds the LSH index and populates nearest-neighbor lists in parallel
    const size_t ns = result.names_.size();
    //const int topk = opts.min_similarity_ > 0. ? -1: opts.num_neighbors_ > 0 ? 1: 0;
    // Make the similarities negative so that the smallest items are the ones with the highest similarities
    std::vector<pqueue> neighbor_lists(ns);
    if(opts.output_kind_ == KNN_GRAPH && opts.num_neighbors_ > 0)
        for(auto &n: neighbor_lists)
            n.reserve(opts.num_neighbors_);
    std::vector<flat_hash_set<LSHIDType>> neighbor_sets(ns);
    std::unique_ptr<std::mutex[]> mutexes(new std::mutex[ns]);
    auto idxstart = std::chrono::high_resolution_clock::now();
    auto idxstop = std::chrono::high_resolution_clock::now();
    // Build neighbor lists
    // Currently parallelizing the outer loop,
    // but the inner might be worth trying

    const bool isdist = distance(opts.measure_);
    const LSHDistType mult = isdist ? 1.: -1.;
    const double simt = opts.min_similarity_ > 0. ? opts.min_similarity_: 0.9;
    OMP_PFOR_DYN
    for(size_t id = 0; id < ns; ++id) {
        auto &nl = neighbor_lists[id];
        for(size_t rhid = 0; rhid < ns; ++rhid) {
            if(rhid == id) continue; // skip self.
            auto sim = mult * compare(opts, result, id, rhid);
            if(opts.output_kind_ == KNN_GRAPH) {
                // Don't include as a nearest-neighbor if the similarity is 0.
                if(!isdist && !sim)
                    continue;
                // If top-k is not filled, keep adding items.
                if(static_cast<std::ptrdiff_t>(nl.size()) < opts.num_neighbors_) {
                    nl.push({sim, rhid});
                } else {
                    const auto oldv = nl.top().first;
                    if(sim < oldv) {
                        nl.push({sim, rhid});
                        if(nl.size() > size_t(opts.num_neighbors_))
                            nl.pop();
                    } else {
                        // If the k-th best item is the same as this item, add it to the list so that we aren't ignoring equally-good-top-k items
                        if(sim == oldv) nl.push({sim, rhid});
                    }
                }
            } else {
                if(sim <= mult * simt)
                    nl.push({sim, rhid});
            }
        }
        nl.sort();
        if(nl.size() > size_t(opts.num_neighbors_)) {
            // In case we maintained a longer top-k than necessary (because items *were* in the top-k, but are no longer
            // we have to remove them.
            // This leaves any items sharing the same distance/similarity present, as x.first == kth_bestv for those items and therefore the lambda returns false.
            const auto kth_bestv = nl[opts.num_neighbors_ - 1].first;
            nl.erase(std::find_if(nl.begin() + opts.num_neighbors_, nl.end(), [kth_bestv](const auto &x) {return x.first > kth_bestv;}), nl.end());
        }
    }
    auto knnstop = std::chrono::high_resolution_clock::now();

    std::fprintf(stderr, "Index building took %Lgs. KNN Generation took %Lgs\n", std::chrono::duration<long double, std::ratio<1, 1>>(idxstop - idxstart).count(), std::chrono::duration<long double, std::ratio<1, 1>>(knnstop - idxstop).count());
    return neighbor_lists;
}

}
