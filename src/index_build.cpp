#include "cmp_main.h"
#include "sketch/ssi.h"
#include "minispan.h"
#include <vector>
#include <mutex>
#include "index_build.h"

namespace dashing2 {


// std::less -> vector.back() == heap.top(), where it is the smallest of the top
// if std::less<void>{}(vector.back(), newitem), then pop heap, add item and push heap
// We store similarities as -x in our nearest neighbor lists so that
// we only compile one heap update step -- for smallest distances

void update(std::vector<PairT> &x, ska::flat_hash_set<LSHIDType> &xset, const PairT &item, const int topk, size_t k, std::mutex &mut) {
    if(xset.find(item.second) != xset.end()) return;
    if(topk <= 0 || x.size() < k) {
        std::lock_guard<std::mutex> lock(mut);
        xset.insert(item.second);
        x.emplace_back(item);
        std::push_heap(x.begin(), x.end());
    } else if(x.back() < item) {
        std::lock_guard<std::mutex> lock(mut);
        std::pop_heap(x.begin(), x.end());
        xset.erase(x.back().second);
        xset.insert(item.second);
        x.back() = item;
        std::push_heap(x.begin(), x.end());
    }
}

std::vector<std::vector<PairT>> build_index(SetSketchIndex<uint64_t, LSHIDType> &idx, Dashing2DistOptions &opts, const SketchingResult &result, const bool index_compressed) {
    // Builds the LSH index and populates nearest-neighbor lists in parallel
    const size_t ns = result.names_.size();
    const int topk = opts.min_similarity_ > 0. ? -1: opts.num_neighbors_ > 0 ? 1: 0;
    const LSHDistType INFLATE_FACTOR = 3.5;
    // Make the similarities negative so that the smallest items are the ones with the highest similarities
    size_t ntoquery = opts.num_neighbors_ <= 0 ? ns: std::min(ns, size_t(opts.num_neighbors_ * INFLATE_FACTOR));
    std::vector<std::vector<PairT>> neighbor_lists(ns);
    std::vector<ska::flat_hash_set<LSHIDType>> neighbor_sets(ns);
    std::unique_ptr<std::mutex[]> mutexes(new std::mutex[ns]);
    for(size_t i  = 0; i < ns; ++i) {
        std::tuple<std::vector<LSHIDType>, std::vector<uint32_t>, std::vector<uint32_t>> idx_res;
        if(index_compressed && opts.fd_level_ >= 1. && opts.fd_level_ < sizeof(RegT) && opts.kmer_result_ < FULL_MMER_SET) {
            switch(int(opts.fd_level_)) {
#define CASE_N(i, TYPE) case i: {minispan<TYPE> mspan((TYPE *)opts.compressed_ptr_ + opts.sketchsize_ * i, opts.sketchsize_);\
    idx_res = idx.update_query(mspan, ntoquery);} break
               CASE_N(8, uint64_t);
               CASE_N(4, uint32_t);
               CASE_N(2, uint16_t);
               CASE_N(1, uint8_t);
                default: __builtin_unreachable();
#undef CASE_N
            }
        } else {
            minispan<RegT> mspan(&result.signatures_[opts.sketchsize_ * i], opts.sketchsize_);
            if(opts.kmer_result_ >= FULL_MMER_SET) idx_res = idx.update_query_bottomk(mspan, ntoquery);
            else idx_res = idx.update_query(mspan, ntoquery);
        }
        auto [ids, counts, perrows] = idx_res;
        const size_t idn = ids.size();
        OMP_PFOR_DYN
        for(size_t j = 0; j < idn; ++j) {
            const LSHIDType id = ids[j];
            const PairT item{id, -LSHDistType(counts[j])};
            update(neighbor_lists[j], neighbor_sets[j], item, topk, ntoquery, mutexes[j]);
            update(neighbor_lists[id], neighbor_sets[id], item, topk, ntoquery, mutexes[id]);
        }
        // Next, compute distances against those in the list
        // res = update_query(,, x)
        // Build a candidate list and the index at the same time
        // When this is complete, we'll have a list of IDs to do comparisons with
    }
    return neighbor_lists;
}

}
