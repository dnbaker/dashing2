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
    const auto [dist, id] = item;
    if(xset.find(id) != xset.end()) {
        std::fprintf(stderr, "id %u is already present\n", int(id));
        return;
    }
    if(topk <= 0 || x.size() < k) {
        std::fprintf(stderr, "Increasing size until it includes the k = %zu\n", k);
        std::lock_guard<std::mutex> lock(mut);
        xset.insert(id);
        x.emplace_back(item);
        std::push_heap(x.begin(), x.end());
        for(size_t i = 0; i < x.size(); ++i) {
            std::fprintf(stderr, "%g:", x[i].first);
        }
        std::fprintf(stderr, "New top: %g. Size: %zu\n", x.front().first, x.size());
        std::fprintf(stderr, "New bottom: %g. Size: %zu\n", x.front().first, x.size());

        assert(*std::max_element(x.begin(), x.end()) == x.front());
        //assert(*std::min_element(x.begin(), x.end()) == x.front());
    } else if(x.front() < item) {
        std::fprintf(stderr, "New top before update: %g/Size %zu\n", x.front().first, x.size());
        std::lock_guard<std::mutex> lock(mut);
        std::pop_heap(x.begin(), x.end());
        xset.erase(x.front().second);
        xset.insert(id);
        x.back() = item;
        std::push_heap(x.begin(), x.end());
        std::fprintf(stderr, "New top after update: %g/Size %zu\n", x.front().first, x.size());
    }
}

std::vector<std::vector<PairT>> build_index(SetSketchIndex<uint64_t, LSHIDType> &idx, Dashing2DistOptions &opts, const SketchingResult &result, const bool index_compressed) {
    // Builds the LSH index and populates nearest-neighbor lists in parallel
    const size_t ns = result.names_.size();
    const int topk = opts.min_similarity_ > 0. ? -1: opts.num_neighbors_ > 0 ? 1: 0;
    const LSHDistType INFLATE_FACTOR = 3.5;
    // Make the similarities negative so that the smallest items are the ones with the highest similarities
    size_t ntoquery = opts.num_neighbors_ <= 0 ? ns - 1: std::min(ns - 1, size_t(opts.num_neighbors_ * INFLATE_FACTOR));
    std::vector<std::vector<PairT>> neighbor_lists(ns);
    std::vector<ska::flat_hash_set<LSHIDType>> neighbor_sets(ns);
    std::unique_ptr<std::mutex[]> mutexes(new std::mutex[ns]);
    std::fprintf(stderr, "About to build index\n");
    for(size_t i  = 0; i < ns; ++i) {
        if(index_compressed && opts.fd_level_ >= 1. && opts.fd_level_ < sizeof(RegT) && opts.kmer_result_ < FULL_MMER_SET) {
            switch(int(opts.fd_level_)) {
#define CASE_N(i, TYPE) case i: {minispan<TYPE> mspan((TYPE *)opts.compressed_ptr_ + opts.sketchsize_ * i, opts.sketchsize_);\
    idx.update(mspan);} break
               CASE_N(8, uint64_t);
               CASE_N(4, uint32_t);
               CASE_N(2, uint16_t);
               CASE_N(1, uint8_t);
                default: __builtin_unreachable();
#undef CASE_N
            }
        } else {
            std::fprintf(stderr, "Indexing with span index %zu/%zu\n", i + 1, ns);
            minispan<RegT> mspan(&result.signatures_[opts.sketchsize_ * i], opts.sketchsize_);
            idx.update(mspan);
        }
        // res = update_query(,, x)
        // Build a candidate list and the index at the same time
        // When this is complete, we'll have a list of IDs to do comparisons with
    }
    for(size_t id = 0; id < ns; ++id) {
        std::fprintf(stderr, "%zu\t%zu\n", id, ns);
        //OMP_PFOR_DYN
        minispan<RegT> mspan(&result.signatures_[opts.sketchsize_ * id], opts.sketchsize_);
        auto [ids, counts, npr] = idx.query_candidates(mspan, ntoquery);
        const size_t idn = ids.size();
        for(size_t j = 0; j < idn; ++j) {
            const LSHIDType oid = ids[j];
            if(id == oid) continue; // Don't track one's self
            const PairT item{-LSHDistType(counts[j]), oid};
            update(neighbor_lists[ids[j]], neighbor_sets[ids[j]], item, topk, ntoquery, mutexes[ids[j]]);
            update(neighbor_lists[id], neighbor_sets[id], item, topk, ntoquery, mutexes[id]);
        }
    }
    
    for(size_t id = 0; id < ns; ++id) {
        std::fprintf(stderr, "ID %zu/%zu has a list of %zu long\n", id + 1, ns, neighbor_lists[id].size());
        for(size_t i = 0; i < neighbor_lists[id].size(); ++i) {
            size_t my_id = neighbor_lists[id][i].second;
            if(my_id > ns) {
                std::fprintf(stderr, "How do I have this id (%d)? d = %g, id = %zu\n", int(i), neighbor_lists[id][i].first, my_id);
                continue;
            }
            std::fprintf(stderr, "First sequence list: %s/%zu\n", result.names_.at(neighbor_lists[id][i].second).data(), my_id);
        }
    }
    return neighbor_lists;
}

}
