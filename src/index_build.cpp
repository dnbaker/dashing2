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
void pqueue::erase(typename std::priority_queue<PairT>::container_type::iterator it, typename std::priority_queue<PairT>::container_type::iterator oit) {
    this->c.erase(it, oit);
}


void update(pqueue &x, ska::flat_hash_set<LSHIDType> &xset, const PairT &item, const int topk, size_t k, std::mutex &mut) {
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
    DBG_ONLY(std::fprintf(stderr, "Updating item %g/%u\n", item.first, item.second);)
    if(x.top() < item) {
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
    }
}

#define ALL_CASE_NS\
               CASE_N(8, uint64_t);\
               CASE_N(4, uint32_t);\
               CASE_N(2, uint16_t);\
               CASE_N(1, uint8_t);\
                default: __builtin_unreachable();

std::vector<pqueue> build_index(SetSketchIndex<uint64_t, LSHIDType> &idx, Dashing2DistOptions &opts, const SketchingResult &result, const bool index_compressed) {
    // Builds the LSH index and populates nearest-neighbor lists in parallel
    const size_t ns = result.names_.size();
    const int topk = opts.min_similarity_ > 0. ? -1: opts.num_neighbors_ > 0 ? 1: 0;
    const LSHDistType INFLATE_FACTOR = 3.5;
    // Make the similarities negative so that the smallest items are the ones with the highest similarities
    size_t ntoquery = opts.num_neighbors_ <= 0 ? ns - 1: std::min(ns - 1, size_t(opts.num_neighbors_ * INFLATE_FACTOR));
    std::vector<pqueue> neighbor_lists(ns);
    if(opts.output_kind_ == KNN_GRAPH && opts.num_neighbors_ > 0)
        for(auto &n: neighbor_lists)
            n.reserve(opts.num_neighbors_);
    std::vector<ska::flat_hash_set<LSHIDType>> neighbor_sets(ns);
    std::unique_ptr<std::mutex[]> mutexes(new std::mutex[ns]);
    // Build the index
    for(size_t i  = 0; i < ns; ++i) {
        if(index_compressed && opts.fd_level_ >= 1. && opts.fd_level_ < sizeof(RegT) && opts.kmer_result_ < FULL_MMER_SET) {
            switch(int(opts.fd_level_)) {
#define CASE_N(digit, TYPE) case digit: idx.update(minispan<TYPE>((TYPE *)opts.compressed_ptr_ + opts.sketchsize_ * i, opts.sketchsize_)); break
            ALL_CASE_NS
#undef CASE_N
            }
        } else {
            idx.update(minispan<RegT>(&result.signatures_[opts.sketchsize_ * i], opts.sketchsize_));
        }
    }
    // Build neighbor lists
    // Currently parallelizing the outer loop,
    // but the inner might be worth trying
    OMP_PFOR_DYN
    for(size_t id = 0; id < ns; ++id) {
        //std::fprintf(stderr, "%zu\t%zu\n", id, ns);
        std::tuple<std::vector<LSHIDType>, std::vector<uint32_t>, std::vector<uint32_t>> query_res;
        if(index_compressed && opts.fd_level_ >= 1. && opts.fd_level_ < sizeof(RegT) && opts.kmer_result_ < FULL_MMER_SET) {
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
            const LSHDistType cd(counts[j]);
            update(neighbor_lists[oid], neighbor_sets[oid], PairT{cd, id}, topk, ntoquery, mutexes[oid]);
            update(neighbor_lists[id], neighbor_sets[id], PairT{cd, oid}, topk, ntoquery, mutexes[id]);
        }
    }

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
    return neighbor_lists;
}

}
