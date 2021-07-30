#include "d2.h"
#include "bwsketch.h"
#ifndef NOCURL
#define NOCURL 1
#endif
#include "bigWig.h"

namespace dashing2 {
std::vector<std::pair<int, bwOverlapIterator_t *>> get_iterators(bigWigFile_t *fp);
template<typename T>
void reduce_pair(T &lhs, T &rhs) {
    const size_t n = lhs.size();
    OMP_PRAGMA("omp simd")
    for(size_t i = 0; i < n; ++i){
        lhs[i] = std::min(lhs[i], rhs[i]);
    }
}
std::vector<RegT> reduce(const flat_hash_map<std::string, std::vector<RegT>> &map) {
    const size_t n = map.size();
    const unsigned int ln = static_cast<int>(std::ceil(n));
    std::vector<std::vector<RegT>> vals(n);
    {
        auto it = vals.data();
        for(auto &pair: map) *it++ = std::move(pair.second);
    }
    for(size_t i = 0; i < ln; ++i) {
        const size_t step_size = 1ull << i;
        const size_t sweep_size = (i + 1) << 1;
        const size_t nsweeps = (n + (sweep_size - 1)) / sweep_size;
        OMP_PFOR
        for(size_t j = 0; j < nsweeps; ++j) {
            const auto lh = j * sweep_size, rh = lh + step_size;
            if(rh < n)
                reduce_pair(vals[lh], vals[rh]);
        }
    }
    return std::move(vals.front());
}

std::vector<std::pair<int, bwOverlapIterator_t *>> get_iterators(bigWigFile_t *fp) {
    auto cp = fp->cl;
    const int nk = cp->nKeys;
    std::vector<std::pair<int, bwOverlapIterator_t *>> ret;
    for(auto i = 0; i < nk; ++i) {
        if(cp->len[i] < 1) continue;
        auto ptr = bwOverlappingIntervalsIterator(fp, fp->cl->chrom[i], 0, fp->cl->len[i], blocks_per_iter);
        if(!ptr->data || (ptr->intervals && ptr->intervals->l == 0)) bwIteratorDestroy(ptr);
        else ret.push_back({i, ptr});
    }
    auto it = std::find_if(ret.begin(), ret.end(), [](auto x) {return x.second == nullptr;});
    if(it != ret.end()) {
        std::fprintf(stderr, "%d/%p\n", it->first, (void *)it->second);
        std::exit(1);
    }
    const size_t npass = std::count_if(ret.begin(), ret.end(), [](auto x) {return x.second != nullptr;});
    if(ret.size() != npass) throw 1;
    //for(auto &p: ret) bwIteratorDestroy(p.second);
    return ret;
}

}
