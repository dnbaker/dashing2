#include "d2.h"


static constexpr uint32_t default_BW_READ_BUFFER = 1<<30;

uint32_t BW_READ_BUFFER = default_BW_READ_BUFFER;

#ifndef BLOCKS_PER_ITER
#define BLOCKS_PER_ITER 4000000
#endif

std::vector<std::pair<int, bwOverlapIterator_t *>> get_iterators(bigWigFile_t *fp) {
    auto cp = fp->cl;
    const int nk = cp->nKeys;
    std::vector<std::pair<int, bwOverlapIterator_t *>> ret(nk, {0, static_cast<bwOverlapIterator_t *>(nullptr)});
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(auto i = 0; i < nk; ++i) {
        if(cp->len[i] < 1) continue;
        auto &ref = ret[i].second;
        ref = bwOverlappingIntervalsIterator(fp, fp->cl->chrom[i], 0, fp->cl->len[i], BLOCKS_PER_ITER);
        ret[i].first = i;
        if(!ref->data) bwIteratorDestroy(ref), ref = nullptr;
    }
    ret.erase(std::remove_if(ret.begin(), ret.end(), [](auto x) {return x.second == nullptr;}));
    return ret;
}
ska::flat_hash_map<std::string, std::vector<RegT>> bw2sketch(std::string path, const ParseOptions &opts) {
    if(opts.count_ != EXACT_COUNTING) {
        throw std::invalid_argument("Counting format must be exact for BigWigs. (No Count-Sketch approximation). This may change in the future.");
    }
    ska::flat_hash_map<std::string, std::vector<RegT>> ret;
    if(opts.sspace_ != SPACE_SET && opts.sspace_ != SPACE_MULTISET && opts.sspace_ != SPACE_PSET)
        throw std::invalid_argument("Can't do edit distance for BigWig files");
    if(bwInit(BW_READ_BUFFER)) {
        std::fprintf(stderr, "Error in initializing bigwig\n");
        return ret;
    }
    bigWigFile_t *fp = bwOpen(path.data(), nullptr, "r");
    if(fp == nullptr) throw std::runtime_error("Could not open bigwigfile");

    const size_t nthreads = opts.nthreads();
    std::vector<std::pair<int, bwOverlapIterator_t *>> itpairs = get_iterators(fp);
    for(const auto &p: itpairs) ret.emplace(fp->cl->chrom[p.first], std::vector<RegT>());
    const size_t nbuffers = std::min(nthreads, itpairs.size());
    const size_t ss = opts.sketchsize_;
    std::vector<FullSetSketch> fss;
    std::vector<OPSetSketch> opss;
    std::vector<BagMinHash> bmhs;
    std::vector<ProbMinHash> pmhs;
    if(opts.sspace_ == SPACE_SET) {
        if(opts.one_perm_) std::generate_n(std::back_inserter(opss), nbuffers, [ss]{return OPSetSketch(ss);});
        else std::generate_n(std::back_inserter(fss), nbuffers, [ss]{return FullSetSketch(ss);});
    }
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(size_t i = 0; i < itpairs.size(); ++i) {
        const int tid = OMP_ELSE(omp_get_thread_num(), 0);
        const int contig_id = itpairs[i].first;
        std::string chrom = fp->cl->chrom[contig_id];
        const uint64_t chrom_hash = std::hash<std::string>{}(chrom);
        auto &rvec = ret[chrom];
        if(unlikely(itpairs[i].second == nullptr)) {
            std::fprintf(stderr, "bwOverlapIterator_t * is null when it should be non-null (tid = %d/contig = %d)\n", tid, contig_id);
            std::exit(1);
        }
        if(opts.sspace_ == SPACE_SET) {
            std::fprintf(stderr, "SPACE_SET is set -- sketching only positions, not counts. If this is unintentional, see usage.\n");
            if(opts.one_perm_) {
                auto &sketch = opss[tid];
                sketch.reset();
                do {
                    const uint32_t numi = itpairs[i].second->intervals->l;
                    for(uint32_t j = 0; j < numi; ++j) {
                        auto istart = itpairs[i].second->intervals->start[j];
                        auto iend = itpairs[i].second->intervals->end[j];
                        for(;istart < iend;sketch.update(chrom_hash ^ istart++));
                    }
                    itpairs[i].second = bwIteratorNext(itpairs[i].second);
                } while(itpairs[i].second->data);
                if(rvec.empty()) rvec = sketch.to_sigs();
                else {
                    // Pair-wise minima
                    std::transform(rvec.begin(), rvec.end(), sketch.data(), rvec.begin(),
                                   [](auto x, auto y) {return std::min(x, y);});
                }
            } else {
                auto &sketch = fss[tid];
                sketch.reset();
                do {
                    const uint32_t numi = itpairs[i].second->intervals->l;
                    for(uint32_t j = 0; j < numi; ++j) {
                        auto istart = itpairs[i].second->intervals->start[j];
                        auto iend = itpairs[i].second->intervals->end[j];
                        for(;istart < iend;sketch.update(chrom_hash ^ istart++));
                    }
                    itpairs[i].second = bwIteratorNext(itpairs[i].second);
                } while(itpairs[i].second->data);
                if(rvec.empty()) rvec = sketch.to_sigs();
                else {
                    // Pair-wise minima
                    std::transform(rvec.begin(), rvec.end(), sketch.data(), rvec.begin(),
                                   [](auto x, auto y) {return std::min(x, y);});
                }
            }
        } else {
            if(bmhs.size()) {
                auto &sketch = bmhs[tid];
                sketch.reset();
                do {
                    const uint32_t numi = itpairs[i].second->intervals->l;
                    float *vptr = itpairs[i].second->intervals->value;
                    for(uint32_t j = 0; j < numi; ++j) {
                        auto istart = itpairs[i].second->intervals->start[j];
                        auto iend = itpairs[i].second->intervals->end[j];
                        for(;istart < iend;sketch.update(chrom_hash ^ istart, vptr[j]));
                    }
                    itpairs[i].second = bwIteratorNext(itpairs[i].second);
                } while(itpairs[i].second->data);
                if(rvec.empty()) rvec = sketch.to_sigs();
                else {
                    // Pair-wise minima
                    std::transform(rvec.begin(), rvec.end(), sketch.data(), rvec.begin(),
                                   [](auto x, auto y) {return std::min(x, y);});
                }
            } else if(pmhs.size()) {
                auto &sketch = pmhs[tid];
                sketch.reset();
                do {
                    const uint32_t numi = itpairs[i].second->intervals->l;
                    float *vptr = itpairs[i].second->intervals->value;
                    for(uint32_t j = 0; j < numi; ++j) {
                        auto istart = itpairs[i].second->intervals->start[j];
                        auto iend = itpairs[i].second->intervals->end[j];
                        for(;istart < iend;sketch.update(chrom_hash ^ istart, vptr[j]));
                    }
                    itpairs[i].second = bwIteratorNext(itpairs[i].second);
                } while(itpairs[i].second->data);
                if(rvec.empty()) rvec = sketch.to_sigs();
                else {
                    // Pair-wise minima
                    std::transform(rvec.begin(), rvec.end(), sketch.data(), rvec.begin(),
                                   [](auto x, auto y) {return std::min(x, y);});
                }
            } else throw std::invalid_argument("Not supported: sketching besides PMH, BMH, SetSketch for BED files");
        }
    }
    for(auto &pair: itpairs) bwIteratorDestroy(pair.second);

    bwClose(fp);
#ifndef NOCURL
    bwCleanup();
#endif
    return ret;
}
std::vector<RegT> reduce(ska::flat_hash_map<std::string, std::vector<RegT>> &map) {
    std::vector<std::vector<RegT> *> ptrs(map.size());
    auto it = ptrs.begin();
    for(auto &pair: map) *it = &pair.second, ++it;
    const size_t n = map.size();
    const unsigned int ln = static_cast<int>(std::ceil(n));
    auto reduce = [&](auto &lhs, auto &rhs) {
        std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), [](auto x, auto y) {return std::min(x, y);});
    };
    for(size_t i = 0; i < ln; ++i) {
        const size_t step_size = 1 << i;
        const size_t sweep_size = (i + 1) << 1;
        const size_t nsweeps = (n + (sweep_size - 1)) / sweep_size;
        OMP_PFOR
        for(size_t j = 0; j < nsweeps; ++j) {
            const auto lh = j * sweep_size, rh = lh + step_size;
            if(rh < n)
                reduce(*ptrs[lh], *ptrs[rh]);
        }
    }
    return std::move(*ptrs[0]);
}
std::vector<RegT> reduce(const ska::flat_hash_map<std::string, std::vector<RegT>> &map) {
    std::vector<std::vector<RegT>> vals(map.size());
    auto it = vals.begin();
    for(auto &pair: map) *it = pair.second, ++it;
    const size_t n = map.size();
    const unsigned int ln = static_cast<int>(std::ceil(n));
    auto reduce = [&](auto &lhs, auto &rhs) {
        std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), [](auto x, auto y) {return std::min(x, y);});
    };
    for(size_t i = 0; i < ln; ++i) {
        const size_t step_size = 1 << i;
        const size_t sweep_size = (i + 1) << 1;
        const size_t nsweeps = (n + (sweep_size - 1)) / sweep_size;
        OMP_PFOR
        for(size_t j = 0; j < nsweeps; ++j) {
            const auto lh = j * sweep_size, rh = lh + step_size;
            if(rh < n)
                reduce(vals[lh], vals[rh]);
        }
    }
    return std::move(vals.front());
}
