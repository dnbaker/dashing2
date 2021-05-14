#include "d2.h"
#include "bwsketch.h"

namespace dashing2 {

static constexpr uint32_t default_BW_READ_BUFFER = 1<<30;

uint32_t BW_READ_BUFFER = default_BW_READ_BUFFER;
std::vector<RegT> reduce(const ska::flat_hash_map<std::string, std::vector<RegT>> &map);
std::vector<std::pair<int, bwOverlapIterator_t *>> get_iterators(bigWigFile_t *fp);

#ifndef BLOCKS_PER_ITER
#define BLOCKS_PER_ITER 4000000
#endif

using std::to_string;

BigWigSketchResult bw2sketch(std::string path, const ParseOptions &opts) {
    BigWigSketchResult ret;
    std::string cache_path = path.substr(0, path.find_last_of('.'));
    if(opts.cssize_) {
        cache_path = cache_path + "countsketch-" + to_string(opts.cssize_);
    }
    if(opts.sspace_ == SPACE_SET) cache_path += opts.one_perm_ ? ".opss": ".fss";
    else if(opts.sspace_ == SPACE_MULTISET) cache_path += ".bmh";
    else if(opts.sspace_ == SPACE_PSET) cache_path += ".pmh";
    else throw std::runtime_error("Unexpected space for BigWig sketching");
    if(opts.cache_sketches_ && !opts.by_chrom_ && bns::isfile(cache_path)) {
        auto nb = bns::filesize(cache_path.data());
        std::vector<RegT> save(nb / sizeof(RegT));
        std::FILE *ifp = std::fopen(cache_path.data(), "rb");
        if(!ifp) throw 1;
        if(std::fread(save.data(), nb, 1, ifp) != 1) {
            std::fprintf(stderr, "Failed to read from disk; instead, sketching from scratch (%s)\n", path.data());
        } else {
            ret.global_.reset(new std::vector<RegT>(std::move(save)));
        }
        std::fclose(ifp);
        return ret;
    }
    if(opts.count() != EXACT_COUNTING) {
        throw std::invalid_argument("Counting format must be exact for BigWigs. (No Count-Sketch approximation). This may change in the future.");
    }
    ska::flat_hash_map<std::string, std::vector<RegT>> retmap;
    std::fprintf(stderr, "Space: %s\n", opts.sspace_ == SPACE_SET ? "Set": opts.sspace_ == SPACE_MULTISET ? "Multist": opts.sspace_ == SPACE_PSET ? "Probdist": "Editdist");
    if(opts.sspace_ != SPACE_SET && opts.sspace_ != SPACE_MULTISET && opts.sspace_ != SPACE_PSET)
        throw std::invalid_argument("Can't do edit distance for BigWig files");
    if(bwInit(BW_READ_BUFFER)) {
        std::fprintf(stderr, "Error in initializing bigwig\n");
        return ret;
    }
    bigWigFile_t *fp = bwOpen(path.data(), nullptr, "r");
    if(fp == nullptr) throw std::runtime_error("Could not open bigwigfile");

    auto ids(get_iterators(fp));
    for(const auto &p: ids) retmap.emplace(fp->cl->chrom[p.first], std::vector<RegT>());
    const size_t ss = opts.sketchsize_;
    std::vector<FullSetSketch> fss;
    std::vector<OPSetSketch> opss;
    std::vector<BagMinHash> bmhs;
    std::vector<ProbMinHash> pmhs;

    for(size_t i = 0; i < std::min(size_t(opts.nthreads()), ids.size()); ++i) {
        if(opts.sspace_ == SPACE_SET) {
            if(opts.one_perm_) opss.emplace_back(ss);
            else fss.emplace_back(ss);
        } else if(opts.sspace_ == SPACE_MULTISET) {
            bmhs.emplace_back(ss);
        } else {
            pmhs.emplace_back(ss);
        }
    }
    OMP_PFOR_DYN
    for(size_t i = 0; i < ids.size(); ++i) {
#ifndef NDEBUG
        std::fprintf(stderr, "Processing contig %zu/%zu\n", i, ids.size());
#endif
        const int tid = OMP_ELSE(omp_get_thread_num(), 0);
        const int contig_id = ids[i].first;
        std::string chrom = fp->cl->chrom[contig_id];
        const uint64_t chrom_hash = std::hash<std::string>{}(chrom);
        auto &rvec = retmap[chrom];
        auto &ptr = ids[i].second;
        if(ptr == nullptr) {
            std::fprintf(stderr, "bwOverlapIterator_t * is null when it should be non-null (tid = %d/contig = %d)\n", tid, contig_id);
            std::exit(1);
        }
        {
            if(fss.size()) {fss[tid].reset();}
            else if(opss.size()) {opss[tid].reset();}
            else if(bmhs.size()) {bmhs[tid].reset();}
            else if(pmhs.size()) {pmhs[tid].reset();}
            else throw std::invalid_argument("Not supported: sketching besides PMH, BMH, SetSketch for BED files");
            for(;ptr->data;ptr = bwIteratorNext(ptr)) {
                const uint32_t numi = ptr->intervals->l;
                float *vptr = ptr->intervals->value;
                for(uint32_t j = 0; j < numi; ++j) {
                    auto istart = ptr->intervals->start[j];
                    auto iend = ptr->intervals->end[j];
                    //std::fprintf(stderr, "%s:%u->%u [%u/%u]\n", chrom.data(), istart, iend, j, numi);
                    for(;istart < iend;++istart) {
                        auto k = chrom_hash ^ istart;
                        auto v = vptr[j];
                        if(fss.size()) fss[tid].update(k);
                        else if(opss.size()) opss[tid].update(k);
                        else if(bmhs.size()) bmhs[tid].update(k, v);
                        else if(pmhs.size()) pmhs[tid].update(k, v);
                    }
                }
#ifndef NDEBUG
                std::fprintf(stderr, "processed %u intervals. Now loading next batch (%s)\n", numi, chrom.data());
#endif
                ptr = bwIteratorNext(ptr);
            } while(ptr->data);

            auto newvec = fss.size() ? fss[tid].to_sigs() :
                          opss.size() ? opss[tid].to_sigs() :
                          bmhs.size() ? bmhs[tid].to_sigs(): pmhs[tid].to_sigs();

            if(rvec.empty()) {
                rvec = newvec;
            } else {
#ifdef _OPENMP
                #pragma omp simd
#endif
                for(size_t i = 0; i < rvec.size(); ++i) {
                    rvec[i] = std::min(rvec[i], newvec[i]);
                }
            }
        }
        bwIteratorDestroy(ptr);
    }

    bwClose(fp);
    bwCleanup();
    ret.global_.reset(new std::vector<RegT>(std::move(reduce(retmap))));
    if(opts.by_chrom_)
        ret.chrmap_.reset(new ska::flat_hash_map<std::string, std::vector<RegT>>(std::move(retmap)));
    return ret;
}
std::vector<std::pair<int, bwOverlapIterator_t *>> get_iterators(bigWigFile_t *fp) {
    auto cp = fp->cl;
    const int nk = cp->nKeys;
    std::vector<std::pair<int, bwOverlapIterator_t *>> ret;
    for(auto i = 0; i < nk; ++i) {
        if(cp->len[i] < 1) continue;
        auto ptr = bwOverlappingIntervalsIterator(fp, fp->cl->chrom[i], 0, fp->cl->len[i], BLOCKS_PER_ITER);
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

} // dashing2
