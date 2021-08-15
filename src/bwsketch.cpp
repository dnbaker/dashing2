#include "d2.h"
#include "bwsketch.h"
#ifndef NOCURL
#define NOCURL 1
#endif
#include "bigWig.h"

namespace dashing2 {

static constexpr uint32_t default_BW_READ_BUFFER = 1<<30;

uint32_t BW_READ_BUFFER = default_BW_READ_BUFFER;
std::vector<RegT> reduce(const flat_hash_map<std::string, std::vector<RegT>> &map);
std::vector<std::pair<int, bwOverlapIterator_t *>> get_iterators(bigWigFile_t *fp);


using std::to_string;



BigWigSketchResult bw2sketch(std::string path, const Dashing2Options &opts, bool parallel_process) {
    BigWigSketchResult ret;
    std::string cache_path = path.substr(0, path.find_last_of('.'));
    cache_path += to_suffix(opts);
    if(opts.trim_folder_paths()) {
        DBG_ONLY(std::fprintf(stderr, "Cached path before trimming: %s\n", cache_path.data());)
        cache_path = trim_folder(path);
        DBG_ONLY(std::fprintf(stderr, "Cached path after trimming: %s\n", cache_path.data());)
        if(opts.outprefix_.size()) {
            cache_path = opts.outprefix_ + '/' + cache_path;
        }
    }
    if(opts.seedseed_ != 0)
        cache_path += ".seed" + std::to_string(opts.seedseed_);
    if(opts.kmer_result_ <= FULL_SETSKETCH)
        cache_path = cache_path + std::string(".sketchsize") + std::to_string(opts.sketchsize_);
    cache_path = cache_path + std::string(".k") + std::to_string(opts.k_);
    if(opts.count_threshold_ > 0) {
        cache_path = cache_path + ".ct_threshold";
        if(std::fmod(opts.count_threshold_, 1.)) cache_path = cache_path + std::to_string(opts.count_threshold_);
        else cache_path = cache_path + std::to_string(int(opts.count_threshold_));
    }
    cache_path += ".";
    cache_path += opts.kmer_result_ <= FULL_SETSKETCH ? to_string(opts.sspace_): to_string(opts.kmer_result_);
    DBG_ONLY(std::fprintf(stderr, "Cache path: %s. isfile: %d\n", cache_path.data(), bns::isfile(cache_path));)
    if(opts.cache_sketches_ && !opts.by_chrom_ && bns::isfile(cache_path)) {
        std::FILE *ifp = xopen(cache_path);
        std::fread(&ret.card_, sizeof(ret.card_), 1, ifp);
        auto res = new std::vector<RegT>;
        while(!std::feof(ifp)) {
            RegT v;
            std::fread(&v, sizeof(v), 1, ifp);
            res->push_back(v);
        }
        ::pclose(ifp);
        ret.global_.reset(res);
        return ret;
    }
    if(opts.count() != EXACT_COUNTING) {
        THROW_EXCEPTION(std::invalid_argument("Counting format must be exact for BigWigs. (No Count-Sketch approximation). This may change in the future."));
    }
    flat_hash_map<std::string, std::vector<RegT>> retmap;
    DBG_ONLY(std::fprintf(stderr, "Space: %s\n", to_string(opts.sspace_).data());)
    if(opts.sspace_ != SPACE_SET && opts.sspace_ != SPACE_MULTISET && opts.sspace_ != SPACE_PSET)
        THROW_EXCEPTION(std::invalid_argument("Can't do edit distance for BigWig files"));
    if(bwInit(BW_READ_BUFFER)) {
        DBG_ONLY(std::fprintf(stderr, "Error in initializing bigwig\n");)
        return ret;
    }
    bigWigFile_t *fp = bwOpen(path.data(), nullptr, "r");
    if(fp == nullptr) THROW_EXCEPTION(std::runtime_error("Could not open bigwigfile"));

    if(parallel_process || opts.by_chrom_) {
        auto ids(get_iterators(fp));
        for(const auto &p: ids) retmap.emplace(fp->cl->chrom[p.first], std::vector<RegT>());
        const size_t ss = opts.sketchsize_;

        std::vector<FullSetSketch> fss;
        std::vector<OPSetSketch> opss;
        std::vector<BagMinHash> bmhs;
        std::vector<ProbMinHash> pmhs;
        for(size_t i = 0; i < std::min(size_t(opts.nthreads()), ids.size()); ++i) {
            if(opts.sspace_ == SPACE_SET) {
                if(opts.one_perm()) {
                    opss.emplace_back(ss);
                    if(opts.count_threshold_ > 1) opss.back().set_mincount(opts.count_threshold_);
                } else
                    fss.emplace_back(opts.count_threshold_, ss);
            } else if(opts.sspace_ == SPACE_MULTISET)
                bmhs.emplace_back(ss);
            else
                pmhs.emplace_back(ss);
        }
        long double total_weight = 0.;
        OMP_PRAGMA("omp parallel for schedule(dynamic) reduction(+:total_weight)")
        for(size_t i = 0; i < ids.size(); ++i) {
            DBG_ONLY(std::fprintf(stderr, "Processing contig %zu/%zu\n", i, ids.size());)
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
                else THROW_EXCEPTION(std::invalid_argument("Not supported: sketching besides PMH, BMH, SetSketch for BigWig files"));
                for(;ptr->data;ptr = bwIteratorNext(ptr)) {
                    const uint32_t numi = ptr->intervals->l;
                    float *vptr = ptr->intervals->value;
#define DO_FOR_SKETCH(item) do {\
                    for(uint32_t j = 0; j < numi; ++j) { \
                        for(auto istart = ptr->intervals->start[j], iend = ptr->intervals->end[j];istart < iend;(item).update(chrom_hash ^ istart++, vptr[j]));\
                    } } while(0)
#define DO_FOR_UNWEIGHTED_SKETCH(item) do {\
                    for(uint32_t j = 0; j < numi; ++j) { \
                        for(auto istart = ptr->intervals->start[j], iend = ptr->intervals->end[j];istart < iend;(item).update(chrom_hash ^ istart++));\
                    } } while(0)
                    if(fss.size()) {
                        DO_FOR_UNWEIGHTED_SKETCH(fss[tid]);
                    } else if(opss.size()) {
                        DO_FOR_UNWEIGHTED_SKETCH(opss[tid]);
                    } else if(bmhs.size()) {
                        DO_FOR_SKETCH(bmhs[tid]);
                    } else if(pmhs.size()) {
                        DO_FOR_SKETCH(pmhs[tid]);
                    }
#ifndef NDEBUG
                    std::fprintf(stderr, "processed %u intervals. Now loading next batch (%s)\n", numi, chrom.data());
#endif
                }

                auto newvec = fss.size() ? fss[tid].to_sigs() :
                              opss.size() ? opss[tid].to_sigs() :
                              bmhs.size() ? bmhs[tid].to_sigs(): pmhs[tid].to_sigs();
                total_weight += (fss.size() ? fss[tid].getcard(): opss.size() ? opss[tid].getcard(): bmhs.size() ? bmhs[tid].total_weight(): pmhs.size() ? pmhs[tid].total_weight(): RegT(-1));

                if(rvec.empty()) {
                    rvec = newvec;
                } else {
                    OMP_PRAGMA("omp simd")
                    for(size_t i = 0; i < rvec.size(); ++i) {
                        rvec[i] = std::min(rvec[i], newvec[i]);
                    }
                }
            }
            bwIteratorDestroy(ptr);
        }
        ret.card_ = total_weight;
        ret.global_.reset(new std::vector<RegT>(std::move(reduce(retmap))));
    } else {
        DBG_ONLY(auto timestart = std::chrono::high_resolution_clock::now();)
        auto cp = fp->cl;
        const int nk = cp->nKeys;
        std::unique_ptr<FullSetSketch> fss;
        std::unique_ptr<OPSetSketch> opss;
        std::unique_ptr<ProbMinHash> pmh;
        std::unique_ptr<BagMinHash> bmh;
        if(opts.sspace_ == SPACE_SET) {
            if(opts.one_perm()) {
                opss.reset(new OPSetSketch(opts.sketchsize_));
                if(opts.count_threshold_ > 1) opss->set_mincount(opts.count_threshold_);
            } else
                fss.reset(new FullSetSketch(opts.count_threshold_, opts.sketchsize_));
        } else if(opts.sspace_ == SPACE_MULTISET)
            bmh.reset(new BagMinHash(opts.sketchsize_));
        else
            pmh.reset(new ProbMinHash(opts.sketchsize_));
        for(auto i = 0; i < nk; ++i) {
            if(cp->len[i] < 1) continue;
            auto cid = fp->cl->chrom[i];
            auto ptr = bwOverlappingIntervalsIterator(fp, cid, 0, fp->cl->len[i], blocks_per_iter);
            if(!ptr->data || (ptr->intervals && ptr->intervals->l == 0)) {
                bwIteratorDestroy(ptr);
                continue;
            }
            const uint64_t chrom_hash = std::hash<std::string>{}(cid);
            const uint32_t numi = ptr->intervals->l;
            float *vptr = ptr->intervals->value;
            if(fss) {
                DO_FOR_UNWEIGHTED_SKETCH(*fss);
            } else if(opss) {
                DO_FOR_UNWEIGHTED_SKETCH(*opss);
            } else if(bmh) {
                DO_FOR_SKETCH(*bmh);
            } else if(pmh) {
                DO_FOR_SKETCH(*pmh);
            }
        }
        ret.card_ = (fss ? fss->getcard(): opss ? opss->getcard(): bmh ? bmh->total_weight(): pmh ? pmh->total_weight(): std::numeric_limits<RegT>::infinity());
        if(ret.card_ == std::numeric_limits<RegT>::infinity()) {
            std::fprintf(stderr, "Warning: infinite cardinality\n");
        }
        ret.global_.reset(new std::vector<RegT>(bmh ? bmh->to_sigs(): pmh ? pmh->to_sigs(): opss ? opss->to_sigs(): fss ? fss->to_sigs(): std::vector<RegT>()));
        DBG_ONLY(auto timestop = std::chrono::high_resolution_clock::now();
                 std::fprintf(stderr, "Took %gms to sketch\n", std::chrono::duration<double>(timestop - timestart).count());)
    }

    bwClose(fp);
    bwCleanup();
    if(opts.by_chrom_)
        ret.chrmap_.reset(new flat_hash_map<std::string, std::vector<RegT>>(std::move(retmap)));
    if(opts.kmer_result_ <= FULL_SETSKETCH) {
        std::FILE *ofp = std::fopen(cache_path.data(), "wb");
        if(!ofp) THROW_EXCEPTION(std::runtime_error(std::string("Could not open file at ") + cache_path + " for writing"));
        std::fwrite(&ret.card_, sizeof(ret.card_), 1, ofp);
        checked_fwrite(ofp, ret.global_->data(), sizeof(RegT) * ret.global_->size());
        std::fclose(ofp);
    } else {
        std::fprintf(stderr, "Warning: only SetSketch, OnePermSetSketch, ProbMinHash, and BagMinHash are cached to disk for BigWigs. Nothing being cached to disk.\n");
    }
    return ret;
}
#undef DO_FOR_SKETCH

} // dashing2
