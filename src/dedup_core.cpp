#include "cmp_main.h"
#include "src/ssi.h"
#include "minispan.h"
#include "fmt/format.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "dedup_core.h"


namespace dashing2 {

static constexpr size_t MINCAND = 12;
int exhaustive_dedup = 0;

#ifndef EARLYSTOP
#define EARLYSTOP 1
#endif

struct GreedyClustering {
    std::vector<LSHIDType> ids_;
    std::vector<std::vector<LSHIDType>> constituents_;
    sketch::lsh::SetSketchIndex<LSHIDType, LSHIDType> idx_;
    double simt, mult;
    bool earlystop;
    size_t mincand;
    const SketchingResult &result;
    const Dashing2DistOptions &opts;
    GreedyClustering(const SketchingResult &rs, const Dashing2DistOptions &opts, const sketch::lsh::SetSketchIndex<LSHIDType, LSHIDType> &idx, bool earlystop=EARLYSTOP, size_t mincand=MINCAND)
        : idx_(idx.clone()),
        earlystop(earlystop), mincand(mincand), result(rs), opts(opts)
    {
        simt = opts.min_similarity_ > 0. ? opts.min_similarity_: 0.9; // 90% is the default cut-off for deduplication
        mult = distance(opts.measure_)  ? 1.: -1.;
    }
    GreedyClustering &operator +=(const GreedyClustering &o) {
        const size_t osz = o.ids_.size();
        for(size_t i = 0; i < osz; ++i) {
            auto &orep = o.ids_[i];
            auto &ocon = o.constituents_[i];
            const minispan<RegT> span(&result.signatures_[opts.sketchsize_ * orep], opts.sketchsize_);
            auto [hits, counts, nper] = idx_.query_candidates(span, mincand, size_t(-1), earlystop);
            std::vector<LSHDistType> vals(hits.size());
            auto vp = vals.data();
            const auto vps = vp;
            for(const auto id: hits) {
                *vp++ = mult * compare(opts, result, orep, ids_[id]);
            }
            auto mv = std::min_element(vps, vp);
            if(hits.empty() || (mv != vp && mult * *mv < simt)) {
                ids_.push_back(orep);
                constituents_.emplace_back(std::move(ocon));
                continue;
            }
            auto pos = mv - vps;
            auto cluster_id = hits[pos];
            auto &cv = constituents_[cluster_id];
            auto &rep = ids_[cluster_id];
            if(result.cardinalities_[orep] > result.cardinalities_[rep]) {
                cv.push_back(rep);
                rep = orep;
            } else cv.push_back(orep);
            if(ocon.size())
                cv.insert(cv.end(), ocon.begin(), ocon.end());
        }
        return *this;
    }
};

template<typename T>
void par_reduce(T *x, size_t n) {
    const unsigned int ln = static_cast<int>(std::ceil(std::log2(n)));
    for(size_t i = 0; i < ln; ++i) {
        const size_t step_size = size_t(1) << i;
        const size_t sweep_size = (i + 1) << 1;
        const size_t nsweeps = (n + (sweep_size - 1)) / sweep_size;
        OMP_PFOR
        for(size_t j = 0; j < nsweeps; ++j) {
            const auto lh = j * sweep_size, rh = lh + step_size;
            if(rh < n)
                x[lh] += x[rh];
        }
    }
}



void update_res_mt(LSHIDType oid, std::vector<LSHIDType> &ids, std::vector<std::vector<LSHIDType>> &constituents,
                   sketch::lsh::SetSketchIndex<LSHIDType, LSHIDType> &idx, const Dashing2DistOptions &opts, const SketchingResult &result,
                   bool earlystop=EARLYSTOP, size_t mincand=MINCAND)
{
    const double simt = opts.min_similarity_ > 0. ? opts.min_similarity_: 0.9; // 90% is the default cut-off for deduplication
    assert(opts.sketchsize_ * oid < result.signatures_.size());
    const minispan<RegT> span(&result.signatures_[opts.sketchsize_ * oid], opts.sketchsize_);
    auto [hits, counts, nper] = idx.query_candidates(span, mincand, size_t(-1), earlystop);
    const size_t nh = hits.size();
    std::vector<LSHDistType> vals(hits.size());
    const LSHDistType mult = distance(opts.measure_) ? 1.: -1.;
    OMP_PFOR_DYN
    for(size_t i = 0; i < nh; ++i) {
        const auto id = hits[i];
        assert(id < ids.size());
        vals[i] = mult * compare(opts, result, oid, ids[id]);
    }
    auto mv = std::min_element(vals.begin(), vals.end());
    if(hits.empty() || (mv != vals.end() && mult * *mv < simt)) {
        //DBG_ONLY(if(mv != vals.end()) std::fprintf(stderr, "mult* mv: %g. simt: %g\n", mult * *mv, simt);)
        ids.push_back(oid);
        constituents.emplace_back();
        idx.update_mt(const minispan<RegT>(&result.signatures_[opts.sketchsize_ * oid], opts.sketchsize_));
    } else {
        auto pos = mv - vals.begin();
        assert(size_t(pos) < hits.size());
        auto cluster_id = hits[pos];
        assert(ids.size() == constituents.size());
        assert(cluster_id < constituents.size());
        assert(cluster_id < ids.size());
        auto &cv = constituents[cluster_id];
        auto &rep = ids[cluster_id];
        cv.push_back(oid);
        //DBG_ONLY(std::fprintf(stderr, "Cluster with rep %s is adding new item named %s\n", result.names_[ids[cluster_id]].data(), result.names_[oid].data());)
        if(result.cardinalities_[cv.back()] > result.cardinalities_[rep]) {
            // In case the items are unsorted with respect to cardinality due to the parallelism, swap it out.
            // That way, we'll keep the highest-cardinality set as the representative
            std::swap(cv.back(), rep);
        }
    }
}

void update_res(LSHIDType oid, std::vector<LSHIDType> &ids, std::vector<std::vector<LSHIDType>> &constituents,
                sketch::lsh::SetSketchIndex<LSHIDType, LSHIDType> &idx, const Dashing2DistOptions &opts, const SketchingResult &result,
                bool earlystop=EARLYSTOP, size_t mincand=MINCAND)
{
    const double simt = opts.min_similarity_ > 0. ? opts.min_similarity_: 0.9; // 90% is the default cut-off for deduplication
    assert(opts.sketchsize_ * oid < result.signatures_.size());
    const minispan<RegT> span(&result.signatures_[opts.sketchsize_ * oid], opts.sketchsize_);
    auto [hits, counts, nper] = idx.query_candidates(span, mincand, size_t(-1), earlystop);
    std::vector<LSHDistType> vals(hits.size());
    auto vp = vals.begin();
    const LSHDistType mult = distance(opts.measure_) ? 1.: -1.;
    for(const auto id: hits) {
        assert(id < ids.size());
        *vp++ = mult * compare(opts, result, oid, ids[id]);
    }
    auto mv = std::min_element(vals.begin(), vals.end());
    if(hits.empty() || (mv != vals.end() && mult * *mv < simt)) {
        ids.push_back(oid);
        constituents.emplace_back();
        const minispan<RegT> mp(&result.signatures_[opts.sketchsize_ * oid], opts.sketchsize_);
        idx.update(mp);
        //DBG_ONLY(std::fprintf(stderr, "Added item %zu/%s; %zu  hits, %zu clusters so far (%%%0.4g)\n", myid, result.names_[oid].data(), hits.size(), ids.size(), ids.size() * 100. / result.names_.size());)
    } else {
        auto pos = mv - vals.begin();
        assert(size_t(pos) < hits.size());
        auto cluster_id = hits[pos];
        assert(ids.size() == constituents.size());
        assert(cluster_id < constituents.size());
        assert(cluster_id < ids.size());
        auto &cv = constituents[cluster_id];
        auto &rep = ids[cluster_id];
        cv.push_back(oid);
        //DBG_ONLY(std::fprintf(stderr, "Cluster with rep %s is adding new item named %s\n", result.names_[ids[cluster_id]].data(), result.names_[oid].data());)
        if(result.cardinalities_[cv.back()] > result.cardinalities_[rep]) {
            // In case the items are unsorted with respect to cardinality due to the parallelism, swap it out.
            // That way, we'll keep the highest-cardinality set as the representative
            std::swap(cv.back(), rep);
        }
    }
}

void cleanup(std::pair<std::vector<LSHIDType>, std::vector<std::vector<LSHIDType>>> &ret, sketch::lsh::SetSketchIndex<LSHIDType, LSHIDType> &retidx, const Dashing2DistOptions &opts, const SketchingResult &result, bool earlystop) {
    DBG_ONLY(std::fprintf(stderr, "%zu clusters before\n", ret.first.size());)
    DBG_ONLY(auto ts = std::chrono::high_resolution_clock::now();)
    std::vector<size_t> indicestorm;
    std::unique_ptr<std::mutex[]> locks(new std::mutex[ret.first.size()]);
    const LSHDistType mult = distance(opts.measure_) ? 1.: -1.;
    auto &constituents = ret.second;
    auto &ids= ret.first;
    OMP_PFOR
    for(size_t i = 0; i < ret.first.size(); ++i) {
        const auto oid = ret.first[i];
        assert(oid < result.names_.size());
        const minispan<RegT> span(&result.signatures_[opts.sketchsize_ * oid], opts.sketchsize_);
        auto [hits, counts, nper] = retidx.query_candidates(span, 5, size_t(-1), earlystop);
        for(size_t j = 0; j < hits.size(); ++j) {
            if(hits[j] == i) {
                std::swap(hits.back(), hits[j]);
                hits.pop_back();
                break;
            }
        }
        if(hits.empty()) continue;
        std::vector<LSHDistType> vals(hits.size());
        auto vp = vals.data();
        for(const auto id: hits) {
            assert(id < ids.size());
            *vp++ = mult * compare(opts, result, oid, ids[id]);
        }
        auto mv = std::min_element(vals.data(), vp);
        auto pos = mv - vals.data();
        auto cluster_id = hits[pos];
        auto &rep = ids[cluster_id];
        if(mv != vp && mult * *mv > opts.min_similarity_ && result.cardinalities_[oid] < result.cardinalities_[rep]) {
            auto cluster_id = hits.at(pos);
            auto &cv = constituents.at(cluster_id);
            {
                std::lock_guard<std::mutex> lock(locks[cluster_id]);
                cv.push_back(oid);
                cv.insert(cv.end(), ret.second[i].begin(), ret.second[i].end());
            }
            {
                std::lock_guard<std::mutex> lock2(locks[i]);
                ret.second[i].clear();
            }
            OMP_CRITICAL
            {
                indicestorm.push_back(i);
            }
        }
    }
    DBG_ONLY(std::fprintf(stderr, "Removing %zu indices out of %zu\n", indicestorm.size(), ret.first.size());)
    for(const auto id: indicestorm) {
        std::swap(ret.first[id], ret.first.back()), ret.first.pop_back();
        std::swap(ret.second[id], ret.second.back()), ret.second.pop_back();
    }
    DBG_ONLY(std::fprintf(stderr, "%zu clusters after %gs\n", ret.first.size(), std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - ts).count());)
#if 0
    ts = std::chrono::high_resolution_clock::now();
    indicestorm.clear();
    retidx.clear();
    for(size_t i = 0; i < ret.first.size() - 1; ++i) {
        const auto oid = ret.first[i];
        const size_t mysz = ret.first.size() - i - 1;
        std::unique_ptr<LSHDistType[]> dists(new LSHDistType[mysz]);
        #pragma omp parallel for schedule(static, 8)
        for(size_t j = 0; j < mysz; ++j) {
            dists[j] = mult * compare(opts, result, oid, ret.first[j + i + 1]);
        }
        auto mv = std::min_element(dists.get(), &dists[mysz]);
        auto pos = mv - dists.get();
        auto cluster_id = pos + i + 1;
        auto &rep = ids[cluster_id];
        //std::fprintf(stderr, "mv %g\n", *mv);
        if(mv != &dists[mysz] && *mv * mult > opts.min_similarity_ && result.cardinalities_[oid] < result.cardinalities_[rep]) {
            assert(pos + i + 1 < ret.first.size());
            auto &cv = constituents[cluster_id];
            std::lock_guard<std::mutex> lock(locks[cluster_id]);
            std::lock_guard<std::mutex> lock22(locks[i]);
            cv.push_back(oid);
            if(result.cardinalities_[cv.back()] > result.cardinalities_[rep])
                std::swap(cv.back(), rep);
            indicestorm.push_back(i);
            cv.insert(cv.end(), ret.second[i].begin(), ret.second[i].end());
            ret.second[i].clear();
        }
    }
    for(const auto id: indicestorm) {
        std::swap(ret.first[id], ret.first.back()), ret.first.pop_back();
        std::swap(ret.second[id], ret.second.back()), ret.second.pop_back();
    }
    std::fprintf(stderr, "After exhaustive removal, %zu clusters after %gs\n", ret.first.size(), std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - ts).count());
#endif
}


std::pair<std::vector<LSHIDType>, std::vector<std::vector<LSHIDType>>> dedup_core(sketch::lsh::SetSketchIndex<LSHIDType, LSHIDType> &retidx, const Dashing2DistOptions &opts, const SketchingResult &result)
{
    std::pair<std::vector<LSHIDType>, std::vector<std::vector<LSHIDType>>> ret;
    const size_t nelem = result.names_.size();
    std::unique_ptr<LSHIDType[]> order(new LSHIDType[nelem]);
#if _OPENMP > 201307L
    #pragma omp parallel for simd schedule(static, 1024)
#endif
    for(size_t i = 0; i < nelem; ++i) {
        order[i] = i;
    }
    std::sort(&order[0], &order[nelem], [&c=result.cardinalities_](auto x, auto y) {return c[x] > c[y];});
    int nt = 1;
#ifdef _OPENMP
    _Pragma("omp parallel")
    {nt = std::max(omp_get_num_threads(), 1);}
#endif
    // General strategy:
    // Use a given similarity threshold to then group items into the cluster
    // to which they are most similar if they are > than
    if(exhaustive_dedup) {
        const LSHDistType mult = distance(opts.measure_)  ? 1.: -1.;
        auto &ids = ret.first;
        auto &constituents = ret.second;
        for(size_t i = 0; i < nelem; ++i) {
            std::pair<LSHDistType, LSHIDType> bestc = {std::numeric_limits<LSHDistType>::max(), -1};
#ifdef _OPENMP
#pragma omp declare reduction(min: std::pair<LSHDistType, LSHIDType>: omp_out = std::min(omp_in, omp_out))
            #pragma omp parallel for schedule(dynamic) reduction(min:bestc)
#endif
            for(size_t j = 0; j < ids.size(); ++j) {
                bestc = std::min(bestc, std::pair<LSHDistType, LSHIDType>{compare(opts, result, i, ids[j]) * mult, j});
            }
            if(bestc.first * mult < opts.min_similarity_ || bestc.second == LSHIDType(-1)) {
                ids.push_back(i);
                constituents.emplace_back();
            } else {
                constituents.at(bestc.second).push_back(i);
            }
            /* exhaustive loading*/
        }
    } else {
#if 1
            auto &ids = ret.first;
            auto &constituents = ret.second;
            auto &idx = retidx;
            auto do_update = [&,st=nt<=1](auto id) __attribute__((always_inline)) {
                st ? update_res(id, ids, constituents, idx, opts, result)
                   : update_res_mt(id, ids, constituents, idx, opts, result);
            };
            for(size_t i = 0; i < nelem;do_update(order[i++]));
#else
        auto &idx = retidx;
        auto &ids = ret.first;
        auto &constituents = ret.second;
        if(nt <= 1) {
            for(size_t i = 0; i < nelem;) {
                 update_res(order[i++], ids, constituents, idx, opts, result);
            }
        } else {
#if 1
            const double simt = opts.min_similarity_ > 0. ? opts.min_similarity_: 0.9; // 90% is the default cut-off for deduplication
            using RetT = std::tuple<std::vector<LSHIDType>, std::vector<uint32_t>, std::vector<uint32_t>>;
            std::vector<RetT> batched_hits(nt);
            const size_t nbatches = (nelem + nt - 1) / nt;
            const LSHDistType mult = distance(opts.measure_) ? 1.: -1.;
            //std::deque<std::mutex> locks;
            std::mutex global_lock;
            for(size_t i = 0; i < nbatches; ++i) {
                const size_t start = nt * i, end = std::min(start + nt, nelem);
                #pragma omp parallel for
                for(size_t j = start; j < end; ++j) {
                    const auto oid = order[j];
                    const minispan<RegT> span(&result.signatures_[opts.sketchsize_ * oid], opts.sketchsize_);
                    const auto bhidx = j - start;
                    //std::fprintf(stderr, "bh9dx: %zu\n", bhidx);
                    auto &rettup = batched_hits.at(bhidx);
                    rettup = idx.query_candidates(span, MINCAND, size_t(-1), EARLYSTOP);
                    auto &[hits, counts, nper] = rettup;
                    std::vector<LSHDistType> vals;
                    typename std::vector<LSHDistType>::iterator mv;
                    if(!hits.empty()) {
                        //DBG_ONLY(std::fprintf(stderr, "Non-empty hits.\n");)
                        vals.resize(hits.size());
                        for(size_t i = 0, e = hits.size(); i < e; ++i) {
                            if(ids.size() < hits[i]) {
                                std::fprintf(stderr, "ids of size %zu yielded a hit index %zu at %zu. Skipping, but this shouldn't happen...\n", ids.size(), size_t(hits[i]), i);
                                continue;
                            }
                            vals[i] = mult * compare(opts, result, oid, ids.at(hits[i]));
                        }
                        mv = std::min_element(vals.begin(), vals.end());
                        if(mult * *mv < simt) {
                            auto cluster_id = hits.at(mv - vals.begin());
                            std::lock_guard<std::mutex> global(global_lock);
                            if(cluster_id >= constituents.size()) std::fprintf(stderr, "constit %zu hitting cluster_id %zu\n", constituents.size(), cluster_id);
                            auto &cv = constituents.at(cluster_id);
                            //std::lock_guard<std::mutex> lock(locks.at(cluster_id));
                            cv.push_back(oid);
                            auto &brep = cv.back();
                            auto &orep = ids[cluster_id];
                            if(result.cardinalities_[orep] < result.cardinalities_[brep]) {
                                DBG_ONLY(std::fprintf(stderr, "Swapping longer entity to be the representative. %g vs %g\n", result.cardinalities_[orep], result.cardinalities_[brep]);)
                                std::swap(orep, brep);
                            }
                            continue;
                        }
                    }
                    std::lock_guard<std::mutex> global(global_lock);
                    ids.push_back(oid);
                    constituents.emplace_back();
                    //locks.emplace_back();
                    idx.update(minispan<RegT>(&result.signatures_[opts.sketchsize_ * oid], opts.sketchsize_));
                    assert(idx.size() == constituents.size());
                    assert(idx.size() == ids.size());
                }
            }
#else
            std::vector<GreedyClustering> subs;
            retidx.unlock();
            subs.reserve(nt);
            while(subs.size() < unsigned(nt))
                subs.emplace_back(result, opts, retidx, true, MINCAND);
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic, 32)
#endif
            for(size_t i = 0; i < nelem; ++i) {
                const int tid = OMP_ELSE(omp_get_thread_num(), 0);
                auto &lres = subs[tid];
                update_res(order[i], lres.ids_, lres.constituents_, lres.idx_, opts, result);
            }
            par_reduce(subs.data(), subs.size());
            retidx = std::move(subs.front().idx_);
            ret = std::make_pair(std::move(subs.front().ids_), std::move(subs.front().constituents_));
            cleanup(ret, retidx, opts, result, true);
#endif
        }
#endif
    }
    return ret;
}
#ifndef NDEBUG
template<typename T, typename Alloc, typename VAlloc>
double medsize(const std::vector<std::vector<T, Alloc>, VAlloc> &v) {
    std::vector<size_t> sizes(v.size());
    for(size_t i = 0; i < v.size(); ++i)
        sizes[i] = v[i].size() + 1;
    std::sort(sizes.begin(), sizes.end());
    double ret;
    if(sizes.size() & 1)
        ret = sizes[sizes.size() / 2];
    else ret = .5 * sizes[sizes.size() / 2 - 1] + .5 * sizes[sizes.size() / 2];
    std::fprintf(stderr, "Median %g from %zu/%zu\n", ret,  sizes[sizes.size() / 2 - 1] , sizes[sizes.size() / 2]);
    return ret;
}
#endif

void dedup_emit(const std::vector<LSHIDType> &ids, const std::vector<std::vector<LSHIDType>> &constituents, const Dashing2DistOptions &opts, const SketchingResult &result) {
    const std::string &outname = opts.outfile_path_;
    std::FILE *ofp = stdout;
    if(outname.size() && (ofp = bfopen(outname.data(), "wb")) == nullptr) {
        THROW_EXCEPTION(std::runtime_error(std::string("Failed to open file ") + outname + " for writing"));
    }
    const size_t nclusters = ids.size();
    const size_t nitems = result.names_.size();
    const double avgsize = double(nitems) / nclusters;
    DBG_ONLY(const double medsz = medsize(constituents);)
    DBG_ONLY(std::fprintf(stderr, "#Clustering %zu items yielded %zu clusters of average size %0.4g (median size %0.4g), separated by minimum similarity %g\n", nitems, ids.size(), avgsize, medsz, opts.min_similarity_);)
    if(opts.output_format_ == HUMAN_READABLE) {
        fmt::print(ofp, "#Clustering {} items yielded {} clusters of average size {}, separated by minimum similarity {}\n", nitems, ids.size(), avgsize, opts.min_similarity_);
        for(size_t cid = 0;cid < ids.size(); ++cid) {
            auto repid = ids[cid];
            fmt::print(ofp, "Cluster-{}\t{}:{}", cid, result.names_[repid], repid);
            for(const auto child: constituents[cid]) {
                const size_t childid = child; // ids.at(child);
                fmt::print(ofp, "\t{}:{}", result.names_[childid], childid);
            }
            fmt::print(ofp, "\n");
        }
    } else {
        std::vector<uint64_t> indptr(nclusters + 1);
        for(size_t i = 0; i < nclusters; ++i) {
            indptr[i + 1] = indptr[i] + constituents[i].size() + 1; // 1 extra because the representative is not counted
        }
        const size_t nnz = indptr.back();
        const uint64_t dims [] {nclusters, nnz};
        std::vector<LSHIDType> indices(nnz);
        checked_fwrite(ofp, dims, sizeof(dims));
        checked_fwrite(ofp, indptr.data(), (indptr.size() * sizeof(uint64_t)));
        for(size_t i = 0; i < nclusters; ++i) {
            auto &v = constituents[i];
            checked_fwrite(ofp, &ids[i], sizeof(LSHIDType));
            checked_fwrite(ofp, v.data(), (v.size() * sizeof(LSHIDType)));
        }
    }
    if(ofp != stdout) std::fclose(ofp);
}
} // namespace dashing2
