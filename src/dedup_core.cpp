#include "cmp_main.h"
#include "src/ssi.h"
#include "minispan.h"
#include "fastiota.h"
#ifdef _OPENMP
#include "omp.h"
#endif

#define checked_fwrite(fp, ptr, nb) \
    do {\
        if(unsigned long long lrc = std::fwrite(static_cast<const void *>(ptr), 1, nb, fp); lrc != static_cast<size_t>(nb)) \
            throw std::runtime_error(std::string("[E:") + __PRETTY_FUNCTION__ + ':' + __FILE__ + std::to_string(__LINE__) + "] Failed to perform buffered write of " + std::to_string(static_cast<size_t>(nb)) + " bytes, instead writing " + std::to_string(lrc) + " bytes");\
    } while(0)

namespace dashing2 {

static constexpr size_t MINCAND = 10;

struct GreedyClustering {
    std::vector<size_t> ids_;
    std::vector<std::vector<size_t>> constituents_;
    sketch::lsh::SetSketchIndex<LSHIDType, LSHIDType> idx_;
    double simt, mult;
    bool earlystop;
    size_t mincand;
    const SketchingResult &result;
    Dashing2DistOptions &opts;
    size_t sketchsize_;
    GreedyClustering(const SketchingResult &rs, Dashing2DistOptions &opts, const sketch::lsh::SetSketchIndex<LSHIDType, LSHIDType> &idx, bool earlystop=false, size_t mincand=MINCAND)
        : idx_(idx.clone()),
        earlystop(earlystop), mincand(mincand), result(rs), opts(opts)
    {
        simt = opts.min_similarity_ > 0. ? opts.min_similarity_: 0.9; // 90% is the default cut-off for deduplication
        mult = distance(opts.measure_)  ? 1.: -1.;
        assert(sketchsize > 0);
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
            } else {
                cv.push_back(orep);
            }
            cv.insert(cv.end(), ocon.begin(), ocon.end());
        }
        return *this;
    }
};

template<typename T>
void par_reduce(T *x, size_t n) {
    const unsigned int ln = static_cast<int>(std::ceil(std::log2(n)));
    for(size_t i = 0; i < ln; ++i) {
        const size_t step_size = 1 << i;
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



void update_res(size_t oid, std::vector<size_t> &ids, std::vector<std::vector<size_t>> &constituents,
                sketch::lsh::SetSketchIndex<LSHIDType, LSHIDType> &idx, Dashing2DistOptions &opts, const SketchingResult &result,
                bool earlystop=false, size_t mincand=MINCAND)
{
    const double simt = opts.min_similarity_ > 0. ? opts.min_similarity_: 0.9; // 90% is the default cut-off for deduplication
    const minispan<RegT> span(&result.signatures_[opts.sketchsize_ * oid], opts.sketchsize_);
    auto [hits, counts, nper] = idx.query_candidates(span, mincand, size_t(-1), earlystop);
    std::vector<LSHDistType> vals(hits.size());
    auto vp = vals.data();
    const LSHDistType mult = distance(opts.measure_) ? 1.: -1.;
    for(const auto id: hits) {
        *vp++ = mult * compare(opts, result, oid, ids[id]);
    }
    auto mv = std::min_element(vals.data(), vp);
    if(hits.empty() || (mv != vp && mult * *mv < simt)) {
#if 0
        if(mv != vp && mult * *mv < simt)
            std::fprintf(stderr, "Max similarity is %g vs %g vs %g\n", mult * *mv, simt, *std::max_element(vals.data(), vp));
#endif
        ids.push_back(oid);
        DBG_ONLY(size_t nids = ids.size();)
        constituents.emplace_back();
        DBG_ONLY(size_t myid = )
            idx.update(minispan<RegT>(&result.signatures_[opts.sketchsize_ * oid], opts.sketchsize_));
        DBG_ONLY(std::fprintf(stderr, "Added item %zu/%s; %zu clusters so far (%%%0.4g)\n", myid, result.names_[oid].data(), nids, ids.size() * 100. / order.size());)
    } else {
        if(mv == vp) return;
        auto pos = mv - &vals[0];
        auto cluster_id = hits[pos];
        auto &cv = constituents[cluster_id];
        auto &rep = ids[cluster_id];
        cv.push_back(oid);
        DBG_ONLY(std::fprintf(stderr, "Cluster with rep %s is adding new item named %s\n", result.names_[ids[cluster_id]].data(), result.names_[oid].data());)
        if(result.cardinalities_[cv.back()] > result.cardinalities_[rep]) {
            // In case the items are unsorted with respect to cardinality due to the parallelism, swap it out.
            // That way, we'll keep the highest-cardinality set as the representative
            //std::fprintf(stderr, "swapped in %g for %g\n", result.cardinalities_[cv.back()], result.cardinalities_[ids[pos]]);
            if(result.cardinalities_[cv.back()] > result.cardinalities_[rep])
                std::swap(cv.back(), rep);
            //std::fprintf(stderr, "swapped out %g for %g\n", result.cardinalities_[cv.back()], result.cardinalities_[ids[pos]]);
        }
    }
}


std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>> dedup_core(sketch::lsh::SetSketchIndex<LSHIDType, LSHIDType> &retidx, Dashing2DistOptions &opts, const SketchingResult &result) {
    std::vector<size_t> order(result.names_.size());
    fastiota::iota(order.data(), order.size());
    assert(std::all_of(order.begin(), order.end(), [&](auto &x) {return x == uint64_t(&x - order.data());}));
    std::sort(order.begin(), order.end(), [&v=result.cardinalities_](auto x, auto y) {return v[x] < v[y];});

    int nt = 1;
#ifdef _OPENMP
    _Pragma("omp parallel")
    {nt = omp_get_num_threads();}
#endif
    // General strategy:
    // Use a given similarity threshold to then group items into the cluster
    // to which they are most similar if they are > than
    if(nt == 1) {
        std::vector<size_t> ids;
        std::vector<std::vector<size_t>> constituents;
        auto &idx = retidx;
        for(size_t i = 0; i < order.size(); ++i) {
            update_res(order[i], ids, constituents, idx, opts, result);
        }
        return std::make_pair(ids, constituents);
    } else {
        std::vector<GreedyClustering> subs;
        subs.reserve(nt);
        while(subs.size() < unsigned(nt))
            subs.emplace_back(result, opts, retidx, false, MINCAND);
        OMP_PFOR
        for(size_t i = 0; i < order.size(); ++i) {
            const int tid = OMP_ELSE(omp_get_num_threads(), 0);
            auto &lres = subs[tid];
            auto oid = order[i];
            update_res(oid, lres.ids_, lres.constituents_, lres.idx_, opts, result);
        }
        par_reduce(subs.data(), subs.size());
        retidx = std::move(subs.front().idx_);
        return std::make_pair(std::move(subs.front().ids_), std::move(subs.front().constituents_));
    }
}

void dedup_emit(const std::vector<size_t> &ids, const std::vector<std::vector<size_t>> &constituents, const Dashing2DistOptions &opts, const SketchingResult &result) {
    const std::string &outname = opts.outfile_path_;
    std::FILE *ofp = stdout;
    if(outname.size() && (ofp = std::fopen(outname.data(), "wb")) == nullptr) {
        THROW_EXCEPTION(std::runtime_error(std::string("Failed to open file ") + outname + " for writing"));
    }
    const size_t nclusters = ids.size();
    std::fprintf(stderr, "#%zu clusters of average size %0.4g, separated by minimum similarity %g\n", ids.size(), double(ids.size()) / result.names_.size(), opts.min_similarity_);
    if(opts.output_format_ == HUMAN_READABLE) {
        std::fprintf(ofp, "#%zu clusters of average size %0.4g, separated by minimum similarity %g\n", ids.size(), double(ids.size()) / result.names_.size(), opts.min_similarity_);
        for(size_t cid = 0;cid < ids.size(); ++cid) {
            std::fprintf(ofp, "Cluster-%zu\t%s:%zu", cid, result.names_[ids[cid]].data(), ids[cid]);
            for(const auto child: constituents[cid]) {
                size_t childid = child; // ids.at(child);
                std::fprintf(ofp, "\t%s:%zu", result.names_[childid].data(), childid);
            }
            std::fputc('\n', ofp);
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
            checked_fwrite(ofp, &ids[i], sizeof(std::decay_t<decltype(ids[i])>));
            checked_fwrite(ofp, v.data(), (v.size() * sizeof(size_t)));
        }
    }
    if(ofp != stdout) std::fclose(ofp);
}
} // namespace dashing2
