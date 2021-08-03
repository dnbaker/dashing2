#include "cmp_main.h"
#include "sketch/ssi.h"
#include "minispan.h"

#define checked_fwrite(fp, ptr, nb) \
    do {\
        if(unsigned long long lrc = std::fwrite(static_cast<const void *>(ptr), 1, nb, fp); lrc != static_cast<size_t>(nb)) \
            throw std::runtime_error(std::string("[E:") + __PRETTY_FUNCTION__ + ':' + __FILE__ + std::to_string(__LINE__) + "] Failed to perform buffered write of " + std::to_string(static_cast<size_t>(nb)) + " bytes, instead writing " + std::to_string(lrc) + " bytes");\
    } while(0)

namespace dashing2 {
std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>> dedup_core(sketch::lsh::SetSketchIndex<LSHIDType, LSHIDType> &idx, Dashing2DistOptions &opts, const SketchingResult &result) {
    std::mutex mut; // Lock for IDs and constituents
    std::vector<size_t> order(result.names_.size());
    std::iota(order.begin(), order.end(), size_t(0));
    std::sort(order.begin(), order.end(), [&v=result.cardinalities_](auto x, auto y) {return v[x] < v[y];});
    const double simt = opts.min_similarity_ > 0. ? opts.min_similarity_: 0.9; // 90% is the default cut-off for deduplication
    std::vector<size_t> ids;
    std::vector<std::vector<size_t>> constituents;
    // General strategy:
    // Use a given similarity threshold to then group items into the cluster
    // to which they are most similar if they are > than
    OMP_PFOR_DYN
    for(size_t i = 0; i < order.size(); ++i) {
        auto oid = order[i];
        auto [hits, counts, nper] = idx.query_candidates(minispan<RegT>(&result.signatures_[opts.sketchsize_ * oid], opts.sketchsize_), 3);
        std::vector<LSHIDType> vals(hits.size());
        const LSHDistType mult = distance(opts.measure_) ? 1.: -1.;
        auto vp = vals.data();
        // These ids are indexes into the vector of results with ids/constiuents, so we access ids
        for(const auto id: hits) {
            const auto repid = ids[id];
            *vp++ = mult * compare(opts, result, oid, repid);
        }
        auto mv = std::min_element(vals.data(), vp);
        std::transform(vals.data(), vp, vals.data(), [mult](auto x) {return x * mult;});
        // auto v = compare(opts, result, lhid, id);
        if(hits.empty() || (mv != vp && *mv < simt)) {
            // The LSH index is thread-safe, so we only need to lock these.
            {
                std::lock_guard<std::mutex> lock(mut);
                ids.push_back(oid);
                constituents.emplace_back();
            }
            idx.update(minispan<RegT>(&result.signatures_[opts.sketchsize_ * oid], opts.sketchsize_));
            std::fprintf(stderr, "Added item %zu; %zu clusters so far (%%%0.4g)\n", idx.size(), ids.size(), double(ids.size()) / idx.size());
        } else {
            if(mv == vp) continue;
            auto pos = mv - &vals[0];
            auto cluster_id = hits[pos];
            std::lock_guard<std::mutex> lock(mut);
            auto &cv = constituents[cluster_id];
            cv.push_back(oid);
            if(result.cardinalities_[cv.back()] > result.cardinalities_[ids[cluster_id]]) {
                // In case the items are unsorted with respect to cardinality due to the parallelism, swap it out.
                // That way, we'll keep the highest-cardinality set as the representative
                std::swap(cv.back(), ids[cluster_id]);
            }
        }
    }
    return std::make_pair(ids, constituents);
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
            for(const auto child: constituents[cid])
                std::fprintf(ofp, "\t%s:%zu", result.names_[ids[child]].data(), ids[child]);

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
