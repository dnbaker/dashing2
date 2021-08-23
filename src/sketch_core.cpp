#include "sketch_core.h"

namespace dashing2 {

INLINE size_t nbytes_from_line(const std::string &line) {
    size_t ret = 0;
    for_each_substr([&ret](const std::string &s) {ret += bns::filesize(s.data());}, line);
    return ret;
}

SketchingResult sketch_core(Dashing2Options &opts, const std::vector<std::string> &paths, std::string &outfile) {
    SketchingResult result;
    const size_t npaths = paths.size();
    if(opts.dtype_ == DataType::FASTX) {
        result = fastx2sketch(opts, paths);
    } else if(opts.dtype_ == DataType::LEAFCUTTER) {
        auto res = lf2sketch(paths, opts);
        result.signatures_ = std::move(res.registers());
        result.names_ = std::move(res.sample_names());
        result.nperfile_.resize(res.nsamples_per_file().size());
        std::copy(res.nsamples_per_file().begin(), res.nsamples_per_file().end(), result.nperfile_.begin());
    } else if(opts.dtype_ == DataType::BED || opts.dtype_ == DataType::BIGWIG) {
        std::vector<std::pair<size_t, uint64_t>> filesizes = get_filesizes(paths);
        result.signatures_.resize(npaths * opts.sketchsize_);
        result.names_.resize(npaths);
        result.cardinalities_.resize(npaths);
        if(opts.dtype_ == DataType::BED) {
            OMP_PFOR_DYN
            for(size_t i = 0; i < npaths; ++i) {
                auto myind = filesizes.size() ? filesizes[i].second: uint64_t(i);
                auto &p(paths[myind]);
                result.names_[i] = p;
                auto [sig, card] = bed2sketch(p, opts);
                result.cardinalities_[myind] = card;
                std::copy(sig.begin(), sig.end(), &result.signatures_[myind * opts.sketchsize_]);
            }
        } else {
            // BigWig sketching is parallelized within files
            if(opts.by_chrom_) {
                std::vector<flat_hash_map<std::string, std::vector<RegT>>> bc(npaths);
                std::vector<flat_hash_map<std::string, double>> dbc(npaths);
                result.nperfile_.resize(npaths);
                for(size_t i = 0; i < npaths; ++i) {
                    auto &p = paths[i];
                    auto res = bw2sketch(p, opts, /*parallel_process=*/true);
                    result.names_[i] = p;
                    result.cardinalities_[i] = res.card_;
                    result.nperfile_[i] = res.chrmap_->size();
                    bc[i] = std::move(*res.chrmap_.get());
                    dbc[i] = std::move(res.cardmap_);
                    DBG_ONLY(std::fprintf(stderr, "Cardinality %g found from path %s/%zu\n", res.card_, p.data(), i);)
                    //std::copy(sigs.begin(), sigs.end(), &result.signatures_[opts.sketchsize_* i]);
                }
                const auto total_n = std::accumulate(result.nperfile_.begin(), result.nperfile_.end(), size_t(0));
                size_t offset = 0;
                result.signatures_.resize(total_n * opts.sketchsize_);
                result.names_.resize(total_n);
                result.cardinalities_.resize(total_n);
                for(size_t i = 0; i < npaths; ++i) {
                    size_t myi = 0;
                    for(const auto &pair: bc[i]) {
                        const auto om = offset + myi++;
                        std::copy(pair.second.begin(), pair.second.end(), &result.signatures_[opts.sketchsize_ * om]);
                        result.names_[om] = paths[i] + ":" + pair.first;
                        result.cardinalities_[om] = dbc[i][pair.first];
                    }
                    offset += bc[i].size();
                }
            } else {
                for(size_t i = 0; i < npaths; ++i) {
                    auto myind = filesizes.size() ? filesizes[i].second: uint64_t(i);
                    auto &p(paths[myind]);
                    result.names_[i] = p;
                    auto res = bw2sketch(p, opts, /*parallel_process=*/true);
                    std::copy(res.global_->begin(), res.global_->end(), &result.signatures_[myind * opts.sketchsize_]);
                    result.cardinalities_[myind] = res.card_;
                }
            }
        }
    }
    if(paths.size() == 1 && outfile.empty()) {
        const std::string suf =
                opts.sspace_ == SPACE_SET ? (opts.kmer_result_ == ONE_PERM ? ".opss": ".ss"):
                opts.sspace_ == SPACE_MULTISET ? ".bmh":
                opts.sspace_ == SPACE_PSET ? ".pmh" :
                opts.sspace_ == SPACE_EDIT_DISTANCE ? ".omh": ".unknown_sketch";
        outfile = paths.front();
        outfile = outfile.substr(0, outfile.find_first_of(' '));
        // In case the first path has multiple entries, trim to just the first
        outfile = outfile + suf;
        if(opts.trim_folder_paths()) {
            outfile = trim_folder(outfile);
            if(opts.outprefix_.size())
                outfile = opts.outprefix_ + '/' + outfile;
        }
    }
    bool even = (opts.kmer_result_ != FULL_MMER_SEQUENCE && (result.nperfile_.empty() || std::all_of(result.nperfile_.begin() + 1, result.nperfile_.end(), [v=result.nperfile_.front()](auto x) {return x == v;})));
    if(outfile.size()) {
        std::fprintf(stderr, "outfile %s\n", outfile.data());
        if(result.signatures_.empty()) THROW_EXCEPTION(std::runtime_error("Can't write stacked sketches if signatures were not generated"));
        std::fprintf(stderr, "Writing stacked sketches to %s\n", outfile.data());
        std::FILE *ofp;
        if(result.signatures_.size()) {
            ofp = std::fopen(outfile.data(), "wb");
            if(!ofp) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open file at ") + outfile));
            uint64_t n = result.cardinalities_.size();
            std::fwrite(&n, sizeof(n), 1, ofp);
            n = opts.sketchsize_;
            std::fwrite(&n, sizeof(n), 1, ofp);
            std::fwrite(result.cardinalities_.data(), sizeof(result.cardinalities_.front()), result.cardinalities_.size(), ofp);
            if(opts.kmer_result_ > FULL_SETSKETCH || even) {
                checked_fwrite(ofp, result.signatures_.data(), sizeof(RegT) * result.signatures_.size());
            } else {
                DBG_ONLY(std::fprintf(stderr, "%zu total summaries, accumulated into signatures of size %zu\n", std::accumulate(result.nperfile_.begin(), result.nperfile_.end(), size_t(0)), result.signatures_.size());)
                size_t offset = 0;
                const uint64_t terminus = uint64_t(-1);
                for(size_t i = 0; i < result.nperfile_.size(); ++i) {
                    if(result.nperfile_[i]) {
                        std::fwrite(&result.signatures_.at(offset), sizeof(RegT), result.nperfile_[i], ofp);
                        offset += result.nperfile_.at(i);
                    }
                    std::fwrite(&terminus, sizeof(terminus), 1, ofp);
                }
            }
            std::fclose(ofp);
        } DBG_ONLY(else std::fprintf(stderr, "Signatures is empty\n");)
        if(result.names_.size()) {
            if((ofp = std::fopen((outfile + ".names.txt").data(), "wb")) == nullptr)
                THROW_EXCEPTION(std::runtime_error(std::string("Failed to open outfile at ") + outfile + ".names.txt"));
            std::fputs("#Name\tCardinality\n", ofp);
            for(size_t i = 0; i < result.names_.size(); ++i) {
                const auto &n(result.names_[i]);
                if(std::fwrite(n.data(), 1, n.size(), ofp) != n.size()) THROW_EXCEPTION(std::runtime_error("Failed to write names to file"));
                if(result.cardinalities_.size()) {
                    std::fprintf(ofp, "\t%0.12g", result.cardinalities_[i]);
                }
                if(auto &kf = result.kmercountfiles_; kf.size())
                    std::fputc('\t', ofp), std::fwrite(kf[i].data(), 1, kf[i].size(), ofp);
                std::fputc('\n', ofp);
            }
            std::fclose(ofp);
        }
        if(result.kmers_.size()) {
            const size_t nb = result.kmers_.size() * sizeof(uint64_t);
            DBG_ONLY(std::fprintf(stderr, "About to write result kmers to %s of size %zu\n", outfile.data(), nb);)
            if((ofp = std::fopen((outfile + ".kmerhashes.u64").data(), "wb"))) {
                if(std::fwrite(result.kmers_.data(), nb, 1, ofp) != size_t(1))
                    std::fprintf(stderr, "Failed to write k-mers to disk properly, silent failing\n");
                std::fclose(ofp);
            } else {
                DBG_ONLY(std::fprintf(stderr, "Failed to write k-mers, failing silently.\n");)
            }
        } DBG_ONLY(else std::fprintf(stderr, "Not saving k-mers because result kmers is empty\n");)
        if(result.kmercounts_.size()) {
            const size_t nb = result.kmercounts_.size() * sizeof(decltype(result.kmercounts_)::value_type);
            DBG_ONLY(std::fprintf(stderr, "Writing kmercounts of size %zu\n", result.kmercounts_.size());)
            if((ofp = std::fopen((outfile + ".kmercounts.f64").data(), "wb"))) {
                if(std::fwrite(result.kmercounts_.data(), nb, 1, ofp) != size_t(1))
                    std::fprintf(stderr, "Failed to write k-mer counts to disk properly, silent failing\n");
                std::fclose(ofp);
            } else {
                DBG_ONLY(std::fprintf(stderr, "Failed to open file at %s to write k-mer counts, failing silently.\n", (outfile + ".kmercounts.f64").data());)
            }
        }
    }
    return result;
}

std::vector<std::pair<size_t, uint64_t>> get_filesizes(const std::vector<std::string> &paths) {
    const size_t npaths = paths.size();
    std::vector<std::pair<size_t, uint64_t>> filesizes(npaths);
    OMP_PFOR
    for(size_t i = 0; i < npaths; ++i) {
        filesizes[i] = {nbytes_from_line(paths[i]), uint64_t(i)};
    }
    std::sort(filesizes.begin(), filesizes.end(), std::greater<>());
    return filesizes;
}

}
