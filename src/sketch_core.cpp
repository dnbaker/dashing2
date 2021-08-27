#include "sketch_core.h"

namespace dashing2 {
extern size_t MEMSIGTHRESH;

INLINE size_t nbytes_from_line(const std::string &line) {
    size_t ret = 0;
    for_each_substr([&ret](const std::string &s) {ret += bns::filesize(s.data());}, line);
    return ret;
}

SketchingResult &sketch_core(SketchingResult &result, Dashing2Options &opts, const std::vector<std::string> &paths, std::string &outfile) {
    if(opts.kmer_result() == FULL_MMER_SEQUENCE && outfile.empty()) {
        THROW_EXCEPTION(std::runtime_error("outfile must be specified for --seq mode."));
    }
    result.signatures_.memthreshold(MEMSIGTHRESH);
    const size_t npaths = paths.size();
    std::string tmpfile;
    bool rmfile = false;
    if(opts.dtype_ == DataType::FASTX) {
        if(opts.parse_by_seq_) {
            if(paths.size() != 1) {
                result.nperfile_.resize(paths.size());
                THROW_EXCEPTION(std::runtime_error("parse-by-seq currently only handles one file at a time. To process multiple files, simply concatenate them into one file, and run dashing2 on that."));
            }
            KSeqHolder kseqs(std::max(opts.nthreads(), 1u));
            fastx2sketch_byseq(result, opts, paths.front(), kseqs.kseqs_, outfile, true);
        } else {
            fastx2sketch(result, opts, paths, outfile);
        }
    } else if(opts.dtype_ == DataType::LEAFCUTTER) {
        if(outfile.empty()) THROW_EXCEPTION(std::runtime_error("Outfile required for LeafCutter input."));
        auto res = lf2sketch(paths, opts);
        result.names_ = std::move(res.sample_names());
        result.nperfile_.resize(res.nsamples_per_file().size());
        std::copy(res.nsamples_per_file().begin(), res.nsamples_per_file().end(), result.nperfile_.begin());
        std::FILE *of = bfopen(outfile.data(), "wb");
        uint64_t ns = result.names_.size();
        std::fwrite(&ns, sizeof(ns), 1, of);
        ns = opts.sketchsize_;
        std::fwrite(&ns, sizeof(ns), 1, of);
        result.signatures_.assign(outfile);
        result.signatures_.resize(res.registers().size());
        std::copy(res.registers().begin(), res.registers().end(), result.signatures_.begin());
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
                OMP_PFOR_DYN
                for(size_t i = 0; i < npaths; ++i) {
                    auto myind = filesizes.size() ? filesizes[i].second: uint64_t(i);
                    auto &p(paths[myind]);
                    result.names_[i] = p;
                    auto res = bw2sketch(p, opts, /*parallel_process=*/false);
                    std::copy(res.global_->begin(), res.global_->end(), &result.signatures_[myind * opts.sketchsize_]);
                    result.cardinalities_[myind] = res.card_;
                }
            }
        }
    }
    std::FILE *ofp;
    if(opts.kmer_result_ == FULL_MMER_SEQUENCE) {
        if((ofp = bfopen(outfile.data(), "r+")) == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to open output file for mmer sequence results."));
        size_t offset = result.names_.size();
        std::fwrite(&offset, sizeof(offset), 1, ofp);
        {
            const uint32_t k = opts.k_, w = opts.w_;
            std::fwrite(&k, sizeof(k), 1, ofp);
            std::fwrite(&w, sizeof(w), 1, ofp);
        }
        if(std::fwrite(result.cardinalities_.data(), sizeof(double), result.cardinalities_.size(), ofp) != result.cardinalities_.size())
            THROW_EXCEPTION(std::runtime_error("Failed to write lengths of sequences to disk.\n"));
        offset = 0;
        for(size_t i = 0; i < result.nperfile_.size(); ++i) {
            if(result.nperfile_[i]) {
                std::fwrite(&result.signatures_.at(offset), sizeof(RegT), result.nperfile_[i], ofp);
                offset += result.nperfile_.at(i);
            }
        }
        std::fclose(ofp);
    } else {
        if(outfile.size() && outfile != "/dev/stdout" && outfile != "-") {
            // This should not overlap with the memory mapped for result.signatures_
            const uint64_t t = result.cardinalities_.size();
            const size_t nb = t * sizeof(double) + 2 * sizeof(uint64_t);
            int fd = ::open(outfile.data(), O_RDWR);
            void *tmpptr = ::mmap(nullptr, nb, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            ::close(fd);
            if(!tmpptr) {
                THROW_EXCEPTION(std::runtime_error("Failed to set metadata for final stacked sketches."));
            }
            ((uint64_t *)tmpptr)[0] = t;
            ((uint64_t *)tmpptr)[1] = opts.sketchsize_;
            std::copy(result.cardinalities_.begin(), result.cardinalities_.end(), ((double *)tmpptr) + 2);
            ::munmap(tmpptr, nb);
        } else {
            if(outfile.size() && (outfile == "-" || outfile == "/dev/stdout")) {
                THROW_EXCEPTION(std::runtime_error("Not yet supported: writing stacked sketches to file streams. This may change."));
            }
        }
    }
    if(result.names_.size()) {
        if((ofp = bfopen((outfile + ".names.txt").data(), "wb")) == nullptr)
            THROW_EXCEPTION(std::runtime_error(std::string("Failed to open outfile at ") + outfile + ".names.txt"));
        std::fputs("#Name\tCardinality\n", ofp);
        for(size_t i = 0; i < result.names_.size(); ++i) {
            const auto &n(result.names_[i]);
            if(std::fwrite(n.data(), 1, n.size(), ofp) != n.size()) THROW_EXCEPTION(std::runtime_error("Failed to write names to file"));
            if(result.cardinalities_.size()) {
                std::fprintf(ofp, tfmt<double>, result.cardinalities_[i]);
            }
            if(auto &kf = result.kmercountfiles_; !kf.empty())
                std::fputc('\t', ofp), std::fwrite(kf[i].data(), 1, kf[i].size(), ofp);
            std::fputc('\n', ofp);
        }
        std::fclose(ofp);
    }
    if(result.kmers_.size()) {
        const size_t nb = result.kmers_.size() * sizeof(uint64_t);
        DBG_ONLY(std::fprintf(stderr, "About to write result kmers to %s of size %zu\n", outfile.data(), nb);)
        if((ofp = bfopen((outfile + ".kmerhashes.u64").data(), "wb"))) {
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
        if((ofp = bfopen((outfile + ".kmercounts.f64").data(), "wb"))) {
            if(std::fwrite(result.kmercounts_.data(), nb, 1, ofp) != size_t(1))
                std::fprintf(stderr, "Failed to write k-mer counts to disk properly, silent failing\n");
            std::fclose(ofp);
        } else {
            DBG_ONLY(std::fprintf(stderr, "Failed to open file at %s to write k-mer counts, failing silently.\n", (outfile + ".kmercounts.f64").data());)
        }
    }
    if(rmfile) {
        const int rc = std::system(("rm "s + outfile).data());
        if(rc)
            std::fprintf(stderr, "Failed to rm %s, but this is perhaps unimportant. rc %d, exit status %d, and signal %d\n", outfile.data(), rc, WEXITSTATUS(rc), WSTOPSIG(rc));
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
