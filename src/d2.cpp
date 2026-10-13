#include "d2.h"
#include <filesystem>

namespace dashing2 {
int cmp_main(int argc, char **argv);
int contain_main(int argc, char **argv);
int wsketch_main(int argc, char **argv);
int sketch_main(int argc, char **argv);
int printmin_main(int argc, char **argv);
std::string Dashing2Options::to_string() const {
    size_t m = 4096;
    std::string ret(m, '\0');
    auto pos = std::sprintf(&ret[0], "Dashing2Options;k:%d", k_);
    if(w_ > 0) pos += std::sprintf(&ret[pos], ";w:%d", w_);
    pos += std::sprintf(&ret[pos], ";%s", parse_by_seq_ ? "parsebyseq": "parsebyfile");
    if(trim_chr_) pos += std::sprintf(&ret[pos], ";trimchr");
    pos += std::sprintf(&ret[pos], ";sketchsize:%zu", sketchsize_);
    if(count_threshold_ > 0)
        pos += std::sprintf(&ret[pos], ";%u", count_threshold_);
    pos += std::sprintf(&ret[pos], ";sketchtype:%s",
            kmer_result_ == ONE_PERM ? "onepermsetsketch"
                  : kmer_result_ == FULL_SETSKETCH ? (sspace_ == SPACE_SET ? "fullsetsketch": sspace_ == SPACE_MULTISET ? "bagminhash": sspace_ == SPACE_PSET ? "probminhash": sspace_ == SPACE_EDIT_DISTANCE ? "orderminhash": "unknown")
                  : kmer_result_ == FULL_MMER_SEQUENCE ? (use128() ? "mmerseq128": "mmerseq64")
                  : (kmer_result_ == FULL_MMER_SET || kmer_result_ == FULL_MMER_COUNTDICT) ? (use128() ? "mmerset128": "mmerset64"): "fullyunknown"
    );
    if(kmer_result_ == FULL_MMER_COUNTDICT)
        pos += std::sprintf(&ret[pos], ",kmercountsf64");
    pos += std::sprintf(&ret[pos], ";%s", ::dashing2::to_string(dtype_).data());
    if(spacing_.size()) {
        pos += std::sprintf(&ret[pos], ";spacing:%s", spacing_.data());
    }
    if(outprefix_.size()) pos += std::sprintf(&ret[pos], ";outprefix:%s", outprefix_.data());
    if(cssize_) pos += std::sprintf(&ret[pos], ";counting=countsketch%zu\n", cssize_);
    if(bed_parse_normalize_intervals_) pos += std::sprintf(&ret[pos], ";normalize_intervals");
    if(by_chrom_) pos += std::sprintf(&ret[pos], ";sketchbychrom");
    if(homopolymer_compress_minimizers_) pos += std::sprintf(&ret[pos], ";hp-compress-minimizers");
    if(canonicalize()) pos += std::sprintf(&ret[pos], ";canon");
    if(fs_) {
        pos += std::sprintf(&ret[pos], ";%s", fs_->to_string().data());
    }
    ret.resize(pos);
    return ret;
}

void Dashing2Options::filterset(const std::string &fsarg) {
    if(fsarg.empty()) return;
    auto i = fsarg.find_last_of(':');
    filterset(fsarg.substr(0, i), (i != std::string::npos && static_cast<char>(fsarg[i + 1] & 0xdf) != 'K'));
}

void Dashing2Options::filterset(const std::string &path, bool is_kmer) {
    fs_.reset(new FilterSet());
    if(is_kmer) {
        std::string cmd;
        if(endswith(path, ".xz")) cmd = "xz -dc ";
        else if(endswith(path, ".gz")) cmd = "gzip -dc ";
        else if(endswith(path, ".bz2")) cmd = "bzip2 -dc ";
        else if(endswith(path, ".zst")) cmd = "zstd -dc ";
        else cmd = "";
        std::FILE *ifp;
        if(cmd.empty()) {
            if((ifp = bfopen(cmd.data(), "rb")) == 0)
                THROW_EXCEPTION(std::runtime_error("Failed to open file "s + path + " for reading"));
        } else {
            if((ifp= ::popen((cmd + path).data(), "r")) == 0)
                THROW_EXCEPTION(std::runtime_error("Failed to run command "s + cmd + path));
        }
        uint64_t u6; u128_t u12;
        void *const ptr = use128() ? (void *)&u12: (void *)&u6;
        for(const auto is(use128() ? 16: 8);std::fread(ptr, is, 1, ifp) == 1;) {
            if(use128()) fs_->add(u12);
            else         fs_->add(u6);
        }
        if(cmd.empty())
            std::fclose(ifp);
        else
            ::pclose(ifp);
    } else {
        for_each_substr([&](const std::string &subpath) {
            //std::fprintf(stderr, "Doing for_each_substr for subpath = %s\n", subpath.data());
            const auto sp = subpath.data();
            auto lfunc = [&](auto x) {fs_->add(maskfn(x));};
            if(use128()) {
                if(unsigned(k_) <= nremperres128()) {
                    auto encoder = make_encoder<bns::score::Lex, u128_t>(sp_, rht_, canonicalize());
                    encoder.for_each(lfunc, sp);
                } else {
                    rh128_.for_each_hash(lfunc, sp);
                }
            } else if(unsigned(k_) <= nremperres64()) {
                enc_.value().for_each(lfunc, sp);
            } else {
                rh_.for_each_hash(lfunc, sp);
            }
        }, path);
    }
    fs_->finalize();
}
void Dashing2Options::validate() const {
    if(canonicalize() && rh_.hashtype() != bns::DNA) {
        canonicalize(false);
    }
    if((dtype_ == BED || dtype_ == BIGWIG) && kmer_result_ >= FULL_MMER_SET) {
        // Exact k-mer set, count dictionary, and minimizer sequence modes are
        // compared from per-input k-mer files, which interval sketching does
        // not produce; comparing would read whatever the cache holds instead.
        THROW_EXCEPTION(std::invalid_argument("Exact k-mer set, count dictionary, and minimizer sequence modes are only available for sequence (FASTX) inputs. Use a sketching mode for BED and BigWig data."));
    }
    const int lim = use128() ? nremperres128() : nremperres64();
    if(entmin && (k_ > lim || !sp_.unspaced())) {
        std::fprintf(stderr, "Warning: entropy minimization is not supported for spaced seeds or rolling hashers; falling back to random minimization.\n");
        entmin = false;
    }
}

bool entmin = false;
int verbosity = 0;

} // dashing2
