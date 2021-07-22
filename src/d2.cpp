#include "d2.h"

namespace dashing2 {
int cmp_main(int argc, char **argv);
int wsketch_main(int argc, char **argv);
int sketch_main(int argc, char **argv);
std::string Dashing2Options::to_string() const {
    size_t m = 4096;
    std::string ret(m, '\0');
    auto pos = std::sprintf(&ret[0], "Dashing2Options;k:%d", k_);
    if(w_ > 0) pos += std::sprintf(&ret[pos], ";w:%d", w_);
    pos += std::sprintf(&ret[pos], ";%s", parse_by_seq_ ? "parsebyseq": "parsebyfile");
    if(trim_chr_) pos += std::sprintf(&ret[pos], ";trimchr");
    pos += std::sprintf(&ret[pos], ";sketchsize:%zu", sketchsize_);
    if(count_threshold_ > 0)
        pos += std::sprintf(&ret[pos], ";%0.8g", count_threshold_);
    pos += std::sprintf(&ret[pos], ";sketchtype:%s",
            kmer_result_ == ONE_PERM ? "onepermsetsketch"
                  : kmer_result_ == FULL_SETSKETCH ? (sspace_ == SPACE_SET ? "fullsetsketch": sspace_ == SPACE_MULTISET ? "bagminhash": sspace_ == SPACE_PSET ? "probminhash": sspace_ == SPACE_EDIT_DISTANCE ? "orderminhash": "unknown")
                  : kmer_result_ == FULL_MMER_SEQUENCE ? (use128() ? "mmerseq128": "mmerseq64")
                  : kmer_result_ == FULL_MMER_SET ? (use128() ? "mmerset128": "mmerset64"): "fullyunknown"
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
    ret += ";command:\"";
    ret += cmd_;
    ret += "\"";
    return ret;
}
void Dashing2Options::filterset(std::string path, bool is_kmer) {
    if(is_kmer) {
        std::string cmd;
        if(endswith(path, ".xz")) cmd = "xz -dc ";
        else if(endswith(path, ".gz")) cmd = "gzip -dc ";
        else if(endswith(path, ".bz2")) cmd = "bzip2 -dc ";
        else cmd = "cat ";
        cmd += path;
        std::FILE *ifp = ::popen(cmd.data(), "r");
        fs_.reset(new FilterSet());
        union {
            uint64_t u6;
            u128_t u12;
        } u;
        for(const auto is(use128() ? 16: 8);std::fread(&u, is, 1, ifp) == 1;) {
            if(use128()) fs_->add(u.u12);
            else         fs_->add(u.u6);
        }
        ::pclose(ifp);
        return;
    }
    auto perf_for_substrs = [&](const auto &func) {
        for_each_substr([&](const std::string &subpath) {
            //std::fprintf(stderr, "Doing for_each_substr for subpath = %s\n", subpath.data());
            auto lfunc = [&](auto x) {if(!fs_ || !fs_->in_set(x)) func(x);};
            if(unsigned(k_) <= enc_.nremperres64() && !use128()) {
                enc_.for_each(lfunc, subpath.data());
            } else if(unsigned(k_) <= enc_.nremperres128()) {
                auto encoder(enc_.to_u128());
                encoder.for_each(lfunc, subpath.data());
            } else {
                use128() ? rh128_.for_each_hash(lfunc, subpath.data()): rh_.for_each_hash(lfunc, subpath.data());
            }
        }, path);
    };
    perf_for_substrs([fs=fs_.get()](auto x) {fs->add(x);});
}
void Dashing2Options::validate() const {
    if(canonicalize() && rh_.hashtype() != bns::DNA) {
        canonicalize(false);
    }
}

uint64_t XORMASK = 0x724526e320f9967dull;
u128_t XORMASK2 = (u128_t(12499408336417088522ull) << 64) | XORMASK;
void seed_mask(uint64_t x) {
    if(x == 0) {
        XORMASK = 0; XORMASK2 = 0;
    } else {
        wy::WyRand<uint64_t> rng(x);
        XORMASK = rng();
        XORMASK2 = (u128_t(rng()) << 64) | XORMASK;
    }
}

} // dashing2

int main_usage() {
    std::fprintf(stderr, "dashing2 has several subcommands\n");
    std::fprintf(stderr, "Usage can be seen in those commands.\n");
    std::fprintf(stderr, "sketch: converts FastX into k-mer sets/sketches, and sketches BigWig and BED files; also contains functionality from cmp, for one-step sketch and comparisons\n"
                         "This is probably the most common subcommand to use.\n"
    );
    std::fprintf(stderr, "cmp, a.k.a. dist: compres previously sketched/decomposed k-mer sets and emits results.\n");
    std::fprintf(stderr, "wsketch: Takes a tuple of [1-3] input binary files [(u32 or u64), (float or double), (u32 or u64)] and performs weighted minhash sketching.\n"
                         "sketch is for parsing and sketching (from Fast{qa}, BED, BigWig) and wsketch is for sketching binary files which have already been summed\n");
    return 1;
}
using namespace dashing2;


int main(int argc, char **argv) {
    if(argc > 1) {
        if(std::strcmp(argv[1], "sketch") == 0)
            return sketch_main(argc - 1, argv + 1);
        if(std::strcmp(argv[1], "cmp") == 0 || std::strcmp(argv[1], "dist") == 0)
            return cmp_main(argc - 1, argv + 1);
        if(std::strcmp(argv[1], "wsketch") == 0)
            return wsketch_main(argc - 1, argv + 1);
    }
    return main_usage();
}
