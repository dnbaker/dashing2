#include "bedsketch.h"

namespace dashing2 {

std::pair<std::vector<RegT>, double> bed2sketch(const std::string &path, const Dashing2Options &opts) {
    if(opts.sspace_ > SPACE_PSET) throw std::invalid_argument("Can't do edit distance for BED files");
    if(opts.bed_parse_normalize_intervals_ && opts.sspace_ == SPACE_SET)
        throw std::invalid_argument("Can't normalize BED rows in set space. Use SPACE_MULTISET or SPACE_PSET");
    std::ifstream ifs(path);
    const bool op = opts.one_perm();
    FullSetSketch ss(opts.count_threshold_, opts.sketchsize_);
    OPSetSketch opss(opts.sketchsize_);
    Counter ctr(opts.cssize_);
    std::pair<std::vector<RegT>, double> ret({std::vector<RegT>(opts.sketchsize_), 0.});
    auto &retvec(ret.first);
    std::string cache_path = path + to_suffix(opts);
    DBG_ONLY(std::fprintf(stderr, "Using %s\n", op ? "oneperm": "fullsetsketch");)

    if(opts.trim_folder_paths_) {
        cache_path = trim_folder(path);
        if(opts.outprefix_.size())
            cache_path = opts.outprefix_ + '/' + cache_path;
    }
    if(opts.cache_sketches_ && bns::isfile(cache_path)) {
        auto nb = bns::filesize(cache_path.data());
        retvec.resize(nb / sizeof(RegT));
        std::FILE *ifp = std::fopen(cache_path.data(), "rb");
        if(!ifp) throw 1;
        auto rc = std::fread(retvec.data(), nb, 1, ifp);
        ret.second = retvec.size() / std::accumulate(retvec.begin(), retvec.end(), 0.L);
        std::fclose(ifp);
        if(rc == 1u) return ret;
        else
            std::fprintf(stderr, "Failed to read from disk; instead, sketching from scratch (%s)\n", path.data());
    }
    for(std::string s;std::getline(ifs, s);) {
        if(s.empty() || s.front() == '#') continue;
        char *p = s.data(), *p2;
        if((p2 = std::strchr(p, '\t')) == nullptr)
            throw std::invalid_argument(std::string("Malformed line: ") + s);
        if(opts.trim_chr_ && ((*p == 'c' || *p == 'C') && p[1] == 'h' && p[2] == 'r'))
            p += 3;
        const uint64_t chrhash = XXH3_64bits(p, p2 - p);
        p = p2 + 1;
        const unsigned long start = std::strtoul(p, &p, 10), stop =  std::strtoul(p, &p, 10);
        const double inc = opts.bed_parse_normalize_intervals_ ? 1. / (stop - start): 1;
        // Consider SIMDifying packing these before adding?
        // If set space, sketch directly.
        if(opts.sspace_ == SPACE_SET) {
            if(op) for(auto i = start; i < stop; opss.update(chrhash ^ i++));
            else   for(auto i = start; i < stop; ss.update(chrhash ^ i++));
        } else {
            // else, we need to compute counts before we sketch
            for(auto i = start; i < stop; ctr.add(chrhash ^ i++, inc));
        }
    }
    std::FILE *ofp = std::fopen(cache_path.data(), "w");
    if(opts.sspace_ > SPACE_SET) {
        if(opts.ct() == EXACT_COUNTING) {
            if(opts.sspace_ == SPACE_MULTISET) {
                BagMinHash bmh(opts.sketchsize_);
                ctr.finalize(bmh);
                std::copy(bmh.data(), bmh.data() + opts.sketchsize_, retvec.data());
                std::fwrite(bmh.data(), opts.sketchsize_, sizeof(RegT), ofp);
                ret.second = bmh.total_weight();
            } else {
                sketch::pmh2_t pmh(opts.sketchsize_);
                ctr.finalize(pmh);
                std::copy(pmh.data(), pmh.data() + opts.sketchsize_, retvec.data());
                std::fwrite(pmh.data(), opts.sketchsize_, sizeof(RegT), ofp);
                ret.second = pmh.total_weight();
            }
        } else {
#define __FS() do {\
    for(size_t i = 0; i < csz; ++i) sketcher.update(i, ctr.count_sketch_[i]);\
    auto p = sketcher.data();\
    std::fwrite(p, opts.sketchsize_, sizeof(RegT), ofp);\
    std::copy(p, p + opts.sketchsize_, retvec.data());\
    ret.second = sketcher.total_weight();\
    } while(0)
            const size_t csz = ctr.count_sketch_.size();
            if(opts.sspace_ == SPACE_MULTISET) {
                BagMinHash sketcher(opts.sketchsize_);
                __FS();
            } else {
                sketch::pmh2_t sketcher(opts.sketchsize_);
                __FS();
            }
#undef __FS
        }
    } else {
        ret.second = op ? opss.getcard(): ss.getcard();
        RegT *sptr = op ? opss.data(): ss.data();
        std::copy(sptr, sptr + opts.sketchsize_, retvec.data());
        std::fwrite(sptr, opts.sketchsize_, sizeof(RegT), ofp);
    }
    std::fclose(ofp);
    return ret;
}

} // namespace dashing2
