#include "bedsketch.h"

namespace dashing2 {

std::vector<RegT> bed2sketch(std::string path, const Dashing2Options &opts) {
    if(opts.sspace_ > SPACE_PSET) throw std::invalid_argument("Can't do edit distance for BED files");
    if(opts.bed_parse_normalize_intervals_ && opts.sspace_ == SPACE_SET)
        throw std::invalid_argument("Can't normalize BED rows in set space. Use SPACE_MULTISET or SPACE_PSET");
    std::ifstream ifs(path);
    FullSetSketch ss(opts.one_perm_ ? size_t(1): opts.sketchsize_);
    OPSetSketch opss(opts.one_perm_ ? opts.sketchsize_: size_t(1));
    Counter ctr(opts.cssize_);
    std::vector<RegT> ret(opts.sketchsize_);
    for(std::string s;std::getline(ifs, s);) {
        if(s.empty() || s.front() == '#') continue;
        char *p = &s[0], *p2 = std::strchr(p, '\t');
        if(!p2) goto err;
        if(opts.trim_chr_ && (std::memcmp(p, "chr", 3) == 0 || std::memcmp(p, "Chr", 3) == 0))
            p += 3;
        const uint64_t chrhash = XXH3_64bits(p, p2 - p);
        p = p2 + 1;
        unsigned long start = std::strtoul(p, &p, 10), stop =  std::strtoul(p, &p, 10);
        const float inc = opts.bed_parse_normalize_intervals_ ? 1. / (stop - start): 1.f;
        // If set space, sketch directly.
        if(opts.sspace_ == SPACE_SET) {
            if(opts.one_perm_) for(auto i = start; i < stop; ss.update(chrhash ^ i++));
            else               for(auto i = start; i < stop; opss.update(chrhash ^ i++));
        } else {
            // else, we need to compute counts before we sketch
            for(auto i = start; i < stop; ctr.add(chrhash ^ i++, inc));
        }
    }
    if(0) {
        err:
        throw std::invalid_argument("Malformed line");
    }
    if(opts.sspace_ > SPACE_SET) {
        if(opts.ct() == EXACT_COUNTING) {
            if(opts.sspace_ == SPACE_MULTISET) {
                BagMinHash bmh(opts.sketchsize_);
                for(const auto &pair: ctr.c64_) {
                    bmh.update(pair.first, pair.second);
                }
                std::copy(bmh.data(), bmh.data() + opts.sketchsize_, ret.data());
            } else {
                sketch::pmh2_t pmh(opts.sketchsize_);
                for(const auto &pair: ctr.c64_) {
                    pmh.update(pair.first, pair.second);
                }
                std::copy(pmh.data(), pmh.data() + opts.sketchsize_, ret.data());
            }
        } else {
            const size_t csz = ctr.count_sketch_.size();
            if(opts.sspace_ == SPACE_MULTISET) {
                BagMinHash bmh(opts.sketchsize_);
                for(size_t i = 0; i < csz; ++i) {
                    bmh.update(i, ctr.count_sketch_[i]);
                }
                std::copy(bmh.data(), bmh.data() + opts.sketchsize_, ret.data());
            } else {
                sketch::pmh2_t pmh(opts.sketchsize_);
                for(size_t i = 0; i < csz; ++i) {
                    pmh.update(i, ctr.count_sketch_[i]);
                }
                std::copy(pmh.data(), pmh.data() + opts.sketchsize_, ret.data());
            }
        }
    } else {
        RegT *sptr = opts.one_perm_ ? opss.data(): ss.data();
        std::copy(sptr, sptr + opts.sketchsize_, ret.data());
    }
    return ret;
}

} // namespace dashing2
