#include "d2.h"
#include "enums.h"
#include <array>

namespace dashing2 {

std::string to_string(KmerSketchResultType t) {
    if(t == ONE_PERM) {return "OnePermutationSetSketch";}
    if(t == FULL_SETSKETCH) return "FullSetSketch";
    if(t == FULL_MMER_SET) return "FullMmerSet";
    if(t == FULL_MMER_SEQUENCE) return "FullMmerSequence";
    return "FullMmerCountdict";
}
std::string to_string(DataType dt) {
    if(dt == FASTX) return "Fastx";
    if(dt == BED) return "BED";
    if(dt == BIGWIG) return "BigWig";
    if(dt == LEAFCUTTER) return "LeafCutter";
    return "Unknown";
}

std::string trim_folder(const std::string &s) {
    auto pos = s.find_last_of("/");
    if(pos == std::string::npos) return s;
    return s.substr(pos + 1, std::string::npos);
}

std::string to_suffix(const Dashing2Options &opts) {
    std::string ret = (opts.kmer_result_ == ONE_PERM || opts.kmer_result_ == FULL_SETSKETCH) ?
        (opts.sspace_ == SPACE_SET ? (opts.kmer_result_ == ONE_PERM ? ".opss": ".ss"):  opts.sspace_ == SPACE_MULTISET ? ".bmh": opts.sspace_ == SPACE_PSET ? ".pmh": ".unknown")
        : (opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_COUNTDICT) ? ".kmerset" : opts.kmer_result_ == FULL_MMER_SEQUENCE ? ".mmerseq": ".unknown_kmer";
    if(opts.kmer_result_ == FULL_MMER_SEQUENCE || opts.kmer_result_ == FULL_MMER_SET || opts.kmer_result_ == FULL_MMER_COUNTDICT) {
        static constexpr std::array<const char *, 2> lut{"64", "128"};
        ret += lut[opts.use128()];
    }
    return ret;
}

std::string to_string(SketchSpace ss) {
    if(ss == SPACE_SET) return "SetSpace";
    if(ss == SPACE_MULTISET) return "MultisetSpace";
    if(ss == SPACE_PSET) return "ProbsetSpace";
    if(ss == SPACE_EDIT_DISTANCE) return "EditDistanceSpace";
    return "UNKNOWN SPACE";
}
std::string to_string(CountingType ct) {
    if(ct == EXACT_COUNTING) return "ExactCounting";
    if(ct == COUNTMIN_COUNTING) return "CountMinCounting";
    if(ct == COUNTSKETCH_COUNTING)
        return "CountSketchCounting";
    if(ct == CQF_COUNTING) return "CQFCounting";
    return "UNKNOWN COUNTING";
}

bool iscomp(const std::string &s) {
    bool ret = false;
    if(s.size() >= 3) {
        if(std::equal(&s[s.size() - 3], &s[s.size()], ".gz")) ret = true;
        else if(std::equal(&s[s.size() - 3], &s[s.size()], ".xz")) ret = true;
        else if(std::equal(&s[s.size() - 4], &s[s.size()], ".bz2")) ret = true;
    }
    std::fprintf(stderr, "%s is %d/%s\n", s.data(), int(ret), ret ? "compressed": "uncompressed");
    return ret;
}

}
