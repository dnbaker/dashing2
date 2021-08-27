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


std::string to_string(OutputFormat of) {
    if(of == HUMAN_READABLE) return "HumanReadable";
    return "MachineReadable";
}
std::string to_string(OutputKind ok) {
    if(ok == SYMMETRIC_ALL_PAIRS) return "UpperTriangularSymmetricAllPairs";
    if(ok == ASYMMETRIC_ALL_PAIRS) return "AllPairs";
    if(ok == KNN_GRAPH) return "KNNGraph";
    if(ok == NN_GRAPH_THRESHOLD) return "ThresholdedNNGraph";
    return "Deduplication (not supported yet)";
}

void checked_fwrite(std::FILE *const fp, const void *const ptr, const size_t nb) {
    unsigned long long lrc = std::fwrite(static_cast<const void *>(ptr), 1, nb, fp);
    if(unlikely(lrc != static_cast<size_t>(nb)))
         throw std::runtime_error(std::string("[E:") + __PRETTY_FUNCTION__ + ':' + __FILE__ + std::to_string(__LINE__) + "] Failed to perform buffered write of " + std::to_string(static_cast<size_t>(nb)) + " bytes, instead writing " + std::to_string(lrc) + " bytes");
}

std::FILE *xopen(const std::string &path) {
    std::FILE *fp;
    if(path.size() > 3 && std::equal(path.data() + path.size() - 3, &path[path.size()], ".xz")) {
        auto cmd = std::string("xz -dc ") + path;
        fp = ::popen(cmd.data(), "r");
    } else if(path.size() > 3 && std::equal(path.data() + path.size() - 3, &path[path.size()], ".gz")) {
        auto cmd = std::string("gzip -dc ") + path;
        fp = ::popen(cmd.data(), "r");
    } else if(path.size() > 4 && std::equal(path.data() + path.size() - 3, &path[path.size()], ".bz2")) {
        auto cmd = std::string("bzip2 -dc ") + path;
        fp = ::popen(cmd.data(), "r");
    } else {
        fp = ::popen((std::string("cat ") + path).data(), "r");
    }
    return fp;
}
long signed int BLKSIZE = -1;
std::mutex blksizelock;

void buffer_to_blksize(std::FILE *fp) {
    if(BLKSIZE < 0) {
        struct stat st;
        ::fstat(STDIN_FILENO, &st);
        BLKSIZE = st.st_blksize;
    }
    std::setvbuf(fp, nullptr, _IOFBF, BLKSIZE);
}

std::FILE *bfopen(const char *s, const char *fmt) {
    std::FILE *ifp = std::fopen(s, fmt);
    if(ifp) buffer_to_blksize(ifp);
    return ifp;
}

uint64_t XORMASK = 0x724526e320f9967dull;
u128_t XORMASK2 = (u128_t(12499408336417088522ull) << 64) | XORMASK;
void seed_mask(uint64_t x) {
    if(x == 0) {
        XORMASK = 0; XORMASK2 = 0;
    } else {
        XORMASK = sketch::hash::WangHash::hash(x);
        XORMASK2 = (XORMASK | (u128_t(sketch::hash::WangHash::hash(XORMASK)) << 64));
    }
}

}
