#include <memory>
#include <vector>
#include "bonsai/encoder.h"
#include "bigWig.h"
#include "flat_hash_map/flat_hash_map.hpp"


using u128_t = __uint128_t;
static inline uint64_t encode_refpos(int tid, int pos) {return uint64_t(tid) << 32 | pos;}
static inline u128_t encode_refpos(int tid, uint64_t pos) {return u128_t(tid) << 64 | pos;}


enum DataType {
    FASTX,
    BIGWIG,
    BED
};

enum SketchSpace {
    SPACE_SET,      // MinHash/SetSketch/HLL
    SPACE_MULTISET, // Weighted MinHash -- e.g., Bag or Tree
    SPACE_PSET,      // ProbMinHash
    SPACE_EDIT_DISTANCE // edit distance -- implies OMH and forces us to use strings as input instead of bags of k-mers
};

enum CountingType {
    EXACT_COUNTING,
    COUNTSKETCH_COUNTING
};

using u128_t = __uint128_t;
namespace hash = sketch::hash;

struct FHasher {
    hash::FusedReversible3<hash::InvMul, hash::RotL33, hash::MultiplyAddXoRot<16>> rhasher_;
    FHasher() {}
    INLINE uint64_t operator()(__uint128_t x) const {
        return rhasher_(uint64_t(x>>64)) ^ rhasher_(uint64_t(x));
    }
};

struct Counter {
    CountingType ct_;
    ska::flat_hash_map<uint64_t, uint32_t> c64_;
    ska::flat_hash_map<u128_t, uint32_t, FHasher> c128_;
    std::vector<int32_t> count_sketch_;
    schism::Schismatic<uint64_t> s64_;
    Counter(CountingType t, size_t cssize=50'000'000): ct_(t), count_sketch_(t == EXACT_COUNTING? size_t(0): cssize), s64_(cssize ? cssize: size_t(1)) {
    }
    static constexpr uint64_t BM64 = 0x8000000000000000ull;
    void add(u128_t x) {
        if(ct_ == EXACT_COUNTING) {
            auto it = c128_.find(x);
            if(it == c128_.end()) c128_.emplace(x, 1);
            else ++it->second;
        } else {
            const auto hv = sketch::hash::WangHash::hash(uint64_t(x)) ^ sketch::hash::WangHash::hash(x >> 64);
            count_sketch_[s64_.mod(hv)] += ((hv & BM64) ? 1: -1);
        }
    }
    void add(uint64_t x) {
        if(ct_ == EXACT_COUNTING) {
            auto it = c64_.find(x);
            if(it == c64_.end()) c64_.emplace(x, 1);
            else ++it->second;
        } else {
            const auto hv = sketch::hash::WangHash::hash(x);
            count_sketch_[s64_.mod(hv)] += ((hv & BM64) ? 1: -1);
        }
    }
};

struct ParseOptions {

    // K-mer options
    int k_;
    bns::Spacer sp_;
    bns::Encoder<> enc_;
    bns::RollingHasher<uint64_t> rh_;
    bns::RollingHasher<u128_t> rh128_;
    bns::RollingHashingType rht_;
    bool parse_by_seq_ = false;

    // Whether to sketch multiset, set, or discrete probability distributions

    SketchSpace sspace_;
    CountingType count_;
    DataType dtype_;
    ParseOptions(int k, int w=-1, std::string spacing="", bns::RollingHashingType rht=bns::DNA, SketchSpace space=SPACE_SET, CountingType count=EXACT_COUNTING, DataType dtype=FASTX):
        k_(k), sp_(k, w > 0 ? w: k, spacing.data()), enc_(sp_), rh_(k), rh128_(k), rht_(rht), sspace_(space), count_(count), dtype_(dtype) {}
    ParseOptions &parse_by_seq() {parse_by_seq_ = true; return *this;}
    ParseOptions &parse_bigwig() {dtype_ = BIGWIG; return *this;}
    ParseOptions &parse_bed() {dtype_ = BED; return *this;}
    ParseOptions &parse_protein() {rh_.enctype_ = rh128_.enctype_ = rht_ = bns::PROTEIN; return *this;}
};

struct OutputKind {
    SYMMETRIC_ALL_PAIRS,
    ASYMMETRIC_ALL_PAIRS,
    KNN_GRAPH,
    NN_GRAPH_THRESHOLD // Variable number of similarities, as given by threshold
};

struct Collection {
    ParseOptions opts_;
    std::vector<std::string> names;
    std::vector<std::vector<uint64_t>> sketch_data;
    std::vector<std::vector<uint64_t>> full_kmers; // Optional
    // For BigWig and BED files, full_kmers will be empty,
    std::vector<std::string> full_strings; // For edit-distance clustering only.
    OutputKind kind_;
    double threshold_; // if topk_ < 0, use -fraction to keep only items above *fraction* similarity
                       // if topk_ <= -1, then emit all comparisons
};

int main() {
}
