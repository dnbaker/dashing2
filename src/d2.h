#ifndef D2_H_H___
#define D2_H_H___
#ifdef _OPENMP
#include <omp.h>
#endif
#include "enums.h"
#include <memory>
#include <vector>
#include "bonsai/encoder.h"
#include "xxHash/xxh3.h"
#include "src/setsketch.h"
#include "sketch/bmh.h"
#include "counter.h"
#include "oph.h"
#include "filterset.h"


namespace dashing2 {

using namespace std::literals::string_literals;
using namespace sketch;

// To allow for 64-bit set identifiers, compile with -DLSHIDTYPE=uint64_t
#ifndef LSHIDTYPE
#define LSHIDTYPE uint32_t
#endif
using LSHIDType = LSHIDTYPE;
#undef LSHIDTYPE

#ifndef DASHING2_INDEX_FLOAT_TYPE
#define DASHING2_INDEX_FLOAT_TYPE float
#endif
using LSHDistType = DASHING2_INDEX_FLOAT_TYPE;
#undef DASHING2_INDEX_FLOAT_TYPE



struct IntervalSketchResult {
    using Map = flat_hash_map<std::string, std::vector<RegT>>;
    std::unique_ptr<Map> chrmap_;
    flat_hash_map<std::string, double> cardmap_;
    std::unique_ptr<std::vector<RegT>> global_;
    double card_;
};

template<typename F>
void for_each_substr(const F &func, const std::string &s, const int sep=' ') {
    const char *p;
    if((p = std::strchr(s.data(), sep)) == nullptr) {
        func(s);
        return;
    }
    const char *p2 = s.data();
    std::string tmp(p2, p);
    for(;;) {
        func(tmp);
        std::swap(p2, ++p);
        if((p = std::strchr(p2, sep)) == nullptr) {
            tmp = p2;
            func(tmp);
            break;
        }
        tmp = std::string(p2, p);
        if(std::all_of(tmp.begin(), tmp.end(), [](auto x) {return std::isspace(x);})) break;
    }
}

static inline bool check_compressed(std::string &path, int &ft) {
    if(bns::isfile(path)) {
        ft = 0;
        return true;
    }
    std::string opath = path + ".gz";
    if(bns::isfile(opath)) {
        ft = 1; path = opath;
        return true;
    }
    opath = path + ".xz";
    if(bns::isfile(opath)) {
        ft = 2; path = opath;
        return true;
    }
    return false;
}

struct Dashing2Options {

    // K-mer options
    int k_, w_;
    bns::Spacer sp_;
    mutable bns::Encoder<> enc_;
    mutable bns::RollingHasher<uint64_t> rh_;
    mutable bns::RollingHasher<u128_t> rh128_;
    bns::RollingHashingType rht_;
    bool parse_by_seq_ = false;
    bool trim_chr_ = true;
    size_t sketchsize_ = 2048;
    uint32_t count_threshold_ = 0;
    KmerSketchResultType kmer_result_;
    bool by_chrom_ = false;
    bool bed_parse_normalize_intervals_ = false;
    size_t cssize_ = 0;
    bool save_kmers_ = false;
    bool save_kmercounts_ = false;
    bool homopolymer_compress_minimizers_ = false;
private:
    bool trim_folder_paths_ = false; // When writing output files, write to cwd instead of the directory the files came from
public:
    bool cache_sketches_ = false;
    std::string outprefix_;
    std::string spacing_;
    double kmer_downsample_frac_ = 1.;
    uint64_t sampler_rng_;
    uint64_t sampler_threshold_;
    uint64_t seedseed_ = 0;
    double fd_level_ = sizeof(RegT); // Determines the number of bytes to which the minhash sketches are compressed
    int truncation_method_ = 0;

    // Whether to sketch multiset, set, or discrete probability distributions

    SketchSpace sspace_;
    DataType dtype_;
    bool use128_ = false;
    unsigned nthreads_;
    mutable long double compressed_b_ = -1.L, compressed_a_ = -1.L;
    bool fasta_dedup_ = false;

    std::unique_ptr<FilterSet> fs_;
    Dashing2Options(int k, int w=-1, bns::RollingHashingType rht=bns::DNA, SketchSpace space=SPACE_SET, DataType dtype=FASTX, size_t nt=0, bool use128=false, std::string spacing="", bool canon=false, KmerSketchResultType kres=ONE_PERM):
        k_(k), w_(w), sp_(k, w > 0 ? w: k, spacing.data()), enc_(sp_, canon), rh_(k, canon, rht, w), rh128_(k, canon, rht, w), rht_(rht), spacing_(spacing), sspace_(space), dtype_(dtype), use128_(use128) {
        kmer_result_ = kres;
#ifndef NDEBUG
        if(dtype_ == FASTX) {
            std::fprintf(stderr, "Dashing2 made with k = %d, w = %d, %s target, space = %s, datatype = %s and result = %s\n", k, w, rht == bns::DNA ? "DNA": "Protein", ::dashing2::to_string(sspace_).data(), ::dashing2::to_string(dtype_).data(), ::dashing2::to_string(kmer_result_).data());
        } else {
            std::fprintf(stderr, "Dashing2 made with space = %s, datatype = %s and result = %s\n", ::dashing2::to_string(sspace_).data(), ::dashing2::to_string(dtype_).data(), ::dashing2::to_string(kmer_result_).data());
        }
#endif
        if(nt <= 0) {
            DBG_ONLY(std::fprintf(stderr, "[%s:%s:%d] num threads < 0, checking OMP_NUM_THREADS\n", __FILE__, __func__, __LINE__);)
            if(char *s = std::getenv("OMP_NUM_THREADS"))
                nt = std::max(std::atoi(s), 1);
        }
        if(nt < 1) nt = 1;
        nthreads(nt);
        ht(rht);
    }
    void w(int neww) {w_ = neww; sp_.resize(k_, w_); rh128_.window(neww); rh_.window(neww);}
    std::string to_string() const;
    void validate() const;
#define D2O(name, oname)\
    Dashing2Options &oname(decltype(name) x) {name = x; return *this;}\
    std::add_const_t<decltype(name)> &oname() const {return name;}\
    decltype(name) &oname() {return name;}
#define D2O2(name) D2O(name##_, name)
    D2O2(outprefix) D2O2(save_kmers)
    D2O2(save_kmercounts) D2O2(homopolymer_compress_minimizers)
    D2O2(kmer_result) D2O2(use128) D2O2(cache_sketches)
    D2O2(sketchsize) D2O2(cssize) D2O2(parse_by_seq)
    D2O2(count_threshold)
    D2O2(fasta_dedup);
#undef D2O
#undef D2O2
    void downsample(double f) {
        if(f < 0. || f > 1.) throw std::runtime_error("Can't downsample to anything > 1 or < 0");
        kmer_downsample_frac_ = f;
        std::memcpy(&sampler_rng_, &f, 8);
        sampler_threshold_ = std::ceil(uint64_t(-1) * f);
    }
    INLINE bool downsample_pass() {
        return kmer_downsample_frac_ == 1. ||
               wy::wyhash64_stateless(&sampler_rng_) < sampler_threshold_;
    }
    // Getters and setters for all of the above
    Dashing2Options &parse_bigwig() {dtype_ = BIGWIG; return *this;}
    Dashing2Options &parse_bed() {dtype_ = BED; return *this;}
    bns::InputType input_mode() const {return rht_;}
    Dashing2Options &ht(bns::InputType rt) {
        rht_ = rt;
        rh_.hashtype(rt); rh128_.hashtype(rt);
        enc_.hashtype(rt);
        assert(rh_.hashtype() == rt);
        assert(enc_.hashtype() == rt);
        return *this;
    }
    bns::InputType hashtype() const {return rht_;}
    Dashing2Options &parse_protein() {return ht(bns::PROTEIN);}
    Dashing2Options &parse_protein3bit() {return ht(bns::PROTEIN_3BIT);}
    Dashing2Options &parse_protein6() {return ht(bns::PROTEIN_6);}
    Dashing2Options &parse_protein14() {return ht(bns::PROTEIN_14);}
    Dashing2Options &parse_dna2() {return ht(bns::DNA2);}
    Dashing2Options &parse_dnac() {return ht(bns::DNAC);}
    Dashing2Options &nthreads(int nt) {
        nt = std::max(nt, 1);
        nthreads_ = nt;
        OMP_ONLY(omp_set_num_threads(nt);)
        return *this;
    }
    void filterset(const std::string &path, bool is_kmer);
    void filterset(const std::string &fsarg);
    CountingType ct() const {return cssize_ > 0 ? COUNTMIN_COUNTING: EXACT_COUNTING;}
    CountingType count() const {return ct();}
    bool trim_folder_paths() const {
        return trim_folder_paths_ || outprefix_.size();
    }
    bool canonicalize() const {return enc_.canonicalize();}
    void canonicalize(bool val) const {
        enc_.canonicalize(val);
        rh_.canonicalize(val);
        rh128_.canonicalize(val);
    }
    auto w() const {return w_;}
    bool one_perm() const {return kmer_result_ == ONE_PERM && sspace_ == SPACE_SET;}
    auto nthreads() const {return nthreads_;}
    size_t nremperres64() const {return enc_.nremperres64();}
    size_t nremperres128() const {return enc_.nremperres128();}
    uint64_t seedseed() const {return seedseed_;}
    Dashing2Options &seedseed(uint64_t seed) {seedseed_ = seed; seed_mask(seedseed_); return *this;}
    bool sketch_compressed() const {
        return std::min(compressed_a_, compressed_b_) > 0.L;
        // Note: this always returns True after make_compressed (cmp_main.{h,cpp}) if --fastcmp is less than 8.
        // This is meant solely to be used during sketching
    }
    int sigshift() const {
        return (fd_level_ == 1. ? 3: fd_level_ == 2. ? 2: fd_level_ == 4. ? 1: fd_level_ == 0.5 ? 4: fd_level_ == 8 ? 0: -1) + (sizeof(RegT) == 16);
    }
};

static INLINE bool endswith(std::string lhs, std::string rhs) {
    return std::equal(rhs.begin(), rhs.end(), &lhs[lhs.size() - rhs.size()]);
}


using KmerSigT = std::conditional_t<(sizeof(RegT) == 8), uint64_t, std::conditional_t<(sizeof(RegT) == 4), uint32_t, u128_t>>;
using FullSetSketch = sketch::setsketch::CountFilteredCSetSketch<RegT>;
using OPSetSketch = LazyOnePermSetSketch<KmerSigT>;
using BagMinHash = sketch::BagMinHash2<RegT>;
using ProbMinHash = sketch::pmh2_t<RegT>;
using OrderMinHash = sketch::omh::OMHasher<RegT>;
template<typename T>
INLINE auto total_weight(const T &x) {return x.total_weight();}
template<>
INLINE auto total_weight(const FullSetSketch &x) {return x.total_updates();}
template<>
INLINE auto total_weight(const OPSetSketch &x) {return x.total_updates();}

static constexpr size_t nregperitem(bns::RollingHashingType it, bool is128=false) {
    using namespace bns;
    switch(it) {
        case DNA: return is128 ? RHTraits<DNA>::nper128: RHTraits<DNA>::nper64;
        case PROTEIN: return is128 ? RHTraits<PROTEIN>::nper128: RHTraits<PROTEIN>::nper64;
        case PROTEIN20: return is128 ? RHTraits<PROTEIN20>::nper128: RHTraits<PROTEIN20>::nper64;
        case PROTEIN_3BIT: return is128 ? RHTraits<PROTEIN_3BIT>::nper128: RHTraits<PROTEIN_3BIT>::nper64;
        case PROTEIN_14: return is128 ? RHTraits<PROTEIN_14>::nper128: RHTraits<PROTEIN_14>::nper64;
        case PROTEIN_6: return is128 ? RHTraits<PROTEIN_6>::nper128: RHTraits<PROTEIN_6>::nper64;
        case DNAC: return is128 ? RHTraits<DNAC>::nper128: RHTraits<DNAC>::nper64;
        case PROTEIN_6_FRAME: return is128 ? RHTraits<PROTEIN_6_FRAME>::nper128: RHTraits<PROTEIN_6_FRAME>::nper64;
        default: ;
    }
    return 0; // Should not ever happen
}

struct KSeqHolder {
    kseq_t *kseqs_;
    size_t n_;
    KSeqHolder(size_t n): kseqs_(static_cast<kseq_t *>(std::calloc(n, sizeof(kseq_t)))), n_(n) {
        if(!kseqs_) THROW_EXCEPTION(std::bad_alloc());
        for(auto p = kseqs_; p < kseqs_ + n_; ++p)
            ks_resize(&p->seq, 1<<18);
    }
    KSeqHolder(const KSeqHolder &o) = delete;
    KSeqHolder(KSeqHolder &&o): kseqs_(o.kseqs_), n_(o.n_) {
        o.n_ = 0;
        o.kseqs_ = 0;
    }
    kseq_t &operator[](size_t i) {return kseqs_[i];}
    const kseq_t &operator[](size_t i) const {return kseqs_[i];}
private:
    void free_item(kseq_t &seq) {
        std::free(seq.name.s);
        std::free(seq.comment.s);
        std::free(seq.seq.s);
        std::free(seq.qual.s);
        ks_destroy(seq.f);
    }
public:
    ~KSeqHolder() {
        for(size_t i = 0; i < n_; free_item(kseqs_[i++]));
        std::free(kseqs_);
        kseqs_ = 0;
        n_ = 0;
    }
};


extern bool entmin;

} // namespace dashing2
//std::vector<RegT> reduce(flat_hash_map<std::string, std::vector<RegT>> &map);

#endif
