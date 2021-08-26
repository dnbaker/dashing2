#include "wsketch.h"
#include "enums.h"
#include "d2.h"
#include "sketch/bmh.h"
#include "mio.hpp"
namespace dashing2 {

template<typename BMH>
BMH construct(size_t sketchsize) {
    if constexpr(std::is_same_v<BMH, BagMinHash> || std::is_same_v<BMH, ProbMinHash>)
        return BMH(sketchsize, true);
    if constexpr(std::is_same_v<BMH, OPSetSketch>)
        return OPSetSketch(sketchsize);
    if constexpr(std::is_same_v<BMH, FullSetSketch>) return BMH(1, sketchsize, true); // FullSetSketch
    THROW_EXCEPTION(std::runtime_error(std::string("Failed to construct: this is the wrong type: ") + __PRETTY_FUNCTION__));
}
template<typename BMH, typename FT, typename IT, typename IndPtrT=uint64_t>
std::vector<SimpleMHRet> minhash_rowwise_csr(const FT *weights, const IT *indices, const IndPtrT *indptr, size_t nr, size_t m) {
    std::vector<SimpleMHRet> ret(nr);
    std::atomic<size_t> total_processed;
    total_processed.store(0);
    auto t1 = std::chrono::high_resolution_clock::now();
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for(size_t i = 0; i < nr; ++i) {
        DBG_ONLY(std::fprintf(stderr, "%zu/%zu\n", i, nr);)
        ++total_processed;
        BMH h(construct<BMH>(m));
        const size_t b = indptr[i], e = indptr[i + 1];
        for(size_t j = b; j < e; ++j) {
            h.update(j - b, weights ? weights[j]: FT(1));
        }
        if constexpr(!std::is_same_v<BMH, FullSetSketch>) h.finalize();
        std::get<3>(ret[i]) = total_weight(h);
        std::vector<uint64_t> ids(m);
        auto &hids = h.ids();
        std::transform(hids.begin(), hids.end(), ids.begin(), [ind=indices + b](auto x) {return ind[x];});
        std::get<0>(ret[i]) = h.template to_sigs<RegT>();
        std::get<1>(ret[i]) = h.template to_sigs<uint64_t>();
        std::get<2>(ret[i]) = ids;
        if((total_processed.load() & 1023u) == 0u) {
            std::chrono::duration<double, std::milli> diff(std::chrono::high_resolution_clock::now() - t1);
            auto left = nr - total_processed.load(), done = nr - left;
            auto sperit = diff.count() / (total_processed.load() + 1);
            auto timeleft = left * sperit;
            std::fprintf(stderr, "%g%% of total done, expected %gms time left. Processed %zu/%zu in %gms; Expected < %zu seconds left...\r\n", 100. * done / nr, timeleft, total_processed.load(), nr, diff.count(), size_t(std::ceil(timeleft / 1e3)));
        }
    }
    return ret;
}


template<typename BMH, typename FT, typename IT>
SimpleMHRet minwise_det(const FT *weights, const IT *indices, size_t n, size_t m) {
    BMH h(construct<BMH>(m));
    for(size_t i = 0; i < n; ++i) {
        // We use the offset as a first ID, and then we convert these to the
         // relevant IDs sampled at the end
        h.update(i, weights ? weights[i]: FT(1));
    }
    std::vector<uint64_t> ids(m);
    std::vector<RegT> regs;
    if constexpr(!std::is_same_v<BMH, FullSetSketch>) h.finalize();
    auto &hids = h.ids();
    std::transform(hids.begin(), hids.end(), ids.begin(), [ind=indices](auto x) {return ind[x];});
    SimpleMHRet ret;
    ret.sigs() = h.template to_sigs<RegT>();
    ret.hashes() = h.template to_sigs<uint64_t>();
    ret.ids() = ids;
    ret.total_weight() = total_weight(h);
    return ret;
}
template<typename FT, typename IT>
SimpleMHRet minhash(const FT *weights, const IT *indices, size_t n, size_t m, int usepmh) {
    std::fprintf(stderr, "Made mw det, w = %p, i = %p, size = %zu, %d pmh\n", (void *)weights, (void *)indices, sizeof(IT), usepmh);
    if(usepmh == 0)
        return minwise_det<ProbMinHash>(weights, indices, n, m);
    if(usepmh == 1) return minwise_det<BagMinHash>(weights, indices, n, m);
    return minwise_det<FullSetSketch>(weights, indices, n, m);
}

#define SIG(T1, T2, P1, P2) \
SimpleMHRet minhash##T1##T2(const P1 *weights, const P2 *indices, size_t n, size_t m, int usepmh) { \
    return minhash(weights, indices, n, m, usepmh);}\

#undef SIG

struct FReader{
    std::string path_;
    std::FILE *fp_;
    bool pclose = true;
    FReader(std::string path): path_(path), fp_(nullptr) {
        std::fprintf(stderr, "Reading from %s\n", path.data());
        std::string c = "cat ";
        if(endswith(path, ".gz"))
            c = "gzip -dc ";
        else if(endswith(path, ".xz"))
            c = "xz -dc ";
        else if(endswith(path, ".bz2"))
            c = "bzip2 -dc ";
        else {
            pclose = false;
            fp_ = bfopen(path_.data(), "r");
            return;
        }
        c += path;
        std::fprintf(stderr, "Using cmd = '%s'\n", c.data());
        if((fp_ = ::popen(c.data(), "r")) == nullptr)
            THROW_EXCEPTION(std::runtime_error(std::string("Failed to run cmd '") + c + "'"));
    }
    template<typename VT>
    std::vector<VT> getvec() {
        //std::fprintf(stderr, "%s Getting type of size %zu from %s\n", __PRETTY_FUNCTION__, sizeof(VT), path_.data());
        std::vector<VT> ret;
        for(VT v; std::fread(&v, sizeof(v), 1, fp_) == 1u;ret.push_back(v));
        //std::fprintf(stderr, "ret size is %zu\n", ret.size());
        return ret;
    }
    static inline bool endswith(std::string lhs, std::string rhs) {
        return std::equal(rhs.begin(), rhs.end(), &lhs[lhs.size() - rhs.size()]);
    }
    void clear() {
        if(fp_) {if(pclose) ::pclose(fp_); else std::fclose(fp_); fp_ = nullptr;}
    }
    ~FReader() {
        if(fp_) {clear();}
    }
};

template<typename T>
std::vector<T> fromfile(std::string path) {
    std::FILE *fp = bfopen(&path[0], "rb");
    std::vector<T> ret;
    for(T v;std::fread(&v, sizeof(T), 1, fp) == 1;ret.push_back(v));
    std::fclose(fp);
    return ret;
}

std::vector<SimpleMHRet> wmh_from_file_csr(std::string idpath, std::string cpath, std::string indptrpath, size_t sksz, int usepmh, int usef32=0, bool wordids=false, bool ip32=false) {
    // default is f64
    // For IDs, default is 64 bits
    //
    std::vector<SimpleMHRet> ret;
#define PERF1(T) minhash_rowwise_csr<T>(lp, rp, ip, nr, sksz)
#define PERF(TL, TR, IR) do {\
            std::unique_ptr<mio::mmap_sink> lf;\
            if(cpath.size() && cpath != "-") lf.reset(new mio::mmap_sink(cpath));\
            mio::mmap_sink rf(idpath);\
            const std::vector<IR> cvec(fromfile<IR>(indptrpath));\
            const TL *lp = lf ? (const TL *)lf->data(): (const TL *)nullptr;\
            const TR *rp = (const TR *)rf.data();\
            const IR *ip = cvec.data();\
            const size_t nr = cvec.size() - 1;\
            assert(!lp || lf->size() / sizeof(TL) == rf.size() / sizeof(TR));\
            ret = usepmh == 1 ? PERF1(ProbMinHash)\
                              : usepmh == 0 ? PERF1(BagMinHash)\
                                            : PERF1(FullSetSketch);\
        } while(0)
#define PERF2(LT) do {\
        if(wordids) {\
            if(ip32) {\
                PERF(LT, uint32_t, uint32_t);\
            } else {\
                PERF(LT, uint32_t, uint64_t);\
            }\
        } else if(ip32) {\
            PERF(LT, uint64_t, uint32_t);\
        } else {\
            PERF(LT, uint64_t, uint64_t);\
        }\
    } while(0)
    if(usef32 == 1) {
        PERF2(float);
    } else if(usef32 == -1) {
        PERF2(uint16_t);
    } else if(usef32 == -2) {
        PERF2(uint32_t);
    } else {
        PERF2(double);
    }
    return ret;
#undef PERF
#undef PERF1
#undef PERF2
}
SimpleMHRet wmh_from_file(std::string idpath, std::string cpath, size_t sksz, int usepmh, int usef32, bool wordids) {

    // usef32 means use float32 instead of float64 for weights
    // default is f64
    // For IDs, default is 64 bits
    //
    FReader rhs(idpath);
#define PERF(TL, TR) do {std::vector<TL> lvec; if(cpath.size()) {FReader lhs(cpath); lvec = lhs.getvec<TL>();} auto rvec = rhs.getvec<TR>();\
            assert(lvec.size() == rvec.size());\
            const TL *ptr = lvec.size() ? lvec.data(): static_cast<const TL *>(nullptr);\
            return minhash(ptr, rvec.data(), rvec.size(), sksz, usepmh);\
        } while(0)
#define PERF2(LT) do {\
        if(wordids) {\
            PERF(LT, uint32_t);\
        } else {\
            PERF(LT, uint64_t);\
        } \
    } while(0)
    if(usef32 == 1) {
        PERF2(float);
    } else if(usef32 == -1) {
        PERF2(uint16_t);
    } else {
        PERF2(double);
    }
    THROW_EXCEPTION(std::runtime_error("This should never happen"));
#undef PERF
#undef PERF2
}

int wsketchusage() {
    std::fprintf(stderr, "Sketch raw IDs, with optional weights added\n"
                         "Usage: dashing2 wsketch [input.bin] <Optional: input.weights.bin> <Optional: indptr.bin for CSR data>\n"
                         "If only one path is provided, it treated as indices, and sketched via SetSketch; IE, everything is sketched with equal weight.\n"
                         "If two paths are provided, the second is treated as a weight vector, and the multiset is sketched via ProbMinHash or BagMinHash.\n"
                         "If three paths are provided, the second is treated as a weight vector, and the last is used as indptr; this yields a stacked set of sketches corresponding to the input matrix.\n"
                         "-S: set sketch size\n"
                         "Identifier size: 64-bit by default.\n"
                         "-u: Read 32-bit identifiers from input.bin rather than 64-bit\n"
                         "If your data has non-integral types, you can still convert them to b-bit signatures using hash functions/RNGs\n"
                         "Weight data types:\n"
                         "If a weights file is provided, double is the expected format.\n"\
                         "-f: Read 32-bit floating point weights from [input.weights.bin] \n"
                         "-H: Read 16-bit data weights from [input.weight.bin] (Default: float64)\n"
                         "-U: Read 32-bit data weights from [input.weight.bin] (Default: float64)\n"
                         "Indptr data types:\n"
                         "If any indptr file is provided, it is expected to use 64-bit integers.\n"
                         "-P: Read 32-bit indptr integers [indptr.bin] (Default: uint64_t)\n"
                         "Sketching options: \n"
                         "By default, sketches with ProbMinHash, which treats datasets as discrete probability distributions\n"
                         "This can be changed with the following flags:\n\n"
                         "-B: Sketch with BagMinHash. This treats datasets as multisets.\n"
                         "-q: Sketch with SetSketch\n"
                         "-o: outprefix. If unset, uses [input.bin]\n"
                         "Runtime options:\n"
                         "-p: Set number of threads (processes) [1]\n"
            );
        std::fprintf(stderr, "If two are provided, then 1-D weighted minhashing is performed on the compressed vector. If three are passed, then this result is treated as a CSR-format matrix and then emits a matrix of sketches.\n");
        std::fprintf(stderr, "To unweighted sparse matrices (omitting the data field), use CSR-style sketching but replace the weights file with '-'");
        std::fprintf(stderr, "Example: 'dashing2 wsketch -S 64  -o g1.k31.k64 g1.fastq.k31.kmerset64 g1.fastq.k31.kmercounts.f64'.\n");
        std::fprintf(stderr, "Example: 'dashing2 wsketch -S 64  -o g1.k31.k64 g1.fastq.k31.kmerset64 # sketches k-mers only'.\n");
        std::fprintf(stderr, "Example: 'dashing2 wsketch -S 64  -o g1.k31.k64.mat g1.fastq.k31.data64 fq.fastq.k31.indices64 # sketches weighted k-mer sets'.\n");
        std::fprintf(stderr, "Example: 'dashing2 wsketch -S 64  -o g1.k31.k64.mat g1.fastq.k31.data64 fq.fastq.k31.indices64 fq.fastq.k31.indptr64 # sketches weighted sets from CSR-format and emits these sketches stacked'.\n");
        std::fprintf(stderr, "Example: 'dashing2 wsketch -S 64  -o g1.k31.k64.mat - fq.fastq.k31.indices64 fq.fastq.k31.indptr64 # sketches sets from packed sets and emits these sketches stacked.");
        std::fprintf(stderr, "The use of '-' causes the weights for all points to be uniform.\n");
    return 1;
}
template<typename T>
void write_container(const T &vec, std::string path) {
    const size_t nb = vec.size();
    std::FILE *ifp = bfopen(path.data(), "wb");
    if(!ifp)
        THROW_EXCEPTION(std::runtime_error(std::string("Failed to open path '") + path + "' for reading"));
    static constexpr size_t vnb = sizeof(std::decay_t<decltype(*vec.begin())>);
    std::fwrite(vec.data(), vnb, nb, ifp);
    std::fclose(ifp);
}
int wsketch_main(int argc, char **argv) {
    uint64_t sketchsize = 1024;
    int sketchtype = 1;
    bool u32 = false;
    int f32 = 0;
    bool ip32 = false;
    std::string outpref;
    int nthreads = 1;
    for(int c;(c = getopt(argc, argv, "p:o:S:UPqBHPufh?")) >= 0;) { switch(c) {
        case 'p': nthreads = std::atoi(optarg); break;
        case 'S': sketchsize = std::strtoull(optarg, nullptr, 10); break;
        case 'B': sketchtype = 0; break;
        case 'q': sketchtype = -1; break;
        case 'u': u32 = true; break;
        case 'f': f32 = true; break;
        case 'H': f32 = -1; break;
        case 'U': f32 = -2; break;
        case 'o': outpref = optarg; break;
        case 'P': ip32 = true; break;
        case '?': case 'h': return wsketchusage();
    }}
    OMP_ONLY(omp_set_num_threads(std::max(nthreads, 1));)
    auto diff = argc - optind;
    if(diff < 1 || diff > 3) {
        std::fprintf(stderr, "Required: between one and three positional arguments. All flags must come before positional arguments. Diff: %d\n", diff);
        return wsketchusage();
    }
    if(outpref.empty()) {
        outpref = argv[optind];
    }
    if(diff == 3) {
        //std::fprintf(stderr, "Getting CSR hashes\n");
        auto mhrs = wmh_from_file_csr(argv[optind], argv[optind + 1], argv[optind + 2], sketchsize, sketchtype, f32, u32, ip32);
        if(mhrs.empty()) {
            THROW_EXCEPTION(std::runtime_error("No sketches found in file; this suggests there was an error."));
        }
        const uint64_t nsketches = mhrs.size();
        std::string of = outpref + ".sampled.indices.stacked." + std::to_string(nsketches) + "." + std::to_string(sketchsize) + ".i64";
        std::FILE *fp = bfopen(of.data(), "wb");
        if(fp == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to open " + of));
        for(size_t i = 0; i < nsketches; ++i) {
            if(std::fwrite(std::get<2>(mhrs[i]).data(), sizeof(uint64_t), std::get<2>(mhrs[i]).size(), fp) != std::get<2>(mhrs[i]).size()) {
                THROW_EXCEPTION(std::runtime_error("Failed to write MH identifiers to disk."));
            }
        }
        std::fclose(fp);
        of =  outpref + ".sampled.regs.stacked." + std::to_string(nsketches) + "." + std::to_string(sketchsize) + ".f" + std::to_string(sizeof(RegT) * 8);
        fp = bfopen(of.data(), "wb");
        if(fp == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to open " + of));
        checked_fwrite(fp, &nsketches, sizeof(nsketches));
        checked_fwrite(fp, &sketchsize, sizeof(sketchsize));
        {
            std::vector<double> cards(nsketches);
            std::transform(mhrs.begin(), mhrs.end(), cards.data(), [](const auto &x) {return x.total_weight();});
            if(size_t nwritten = std::fwrite(cards.data(), sizeof(double), cards.size(), fp); nwritten != cards.size()) {
                THROW_EXCEPTION(std::runtime_error(std::string("Failed to write cards to disk; wrote ") + std::to_string(nwritten) + " instead of " + std::to_string(cards.size())));
            }
        }
        const auto sz = std::get<0>(mhrs.front()).size();
        for(size_t i = 0; i < nsketches; ++i) {
            if(const size_t fc = std::fwrite(std::get<0>(mhrs[i]).data(), sizeof(RegT), sz, fp); fc != sz) {
                std::fprintf(stderr, "Failed to write total of %zu to disk (wrote %zu)\n", sz, fc);
                THROW_EXCEPTION(std::runtime_error("Failed to write MH registers to disk."));
            }
        }
        std::fclose(fp);
        of = outpref + ".sampled.hashes.stacked." + std::to_string(nsketches) + "." + std::to_string(sketchsize) + ".i64";
        fp = bfopen(of.data(), "wb");
        if(fp == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to open " + of));
        for(size_t i = 0; i < nsketches; ++i) {
            auto &regs = std::get<1>(mhrs[i]);
            auto sz = regs.size();
            if(size_t fc = std::fwrite(regs.data(), sizeof(uint64_t), sz, fp); fc != sz) {
                std::fprintf(stderr, "Failed to write total of %zu to disk (wrote %zu)\n", sz, fc);
                THROW_EXCEPTION(std::runtime_error("Failed to write MH registers to disk."));
            }
        }
        std::fclose(fp);
        fp = bfopen((outpref + ".sampled.info.txt").data(), "wb");
        if(fp == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to open " + of));
        static constexpr const char *fmtstr = nlfmt<std::decay_t<decltype(mhrs.front().total_weight())>>;
        for(size_t i = 0; i < nsketches; std::fprintf(fp, fmtstr, mhrs[i++].total_weight()));
        std::fclose(fp);
        return 0;
    }
    SimpleMHRet mh;
#ifndef NDEBUG
    for(int i = 0; i < argc; ++i) std::fprintf(stderr, "Arg %d/%d is %s (optind = %d)\n", i, argc, argv[i], optind);
#endif
    std::string destpath;
    if(diff == 2) destpath = argv[optind + 1];
    auto [sigs, indices, ids, total_weight] = std::move(wmh_from_file(argv[optind], destpath, sketchsize, sketchtype, f32, u32).tup());

    write_container(indices, outpref + ".sampled.indices.u64");
    const std::string suffix = (sizeof(RegT) == 4  ? "32" : sizeof(RegT) == 8 ? "64": sizeof(RegT) == 16 ? "128": "UNKNOWN");
    if(suffix == "UNKNOWN") THROW_EXCEPTION(std::runtime_error("RegT is not float32, float64, or long double...."));
    std::string sigpath = outpref + ".sampled.hashes.f" + suffix;
    std::FILE *sigfp = bfopen(sigpath.data(), "wb");
    if(sigfp == nullptr) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open sigpath ") + sigpath));
    double tw = total_weight;
    std::fwrite(&tw, sizeof(tw), 1, sigfp);
    std::fwrite(sigs.data(), sizeof(RegT), sigs.size(), sigfp);
    std::fclose(sigfp);
    write_container(ids, outpref + ".sampled.ids.u64");
    std::ofstream ofs(outpref + ".sampled.tw.txt");
    std::string msg = std::string("Total weight: ") + std::to_string(total_weight) + ";" + argv[optind];
    if(optind + 1 < argc) msg += std::string(";") + argv[optind + 1];
    msg += ';' + (f32 == 1 ? 'f' : f32 == 0 ?'d': 'H') + ';' + (u32 ? 'W': 'L');
    msg += '\n';
    ofs <<  msg;
    std::cerr << msg;
    return 0;
}


} // namespace wmh
