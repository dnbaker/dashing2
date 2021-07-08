#include "wsketch.h"
#include "enums.h"
#include "d2.h"
#include "sketch/bmh.h"
namespace dashing2 {

template<typename BMH, typename FT, typename IT>
SimpleMHRet minwise_det(const FT *weights, const IT *indices, size_t n, size_t m) {
    BMH h(m, true);
    for(size_t i = 0; i < n; ++i) {
        // We use the offset as a first ID, and then we convert these to the
        // relevant IDs sampled at the end
        h.add(i, weights[i]);
    }
    std::vector<uint64_t> ids(m);
    std::vector<RegT> regs;
    h.finalize();
    auto &hids = h.ids();
    std::transform(hids.begin(), hids.end(), ids.begin(), [ind=indices](auto x) {return ind[x];});
    return std::make_tuple(h.template to_sigs<RegT>(), hids, ids, double(h.total_weight()));
}
template<typename FT, typename IT>
SimpleMHRet minhash(const FT *weights, const IT *indices, size_t n, size_t m, bool usepmh=true) {
    if(usepmh)
        return minwise_det<ProbMinHash>(weights, indices, n, m);
    return minwise_det<BagMinHash>(weights, indices, n, m);
}

#define SIG(T1, T2, P1, P2) \
SimpleMHRet minhash##T1##T2(const P1 *weights, const P2 *indices, size_t n, size_t m, bool usepmh=true);

#define DMH(T1, T2, P1, P2) \
SimpleMHRet minhash##T1##T2(const P1 *weights, const P2 *indices, size_t n, size_t m, bool usepmh) { \
    return minhash(weights, indices, n, m, usepmh);}\

#define BODY(MACRO) \
MACRO(f64, u64, double, uint64_t)\
MACRO(f64, u32, double, uint32_t)\
MACRO(f32, u64, float, uint64_t)\
MACRO(f32, u32, float, uint32_t)

BODY(SIG) BODY(DMH)
#undef SIG
#undef DMH

struct FReader{
    std::string path_;
    std::FILE *fp_;
    FReader(std::string path): path_(path), fp_(nullptr) {
        std::string c = "cat ";
        if(endswith(path, ".gz"))
            c = "gzip -dc ";
        else if(endswith(path, ".xz"))
            c = "xz -dc ";
        else if(endswith(path, ".bz2"))
            c = "bzip2 -dc ";
        if((fp_ = ::popen((c + path).data(), "r")) == nullptr)
            THROW_EXCEPTION(std::runtime_error(std::string("Failed to run cmd '") + c + path + "'"));
    }
    template<typename VT>
    std::vector<VT> getvec() {
        std::vector<VT> ret;
        int rc = 0;
        while(!std::feof(fp_)) {
            VT v;
            if(std::fread(&v, sizeof(v), 1, fp_) != 1u) {rc = 1; break;}
            ret.push_back(v);
        }
        if(rc) throw std::runtime_error("Failed to read from file in expected increments");
        return ret;
    }
    static inline bool endswith(std::string lhs, std::string rhs) {
        return std::equal(rhs.begin(), rhs.end(), &lhs[lhs.size() - rhs.size()]);
    }
    void clear() {
        if(fp_) {::pclose(fp_); fp_ = nullptr;}
    }
    ~FReader() {
        if(fp_) {clear();}
    }
};

SimpleMHRet wmh_from_file(std::string idpath, std::string cpath, size_t sksz, bool usepmh, bool usef32=false, bool wordids=false) {
    // usef32 means use float32 instead of float64 for weights
    // default is f64
    // For IDs, default is 64 bits
    // 
    FReader lhs(cpath), rhs(idpath);
#define PERF(TL, TR) do {auto lvec = lhs.getvec<TL>(); auto rvec = rhs.getvec<TR>();\
            assert(lvec.size() == rvec.size());\
            return minhash(lvec.data(), rvec.data(), lvec.size(), sksz, usepmh);\
        } while(0)
#define PERF2(LT) do {\
        if(wordids) {\
            PERF(LT, uint32_t);\
        } else {\
            PERF(LT, uint64_t);\
        } \
    } while(0)
    if(usef32) {
        PERF2(float);
    } else {
        PERF2(double);
    }
    throw std::runtime_error("This should never happen");
#undef PERF
#undef PERF2
}

int wsketchusage() {
    std::fprintf(stderr, "Sketch raw IDs, with optional weights added\n"
                         "Usage: dashing2 binsketch [input.bin] [input.weights.bin]\n"
                         "-u: Read 32-bit identifiers from input.bin rather than 64-bit (Default: 64-bit integers)\n"
                         "-f: Read 32-bit floating point weights from [input.weights.bin] (Default: double)\n"
                         "-B: Sketch with BagMinHash (Default: Uses ProbMinhash)\n"
                         "-o: outprefix. If unset, uses [input.bin]\n"
            );
    return 1;
}
template<typename T>
void write_container(const T &vec, std::string path, size_t nb = 0) {
    if(nb == 0) nb = vec.size();
    std::FILE *ifp = std::fopen(path.data(), "wb");
    std::fwrite(vec.data(), sizeof(T), nb, ifp);
    std::fclose(ifp);
}
int wsketch_main(int argc, char **argv) {
    int sketchsize = 1024;
    bool bmhsketch = true;
    bool u32 = false, f32 = false;
    std::string outpref;
    for(int c;(c = getopt(argc, argv, "o:B:S:uf")) >= 0;) { switch(c) {
        case 'S': sketchsize = std::atoi(optarg); break;
        case 'B': bmhsketch = true; break;
        case 'u': u32 = true; break; case 'f': f32 = true; break;
        case 'o': outpref = optarg; break;
        case 'h': return wsketchusage();
    }}
    if(argc + 2 != optind) {
        std::fprintf(stderr, "Required: two positional arguments. All flags must come before positional arguments.\n");
        std::fprintf(stderr, "Example: 'dashing2 wsketch -o g1.k31.k64 g1.fastq.k31.kmerset64 g1.fastq.k31.kmercounts.f64'.\n");
        return wsketchusage();
    }
    if(outpref.empty()) outpref = argv[optind];
    SimpleMHRet mh(wmh_from_file(argv[optind], argv[optind + 1], sketchsize, !bmhsketch, f32, u32));
    auto &[sigs, indices, ids, total_weight] = mh;
    std::string isamples = outpref + ".sampled.indices.u64";
    
    write_container(indices, outpref + ".sampled.indices.u64");
    write_container(sigs, outpref + std::string(".sampled.hashes.f") + (sizeof(RegT) == 4  ? "32" : sizeof(RegT) == 8 ? "64": sizeof(RegT) == 16 ? "128": "UNKNOWN"));
    write_container(ids, outpref + ".sampled.ids.u64");
    write_container(std::string("Total weight: ") + std::to_string(total_weight),
                    outpref + ".sampled.tw.txt");
    return 0;
}


} // namespace wmh
