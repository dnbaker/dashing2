#include "d2.h"
#include <mio.hpp>
#include "hash.h"
#include "bonsai/encoder.h"
#include "fmt/format.h"
#include "FastxParser.hpp"

namespace dashing2 {
template<typename Key, typename V, typename Hash>
using flat_hash_map = robin_hood::unordered_flat_map<Key, V, Hash>;

int contain_usage() {
    std::fprintf(stderr, "Usage: dashing2 contain <flags> database.kmers <input.fq> <input2.fq>...\n"
                         "This application is inspired by mash screen.\n"
                         "Flags:\n"
                         "-h: help\n"
                         "-p: set number of threads. [1]\n"
                         "-o: Set output [stdout]\n"
                         "-b: Emit binary output instead of human-readable.\n"\
                         "-F: read input paths from file at <arg>\n"\
                );
    return EXIT_FAILURE;
}

INLINE flat_hash_map<uint64_t, uint64_t> & operator+=(flat_hash_map<uint64_t, uint64_t> &lhs, const flat_hash_map<uint64_t, uint64_t> &rhs) {
    for(const auto &pair: rhs) {
        auto lit = lhs.find(pair.first);
        if(lit != lhs.end()) lit->second += pair.second;
        else lhs.emplace(pair.first, pair.second);
    }
    return lhs;
}

std::vector<flat_hash_map<uint64_t, uint64_t>> get_results(bns::Encoder<bns::score::Lex, uint64_t> &eenc, bns::RollingHasher<uint64_t> &renc, std::vector<std::string> input_files, const flat_hash_map<uint64_t, std::vector<uint64_t>> &kmer2ids, const uint64_t maxkmer, const uint64_t minkmer) {
    std::vector<flat_hash_map<uint64_t, uint64_t>> res(input_files.size());
    //KSeqHolder kseqs(nthreads);
    OMP_PFOR_DYN
    for(size_t i = 0; i < input_files.size(); ++i) {
        auto &myres = res[i];
        auto func = [&](auto kmer) {
            kmer = maskfn(kmer);
            if(kmer < minkmer || kmer > maxkmer) return;
            auto kmeridit = kmer2ids.find(kmer);
            if(kmeridit == kmer2ids.end()) return;
            auto it = myres.find(kmer);
            if(it == myres.end()) myres.emplace(kmer, 1);
            else ++it->second;
        };
        auto path = input_files[i].data();
        if(eenc.k() <= eenc.nremperres64()) {
            bns::Encoder<bns::score::Lex, uint64_t> mye(eenc);
            mye.for_each(func, path);
        } else {
            bns::RollingHasher<uint64_t> myr(renc);
            myr.for_each_hash(func, path);
        }
    }
    return res;
}

template<typename T>
void par_reduce(T *x, size_t n) {
    const unsigned int ln = static_cast<int>(std::ceil(std::log2(n)));
    for(size_t i = 0; i < ln; ++i) {
        const size_t step_size = size_t(1) << i;
        const size_t sweep_size = (i + 1) << 1;
        const size_t nsweeps = (n + (sweep_size - 1)) / sweep_size;
        OMP_PFOR
        for(size_t j = 0; j < nsweeps; ++j) {
            if(const size_t lh = j * sweep_size, rh = lh + step_size; rh < n)
                x[lh] += x[rh];
        }
    }
}

flat_hash_map<uint64_t, uint64_t> get_results_sf(bns::Encoder<bns::score::Lex, uint64_t> &eenc, bns::RollingHasher<uint64_t> &renc, std::string input_file, const flat_hash_map<uint64_t, std::vector<uint64_t>> &kmer2ids, const uint64_t maxkmer, const uint64_t minkmer, const int nthreads) {
    std::vector<flat_hash_map<uint64_t, uint64_t>> res(nthreads);
    std::vector<std::string> sf;
    for_each_substr([&sf](const auto &x) {sf.push_back(x);}, input_file);
    std::vector<std::thread> threads;
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(sf, nthreads, 1);
    const bool use_direct_encoding = eenc.k() <= eenc.nremperres64();
    parser.start();
    for(size_t i = 0; i < size_t(nthreads); ++i) {
        threads.emplace_back([&,i]() {
            auto func = [minkmer,maxkmer,&kmer2ids,&myres=res[i]](auto kmer) __attribute__((always_inline)) {
                kmer = maskfn(kmer);
                if(kmer < minkmer || kmer > maxkmer) return;
                auto kmeridit = kmer2ids.find(kmer);
                if(kmeridit == kmer2ids.end()) return;
                auto it = myres.find(kmer);
                if(it == myres.end()) myres.emplace(kmer, 1);
                else ++it->second;
            };
            bns::Encoder<bns::score::Lex, uint64_t> mye(eenc);
            bns::RollingHasher<uint64_t> myr(renc);
            for(auto rg = parser.getReadGroup();parser.refill(rg);) {
                if(use_direct_encoding) {
                    for(const auto &seq: rg) mye.for_each(func, seq.seq.data(), seq.seq.size());
                } else {
                    for(const auto &seq: rg) myr.for_each_hash(func, seq.seq.data(), seq.seq.size());
                }
            }
        });
    }
    for(auto &t: threads) t.join();
    parser.stop();
    if(nthreads > 1) {
        OMP_ONLY(omp_set_num_threads(nthreads));
    }
    par_reduce(res.data(), res.size());
    OMP_ONLY(omp_set_num_threads(1));
    return std::move(res.front());
}


#ifdef __AVX512F__
INLINE void interleave_512_ps(__m512 &lhs, __m512 &rhs) {
    const __m512 lhret = _mm512_permutex2var_ps(lhs, _mm512_setr_epi32(0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23), rhs);
    const __m512 rhret = _mm512_permutex2var_ps(lhs, _mm512_setr_epi32(8, 24, 9, 9 + 16, 10, 10 + 16, 11, 11 + 16, 12, 12 + 16, 13, 13 + 16, 14, 14 + 16, 15, 15 + 16), rhs);
    lhs = lhret;
    rhs = rhret;
}
#endif
#ifdef __AVX2__
INLINE void interleave_256_ps(__m256 &lhs, __m256 &rhs) {
    const __m256i select = _mm256_setr_epi32(0, 4, 1, 5, 2, 6, 3, 7);
    lhs = _mm256_permutevar8x32_ps(_mm256_permute2f128_ps(lhs, rhs, 0b00100000), select);
    rhs = _mm256_permutevar8x32_ps(_mm256_permute2f128_ps(lhs, rhs, 0b00110001), select);
}
#endif

int contain_main(int argc, char **argv) {
    int nthreads = 1;
    bool binary_output = false;
    char *outpath = 0;
    std::vector<std::string> streamfiles;
    for(int c;(c = getopt(argc, argv, "bh?p:o:F:")) >= 0;) switch(c) {
        case 'h': case '?': return contain_usage();
        case 'p': nthreads = std::max(std::atoi(optarg), 1); break;
        case 'b': binary_output = true; break;
        case 'o': outpath = optarg; break;
        case 'F': {
            std::ifstream ifs(optarg);
            for(std::string line;std::getline(ifs, line);)
                streamfiles.emplace_back(line);
            break;
        }
    }
    nthreads = std::max(nthreads, 1);
    auto diff = argc - optind;
    if(diff == 0) return contain_usage();
    std::string databasefile = argv[optind];
    streamfiles.insert(streamfiles.end(), argv + optind + 1, argv + argc);
    if(streamfiles.empty()) {
        streamfiles.push_back("/dev/stdin");
        std::fprintf(stderr, "No input files provided; this defaults to stdin. If this is in error, cancel dashing2.\n");
    }
    const size_t nq = streamfiles.size();
    static constexpr size_t headerlen = 24;
    mio::mmap_source db(databasefile);
    assert(db.size() >= headerlen);
    const void *dbptr = (void *)db.data();
    bns::InputType rht = static_cast<bns::InputType>(*((uint32_t *)dbptr) & 0xFF);
    const bool canon = *(uint32_t *)dbptr & 0x100;

    const uint32_t sketchsize = ((const uint32_t *)dbptr)[1];
    const uint32_t k = ((const uint32_t *)dbptr)[2];
    const uint32_t w = ((const uint32_t *)dbptr)[3];
    const uint64_t seed = ((const uint64_t *)dbptr)[2];
    seed_mask(seed);
    if((db.size() - headerlen) % sketchsize) THROW_EXCEPTION(std::runtime_error("Database corrupted (not a multiple of uint64_t size). Regenerate?"));
    std::vector<std::string> names;
    {
        if(bns::isfile(databasefile + ".names.txt")) {
            std::ifstream ifs(databasefile + ".names.txt");
            for(std::string line;std::getline(ifs, line); names.emplace_back(line));
        } else {
            const size_t v = ((db.size() - headerlen) / sketchsize);
            while(names.size() < v) names.push_back(std::to_string(v));
        }
    }
    const size_t nitems = names.size();
    if(nitems != ((db.size() - headerlen) / sketchsize / sizeof(uint64_t))) THROW_EXCEPTION(std::runtime_error("Database corrupted; wrong number of names."));
    bns::Spacer sp(k, w);
    bns::Encoder<bns::score::Lex, uint64_t> e64(sp, nullptr, canon);
    bns::RollingHasher<uint64_t> rh64(k, canon, rht, w);
    flat_hash_map<uint64_t, std::vector<uint64_t>> kmer2ids;
    for(size_t i = 0; i < nitems; ++i) {
        uint64_t *ptr = ((uint64_t *)dbptr + 3 + sketchsize * i);
        assert((const char *)ptr < &db.data()[db.size()]);
#if _OPENMP >= 201307L
        #pragma omp simd
#endif
        for(size_t j = 0; j < sketchsize; ++j) {
            kmer2ids[ptr[j]].push_back(i);
        }
    }
    // We can quickly filter out k-mers by using min/max k-mers
    // to eliminate anything > or <
    // These only take a comparison instruction and a jump, which means there's no cache miss.
    // We also compute the minimum over the set of keys, so it's on the order of the number of unique sampled k-mers
    // instead of the total number of sampled k-mers (num_samples * sample_size).
    uint64_t maxkmer = 0, minkmer = std::numeric_limits<uint64_t>::max();
    for(const auto &pair: kmer2ids) {
        maxkmer = std::max(maxkmer, pair.first);
        minkmer = std::min(minkmer, pair.first);
    }
    std::vector<flat_hash_map<uint64_t, uint64_t>> res;
    if(nthreads > 1 && nq < size_t(nthreads)) {
        for(const auto &sf: streamfiles) {
            res.emplace_back(get_results_sf(e64, rh64, sf, kmer2ids, maxkmer, minkmer, nthreads));
        }
    } else {
        res = get_results(e64, rh64, streamfiles, kmer2ids, maxkmer, minkmer);
    }
    const size_t tablesize = nitems * streamfiles.size();
    const size_t table2size = tablesize * 2;
    std::vector<float, sketch::Allocator<float>> coverage_mat(table2size);
    float *const coverage_stats = coverage_mat.data() + tablesize;
    const double ssiv = 1. / sketchsize;
    OMP_PFOR_DYN
    for(size_t i = 0; i < res.size(); ++i) {
        uint32_t *matches = static_cast<uint32_t *>(std::calloc(nitems * 2, sizeof(uint32_t)));
        uint32_t *matchsums = matches + nitems;
        for(const auto &[key, value]: res[i]) {
            for(const auto refid: kmer2ids[key]) {
                ++matches[refid];
                matchsums[refid] += value;
            }
        }
        const auto cmatptr = &coverage_mat[nitems * i];
        const auto cstatsptr = &coverage_stats[nitems * i];
#ifdef _OPENMP
        #pragma omp simd
#endif
        for(size_t j = 0; j < nitems; ++j) {
            if(matches[j]) {
                cmatptr[j] = ssiv * matches[j];
                cstatsptr[j] = matchsums[j] / matches[j];
            }
        }
        std::free(matches);
    }
    std::FILE *ofp = outpath ? bfopen(outpath, "w"): stdout;
    if(binary_output) {
        uint64_t tmparr[2] {nitems, res.size()};
        checked_fwrite(tmparr, sizeof(tmparr), 1, ofp);
        checked_fwrite(coverage_mat.data(), sizeof(float), coverage_mat.size(), ofp);
    } else {
        // TODO: Parallize results formatting
        fmt::print(ofp, "#Dashing2 contain - a list of coverage %%s for the set of references, + mean coverage levels.\n"
                        "#Each matrix entry consists of <coverage%%:mean depth of coverage>\n"
                        "##References:");
        for(size_t i = 0; i < nitems; ++i)
            fmt::print(ofp, "\t{}", names[i]);
        fmt::print(ofp, "\n");
        for(size_t i = 0; i < nq; ++i) {
            fmt::print(ofp, streamfiles[i]);
            const float *cmatptr = &coverage_mat[nitems * i];
            const float *cstatsptr = &coverage_stats[nitems * i];
            size_t j = 0;
#if __AVX512F__
            for(;nitems - j >= 16;j += 16) {
                __m512 matd = _mm512_mul_ps(_mm512_loadu_ps(cmatptr + j), _mm512_set1_ps(100.f));
                __m512 statd = _mm512_loadu_ps(cstatsptr + j);
                interleave_512_ps(matd, statd);
                fmt::print("\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}",
                           matd[0], matd[1], matd[2], matd[3], matd[4], matd[5], matd[6], matd[7], matd[8], matd[9], matd[10], matd[11], matd[12], matd[13], matd[14], matd[15], statd[0], statd[1], statd[2], statd[3], statd[4], statd[5], statd[6], statd[7], statd[8], statd[9], statd[10], statd[11], statd[12], statd[13], statd[14], statd[15]
                );
            }
#elif __AVX2__
            for(;nitems - j >= 8;j += 8) {
                __m256 matd = _mm256_mul_ps(_mm256_loadu_ps(cmatptr + j), _mm256_set1_ps(100.f));
                __m256 statd = _mm256_loadu_ps(cstatsptr + j);
                interleave_256_ps(matd, statd);
                fmt::print("\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}\t{:0.6g}%:{}",
                        matd[0], matd[1], matd[2], matd[3], matd[4], matd[5], matd[6], matd[7], statd[0], statd[1], statd[2], statd[3], statd[4], statd[5], statd[6], statd[7]
                );
            }
#endif
            for(;j < nitems; ++j) {
                fmt::print(ofp, "\t{:0.6g}%:{}", 100.f * cmatptr[j], cstatsptr[j]);
            }
            std::fputc('\n', ofp);
        }
    }
    if(ofp != stdout) std::fclose(ofp);
    return 0;
}

}
