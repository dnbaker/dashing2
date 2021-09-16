#include "d2.h"
#include <mio.hpp>
#include "hash.h"
#include "bonsai/encoder.h"

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

std::vector<flat_hash_map<uint64_t, uint64_t>> get_results(bns::Encoder<bns::score::Lex, uint64_t> &eenc, bns::RollingHasher<uint64_t> &renc, std::vector<std::string> input_files, const flat_hash_map<uint64_t, std::vector<uint64_t>> &kmer2ids) {
    std::vector<flat_hash_map<uint64_t, uint64_t>> res(input_files.size());
    //KSeqHolder kseqs(nthreads);
    OMP_PFOR_DYN
    for(size_t i = 0; i < input_files.size(); ++i) {
        auto &myres = res[i];
        auto func = [&](auto kmer) {
            kmer = maskfn(kmer);
            auto kmeridit = kmer2ids.find(kmer);
            if(kmeridit == kmer2ids.end()) {
                //std::fprintf(stderr, "Not in map, ignoring k-mer\n");
                // We don't have to process any k-mers not in this set!
                return;
            }
            auto it = myres.find(kmer);
            if(it == myres.end()) it = myres.emplace(kmer, 1).first;
            else ++it->second;
        };
        auto path = input_files[i].data();
        if(eenc.k() <= eenc.nremperres64()) {
            eenc.for_each(func, path);
        } else {
            renc.for_each_hash(func, path);
        }
    }
    return res;
}

int contain_main(int argc, char **argv) {
    int nthreads = 1;
    bool binary_output = false;
    char *outpath = 0;
    std::vector<std::string> streamfiles;
    for(int c;(c = getopt(argc, argv, "h?p:")) >= 0;) switch(c) {
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
#ifdef _OPENMP
    if(nthreads > 1) omp_set_num_threads(nthreads);
#endif
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
    dashing2::flat_hash_map<uint64_t, uint32_t> fhs;
    for(size_t i = 0; i < nitems; ++i) {
        uint64_t *ptr = ((uint64_t *)dbptr + 3 + sketchsize * i);
        assert((const char *)ptr < &db.data()[db.size()]);
        for(size_t j = 0; j < sketchsize; ++j) {
            const auto v = ptr[j];
            kmer2ids[v].push_back(i);
            auto it = fhs.find(v);
            if(it == fhs.end()) fhs.emplace(v, 1u);
            else ++it->second;
        }
        fhs.clear();
    }
    std::vector<flat_hash_map<uint64_t, uint64_t>> res = get_results(e64, rh64, streamfiles, kmer2ids);
    std::vector<float> coverage_mat(nitems * streamfiles.size());
    std::vector<float> coverage_stats(nitems * streamfiles.size());
    OMP_PFOR_DYN
    for(size_t i = 0; i < res.size(); ++i) {
        std::unique_ptr<uint32_t[]> matches(new uint32_t[nitems]);
        std::unique_ptr<uint32_t[]> matchsums(new uint32_t[nitems]);
        for(const auto &[key, value]: res[i]) {
            for(const auto refid: kmer2ids.at(key)) {
                ++matches[refid];
                matchsums[refid] += value;
            }
        }
        const double ssiv = 1. / sketchsize;
        const auto cmatptr = &coverage_mat[nitems * i];
        const auto cstatsptr = &coverage_stats[nitems * i];
        for(size_t j = 0; j < nitems; ++j) {
            if(matches[j]) {
                cmatptr[j] = ssiv * matches[j];
                cstatsptr[j] = matchsums[j] / matches[j];
            }
        }
    }
    std::FILE *ofp = outpath ? bfopen(outpath, "w"): stdout;
    if(binary_output) {
        std::fwrite(coverage_mat.data(), sizeof(float), coverage_mat.size(), ofp);
        std::fwrite(coverage_stats.data(), sizeof(float), coverage_mat.size(), ofp);
    } else {
        std::fprintf(ofp, "#Dashing2 contain - a list of coverage %%s for the set of references, + mean coverage levels.\n");
        std::fprintf(stderr, "#Each matrix entry consists of <coverage%%:mean depth of coverage>\n");
        std::fprintf(ofp, "##References:");
        for(size_t i = 0; i < nitems; ++i) {
            std::fputc('\t', ofp);
            std::fwrite(names[i].data(), 1, names[i].size(), ofp);
        }
        std::fputc('\n', ofp);
        for(size_t i = 0; i < nq; ++i) {
            std::fwrite(streamfiles[i].data(), 1, streamfiles[i].size(), ofp);
            const float *cmatptr = &coverage_mat[nitems * i];
            const float *cstatsptr = &coverage_stats[nitems * i];
            for(size_t j = 0; j < nitems; ++j) {
                assert(cmatptr + j < &*coverage_mat.end());
                assert(cstatsptr + j < &*coverage_stats.end());
                std::fprintf(ofp, "\t%0.8g%%:%0.4g", 100. * cmatptr[j], cstatsptr[j]);
            }
            std::fputc('\n', ofp);
        }
    }
    if(ofp != stdout) std::fclose(ofp);
    return 0;
}

}
