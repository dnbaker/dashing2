#include "d2.h"
#include <mio.hpp>
#include "robin_hood.h"
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

std::vector<flat_hash_map<uint64_t, uint64_t>> get_results(bns::Encoder<bns::score::Lex, uint64_t> &eenc, bns::RollingHasher<uint64_t> &renc, std::vector<std::string> input_files, const flat_hash_map<uint64_t, std::vector<uint64_t>> &kmer2ids, int nthreads) {
    std::vector<flat_hash_map<uint64_t, uint64_t>> res(input_files.size());
    KSeqHolder kseqs(nthreads);
    OMP_PFOR_DYN
    for(size_t i = 0; i < input_files.size(); ++i) {
        auto &myres = res[i];
        auto func = [&](auto kmer) {
            auto kmeridit = kmer2ids.find(kmer);
            if(kmeridit == kmer2ids.end()) {
                // We don't have to process any k-mers not in this set!
                return;
            }
            auto it = myres.find(kmer);
            if(it == myres.end()) myres.emplace(kmer, 1);
            else ++it->second;
        };
        const int tid = OMP_ELSE(omp_get_thread_num(), 1);
        if(eenc.k() <= eenc.nremperres64()) {
            eenc.for_each(func, input_files[i].data(), &kseqs[tid]);
        } else {
            renc.for_each_hash(func, input_files[i].data(), &kseqs[tid]);
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
    std::string databasefile;
    switch(diff) {
        case 0: return contain_usage(); break;
        case 1: databasefile = argv[optind];
            if(streamfiles.empty()) {
                streamfiles.push_back("/dev/stdin");
                std::fprintf(stderr, "No input files provided; this defaults to stdin. If this is in error, cancel dashing2.\n");
            }
        break;
        default: streamfiles.insert(streamfiles.end(), argv + optind + 1, argv + argc); break;
    }
    const size_t nq = streamfiles.size();
    mio::mmap_source db(databasefile);
    const void *dbptr = (void *)db.data();
    bns::InputType rht = static_cast<bns::InputType>(*((uint32_t *)dbptr) & 0xFF);
    const bool canon = *(uint32_t *)dbptr & 0x100;

    const uint32_t sketchsize = ((const uint32_t *)dbptr)[1];
    const uint32_t k = ((const uint32_t *)dbptr)[2];
    const uint32_t w = ((const uint32_t *)dbptr)[3];
    if((db.size() - 16) % sketchsize) THROW_EXCEPTION(std::runtime_error("Database corrupted (not a multiple of uint64_t size). Regenerate?"));
    std::vector<std::string> names;
    {
        if(bns::isfile(databasefile + ".names.txt")) {
            std::ifstream ifs(databasefile + ".names.txt");
            for(std::string line;std::getline(ifs, line); names.emplace_back(line));
        } else {
            const size_t v = ((db.size() - 16) / sketchsize);
            while(names.size() < v) names.push_back(std::to_string(v));
        }
    }
    const size_t nitems = names.size();
    if(nitems != ((db.size() - 16) / sketchsize)) THROW_EXCEPTION(std::runtime_error("Database corrupted; wrong number of names."));
    bns::Spacer sp(k, w);
    bns::Encoder<bns::score::Lex, uint64_t> e64(sp, nullptr, canon);
    bns::RollingHasher<uint64_t> rh64(k, canon, rht, w);
    flat_hash_map<uint64_t, std::vector<uint64_t>> kmer2ids;
    for(size_t i = 0; i < nitems; ++i) {
        uint64_t *ptr = ((uint64_t *)dbptr + 2 + sketchsize * i);
        for(size_t j = 0; j < sketchsize; ++j)
            kmer2ids[ptr[j]].push_back(i);
    }
    std::vector<flat_hash_map<uint64_t, uint64_t>> res = get_results(e64, rh64, streamfiles, kmer2ids, nthreads);
    std::vector<float> coverage_mat(nitems * streamfiles.size());
    std::vector<float> coverage_stats(nitems * streamfiles.size());
    OMP_PFOR_DYN
    for(size_t i = 0; i < res.size(); ++i) {
        std::vector<uint32_t> matches(nitems);
        std::vector<uint32_t> matchsums(nitems);
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
        std::fprintf(stderr, "Each matrix entry consists of <coverage%%:mean depth of coverage>\n");
        std::fprintf(ofp, "#References:");
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
                std::fprintf(ofp, "\t%0.8g%%:%0.4g", cmatptr[j], cstatsptr[j]);
            }
            std::fputc('\n', ofp);
        }
    }
    if(ofp != stdout) std::fclose(ofp);
    return 0;
}

}
