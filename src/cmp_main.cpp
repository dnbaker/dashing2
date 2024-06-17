#include "cmp_main.h"
#include "sketch_core.h"
#include "options.h"
#include "refine.h"
#include <type_traits>

namespace dashing2 {

#define CMP_OPTS(name) \
static option_struct name[] = {\
    LO_FLAG("presketched", OPTARG_PRESKETCHED, presketched, true)\
    SHARED_OPTS\
}



void cmp_usage() {
    std::fprintf(stderr, "dashing2 cmp <opts> [fastas... (optional)]\n"
                         "We use only m-mers; if w <= k, however, this reduces to k-mers if the -w/--window-size is unspecified.\n"
                         "--presketched\t To compute distances using a pre-sketched method (e.g., dashing2 sketch -o path), use this flag and pass in a single positional argument.\n"
                         SHARED_DOC_LINES
    );
}
void load_results(Dashing2DistOptions &opts, SketchingResult &result, const std::vector<std::string> &paths) {
    DBG_ONLY(std::fprintf(stderr, "Loading results using Dashing2Options: %s\n", opts.to_string().data());)
    if(verbosity >= Verbosity::INFO) {
        std::fprintf(stderr, "Loading results. Paths of size %zu", paths.size());
        for(const std::string& path: paths) {
            std::fprintf(stderr, "Path %s/%zd\n", path.data(), &path - paths.data());
        }
    }
    auto &pf = paths.front();
    struct stat st;
    ::stat(pf.data(), &st);
    if(verbosity >= DEBUG) {
        std::fprintf(stderr, "Size of file: %d\n", int(st.st_size));
    }
    const std::string namesf = pf + ".names.txt";
    if(paths.size() == 1) {
        if(verbosity >= INFO) {
            std::fprintf(stderr, "Reading stacked sketches from %s. Names files is %s\n", paths.front().data(), namesf.data());
        }
        if(bns::isfile(namesf)) {
            std::string l;
            for(std::ifstream ifs(namesf);std::getline(ifs, l);) {
                if(l.empty() || l.front() == '#') continue;
                const typename std::string::size_type it = l.find_first_of('\t');
                result.names_.emplace_back(l.substr(0, it));
                if(it != std::string::npos) {
                    char *s = &l[it + 1];
                    result.cardinalities_.emplace_back(std::strtod(s, &s));
                    if(result.cardinalities_.back() <= 0.) result.cardinalities_.back() = 1.; // Set size to 1 if not available
                    if(*s) {
                        DBG_ONLY(std::fprintf(stderr, "Saving kmer count file %s\n", s + 1);)
                        result.kmercountfiles_.push_back(s + 1);
                    }
                }
            }
        }
        //size_t nregs = st.st_size / result.names_.size() / sizeof(RegT);
        std::FILE *fp = bfopen(pf.data(), "rb");
        if(!fp) THROW_EXCEPTION(std::runtime_error(std::string("Failed to open ") + pf));
        uint64_t num_entities;
        if(st.st_size < 8 || std::fread(&num_entities, sizeof(num_entities), 1, fp) != 1u) {
            THROW_EXCEPTION(std::runtime_error(std::string("Failed to read num_entities from file ") + pf + " of size " + std::to_string(st.st_size)));
        }
        // Read in --
        // # of entities (64-bit integer)
        if(verbosity >= DEBUG) {
            std::fprintf(stderr, "%zu entities in file at %s\n", size_t(num_entities), pf.data());
        }
        // sketch size (64-bit integer)
        uint64_t sketchsize;
        if(std::fread(&sketchsize, sizeof(sketchsize), 1, fp) != 1) {
            THROW_EXCEPTION(std::runtime_error(std::string("Failed to read sketch size from file ") + pf + " of size " + std::to_string(st.st_size)));
        }
        opts.sketchsize_ = sketchsize;
        if(result.names_.empty()) {
            result.names_.resize(num_entities);
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic, 64)
#endif
            for(size_t i = 0; i < num_entities; ++i)
                result.names_[i] = std::to_string(i);
        }
        assert(result.cardinalities_.empty() || result.cardinalities_.size() == num_entities || !std::fprintf(stderr, "card size: %zu. Number expected: %zu\n", result.cardinalities_.size(), size_t(num_entities)));
        result.cardinalities_.resize(num_entities);
        // num_entities * sizeof(double) for the cardinalitiy
        if(std::fread(result.cardinalities_.data(), sizeof(double), result.cardinalities_.size(), fp) != result.cardinalities_.size())
            THROW_EXCEPTION(std::runtime_error("Failed to read cardinalities from disk"));
        std::fclose(fp);
        const size_t offset = (num_entities + 2) * sizeof(uint64_t);
        assert((st.st_size - offset) % sizeof(RegT) == 0);
        result.signatures_.assign(pf, offset, (st.st_size - offset) / sizeof(RegT));
    } else { // Else, we have to load sketches from each file -> The vector paths stores the paths to the precomputed sketches
        if(verbosity >= Verbosity::INFO ) {
            std::fprintf(stderr, "Parsing in data from file\n");
        }
        const size_t npaths = paths.size();
        result.nperfile_.resize(npaths);
        auto &fsizes = result.nperfile_;
        std::vector<size_t> csizes(fsizes.size() + 1);
        struct stat st;
        for(size_t i = 0; i < npaths; ++i) {
            if(::stat(paths[i].data(), &st)) {
                std::fprintf(stderr, "File does not exist at %s/%zu\n", paths[i].data(), i);
                std::exit(EXIT_FAILURE);
            }
            const size_t mysz = st.st_size - 8; // 8 bytes for the cardinality (for sketches/kmer sets), length of the sequence for minimizer sequences
            if(verbosity >= Verbosity::DEBUG) {
                std::fprintf(stderr, "Checking file size for %zu/%zu: %zd\n", i, npaths, mysz);
            }
            assert(mysz % sizeof(RegT) == 0);
            const auto nelem = mysz / sizeof(RegT);
            fsizes[i] = nelem;
            csizes[i + 1] = csizes[i] + nelem;
        }
        result.nqueries(npaths);
        // It's even if the items are actually sketches
        // And there are the same number of them per file
        const bool even = opts.kmer_result_ <= FULL_SETSKETCH &&
                 std::all_of(fsizes.begin() + 1, fsizes.end(), [r=fsizes.front()](auto x) {return x == r;});
        if(verbosity >= Verbosity::INFO) {
            std::fprintf(stderr, "[%s:%d] File sizes are %s\n", __FILE__, __LINE__, even ? "even": "uneven");
        }
        if(even) {
            if(fsizes.empty()) {
                throw std::runtime_error("fsizes are empty but should not be.");
            }
            opts.sketchsize_ = fsizes.front();
        }
        std::fprintf(stderr, "Sketchsize is now %zd\n", size_t(opts.sketchsize_));
        const size_t totalsize = csizes.back();
        std::fprintf(stderr, "Resizing signatures to size %zu\n", totalsize);
        result.signatures_.resize(totalsize); // Account for the size of the sketch registers
        if(even) {
            if(totalsize % opts.sketchsize_) {
                throw std::runtime_error(std::string("sanity check: totalsize should be divisible if we get here. Totalsize ") + std::to_string(totalsize) + ", " + std::to_string(opts.sketchsize_));
            }
            result.cardinalities_.resize(totalsize / opts.sketchsize_);
        } else {
            std::fprintf(stderr, "Warning: uneven file sizes. This is expected for hash sets but not sketches.");
        }
        if(verbosity >= Verbosity::INFO) {
            std::fprintf(stderr, "[%s:%d] Resized signatures\n", __FILE__, __LINE__);
        }
        if(bns::isfile(pf + ".kmerhashes.u64")) {
            DBG_ONLY(std::fprintf(stderr, "Loading k-mer hashes, too\n");)
            result.kmers_.resize(result.signatures_.size());
        }
        if(bns::isfile(pf + ".kmercounts.f64")) {
            DBG_ONLY(std::fprintf(stderr, "Loading k-mer counts, too\n");)
            result.kmercounts_.resize(result.signatures_.size());
        }
        assert(npaths == result.cardinalities_.size());
        if(verbosity >= Verbosity::INFO) {
            std::fprintf(stderr, "[%s:%d] About to load from file\n", __FILE__, __LINE__);
        }
        OMP_PFOR_DYN
        for(size_t i = 0; i < npaths; ++i) {
            if(verbosity >= Verbosity::INFO) {
                std::fprintf(stderr, "Loading from file %zd/%s\n", i, paths[i].data());
            }
            auto &path = paths[i];
            std::FILE *ifp = bfopen(path.data(), "rb");
            if(even && result.cardinalities_.size() > i) {
                std::fread(&result.cardinalities_[i], sizeof(double), 1, ifp);
            }
            if(std::fread(&result.signatures_[csizes[i]], sizeof(RegT), fsizes[i], ifp) != fsizes[i]) {
                std::fprintf(stderr, "Failed to read at path %s\n", path.data());
                std::fclose(ifp);
                std::exit(1);
            }
            std::fclose(ifp);
            if(std::string p(path + ".kmerhashes.u64"); bns::isfile(p)) {
                if((ifp = bfopen(p.data(), "rb")) == nullptr) THROW_EXCEPTION(std::runtime_error(std::string("Could not open path ") + p));
                if(std::fread(&result.kmers_[csizes[i]], sizeof(uint64_t), fsizes[i], ifp) != fsizes[i]) {
                    std::fprintf(stderr, "Failed to read at path %s\n", p.data());
                    std::fclose(ifp);
                    std::exit(1);
                }
                std::fclose(ifp);
            }
            if(std::string p(path + ".kmercounts.f64"); bns::isfile(p)) {
                if((ifp = bfopen(p.data(), "rb")) == nullptr) THROW_EXCEPTION(std::runtime_error("Failed to open kmer count file"));
                if(std::fread(&result.kmercounts_[csizes[i]], sizeof(double), fsizes[i], ifp) != fsizes[i]) {
                    std::fprintf(stderr, "Failed to read k-mer counts at path %s\n", p.data());
                    std::fclose(ifp);
                    std::exit(1);
                }
                std::fclose(ifp);
            }
        }
        if(verbosity >= Verbosity::INFO) {
            std::fprintf(stderr, "Loaded all sketches from file.\n");
        }
    }
}

int cmp_main(int argc, char **argv) {
    int c;
    int k = -1, w = 0, nt = -1;
    SketchSpace sketch_space = SPACE_SET;
    KmerSketchResultType res = ONE_PERM;
    bool save_kmers = false, save_kmercounts = false, cache = false, use128 = false, canon = true, presketched = false;
    bool exact_kmer_dist = false;
    bool refine_exact = false; // This uses sketching for K-NN graph generation, then uses exact distances for NN refinement
    long double compressed_a = -1.L, compressed_b = -1.L;
    unsigned int count_threshold = 0.;
    double similarity_threshold = -1.;
    size_t cssize = 0, sketchsize = 1024;
    std::string ffile, outfile, qfile;
    int option_index = 0;
    bns::RollingHashingType rht = bns::DNA;
    DataType dt = DataType::FASTX;
    std::string outprefix;
    OutputKind ok = SYMMETRIC_ALL_PAIRS;
    std::string cmpout; // Only used if distances are also requested -> this should hold the output path to file where distances are written to
    bool normalize_bed = false;
    int topk_threshold = -1;
    int truncate_mode = 0;
    int nLSH = 2;
    int by_chrom = false;
    double nbytes_for_fastdists = sizeof(RegT);
    double downsample_frac = 1.;
    bool parse_by_seq = false;
    bool hpcompress = false;
    std::string fsarg;
    Measure measure = SIMILARITY;
    uint64_t seedseed = 0;
    size_t batch_size = 0;
    bool fasta_dedup = false;
    std::string spacing;
    std::vector<std::pair<uint32_t, uint32_t>> compareids; // TODO: consider a sparse mode comparing only pairs of presented genomes.
    // By default, use full hash values, but allow people to enable smaller
    OutputFormat of = OutputFormat::HUMAN_READABLE;
    if(verbosity >= Verbosity::DEBUG) {
        std::fprintf(stderr, "output format should be %s based on value at start\n", to_string(of).data());
    }
    validate_options(argv, std::vector<std::string>{{"presketched"}}); 
    CMP_OPTS(cmp_long_options);
    std::vector<std::string> paths;
    if(verbosity >= Verbosity::INFO) {
        std::fprintf(stderr, "output format should be %s before parsing options \n", to_string(of).data());
    }
    for(;(c = getopt_long(argc, argv, "m:p:k:w:c:f:S:F:Q:o:L:vNs2BPWh?ZJGH", cmp_long_options, &option_index)) >= 0;) {switch(c) {
        case OPTARG_HELP: case '?': case 'h': cmp_usage(); return 1;
        SHARED_FIELDS
    }}
    if(verbosity >= INFO) {
        std::fprintf(stderr, "output format should be %s after parsing options \n", to_string(of).data());
    }
    if(k < 0) k = nregperitem(rht, use128);
    if(compareids.empty()) {
        paths.insert(paths.end(), argv + optind, argv + argc);
    } else if(optind != argc) throw std::runtime_error("CLI paths must be empty to use pairlist mode.");
    std::unique_ptr<std::vector<std::string>> qup;
    if(nt < 0) {
        char *s = std::getenv("OMP_NUM_THREADS");
        if(s) nt = std::max(std::atoi(s), 1);
    }
    OMP_ONLY(omp_set_num_threads(nt));
    if(ffile.size()) {
        std::ifstream ifs(ffile);
        static constexpr size_t bufsize = 1<<18;
        std::unique_ptr<char []> buf(new char[bufsize]);
        ifs.rdbuf()->pubsetbuf(buf.get(), bufsize);
        for(std::string l;std::getline(ifs, l);paths.push_back(l));
    }
    size_t nref = paths.size();
    if(qfile.size()) {
        std::ifstream ifs(qfile);
        for(std::string l;std::getline(ifs, l);paths.push_back(l));
    }
    size_t nq = paths.size() - nref;
    Dashing2Options opts(k, w, rht, sketch_space, dt, nt, use128, spacing, canon, res);
    opts
        .cache_sketches(cache)
        .cssize(cssize)
        .sketchsize(sketchsize)
        .outprefix(outprefix)
        .save_kmercounts(save_kmercounts)
        .save_kmers(save_kmers)
        .parse_by_seq(parse_by_seq)
        .count_threshold(count_threshold)
        .homopolymer_compress_minimizers(hpcompress)
        .seedseed(seedseed)
        .fasta_dedup(fasta_dedup);
    opts.by_chrom_ = by_chrom;
    opts.compressed_a_ = compressed_a;
    opts.compressed_b_ = compressed_b;
    opts.set_sketch_compressed();
    if(hpcompress) {
        if(!opts.homopolymer_compress_minimizers_) THROW_EXCEPTION(std::runtime_error("Failed to hpcompress minimizers"));
    }
    opts.filterset(fsarg);
    // Ensure we pad the number of registers to a multiple of 64 bits.
    opts.bed_parse_normalize_intervals_ = normalize_bed;
    opts.downsample(downsample_frac);
    Dashing2DistOptions distopts(opts, ok, of, nbytes_for_fastdists, truncate_mode, topk_threshold, similarity_threshold, cmpout, exact_kmer_dist, refine_exact, nLSH);
    default_batchsize(batch_size, distopts);
    distopts.measure_ = measure;
    distopts.cmp_batch_size_ = default_batchsize(batch_size, distopts);
    SketchingResult result;
    if(presketched) {
        std::set<std::string> suffixset;
        for(const auto &p: paths) {
            suffixset.insert(p.substr(p.find_last_of('.'), std::string::npos));
        }
        std::string suf;
        if(suffixset.size() != 1) {
            auto joinstrings =  [](const auto& strings) {
                return std::accumulate(std::cbegin(strings), std::cend(strings), std::string{},[](auto x, const auto&y) {return x + ',' + y;});
            };
            std::fprintf(stderr, "Multiple suffixes in set (%s). Picking the first one. Paths: %s\n", joinstrings(suffixset).data(), joinstrings(paths).data());
            suf = paths.front().substr(paths.front().find_last_of('.'), std::string::npos);
        } else suf = *suffixset.begin();
        if(suf == ".bmh") {
            distopts.sspace_ = SPACE_MULTISET;
            distopts.kmer_result(FULL_SETSKETCH);
        } else if(suf == ".pmh") {
            distopts.sspace_ = SPACE_PSET;
            distopts.kmer_result(FULL_SETSKETCH);
        } else if(suf == ".ss" || suf == ".opss") {
            distopts.sspace_ = SPACE_SET;
            distopts.kmer_result(suf == ".opss" ? ONE_PERM: FULL_SETSKETCH);
        } else if(suf == ".kmerset64") {
            distopts.sspace_ = SPACE_SET;
            distopts.kmer_result(FULL_MMER_SET);
            distopts.use128(false);
        } else if(suf == ".kmerset128") {
            distopts.sspace_ = SPACE_SET;
            distopts.kmer_result(FULL_MMER_SET);
            distopts.use128(true);
        } else if(suf == "mmerseq64") {
            distopts.sspace_ = SPACE_SET;
            auto &path = paths.front();
            std::string countg = path.substr(0, path.find_last_of('.')) + "kmercounts.f64";
            if(bns::isfile(countg)) {
                distopts.kmer_result(FULL_MMER_COUNTDICT);
            } else distopts.kmer_result(FULL_MMER_SEQUENCE);
            distopts.use128(false);
        } else if(suf == "mmerseq128") {
            distopts.sspace_ = SPACE_SET;
            auto &path = paths.front();
            std::string countg = path.substr(0, path.find_last_of('.')) + "kmercounts.f64";
            if(bns::isfile(countg)) {
                distopts.kmer_result(FULL_MMER_COUNTDICT);
            } else distopts.kmer_result(FULL_MMER_SEQUENCE);
            distopts.use128(true);
        }
        load_results(distopts, result, paths);
    } else {
        sketch_core(result, distopts, paths, outfile);
        result.nqueries(nq);
    }
    if(verbosity >= Verbosity::EXTREME) {
        std::fprintf(stderr, "before cmp_core, output format is %s\n", to_string(distopts.output_format_).data());
        for(uint32_t i = 0; i < result.total_seqs(); ++i) {
            for(uint32_t j = 0; j < distopts.sketchsize_; ++j) {
                std::fprintf(stderr, "Item %u has %g at idx %u\n", i, result.signatures_[i * distopts.sketchsize_ + j], j);
            }
        }
    }
    cmp_core(distopts, result);
    return 0;
}


size_t default_batchsize(size_t &batch_size, const Dashing2DistOptions &opts) {
    if(batch_size == 0) {
        if(opts.kmer_result_ <= FULL_SETSKETCH) {
            size_t expl2csz =
#ifdef D2_CACHE_SIZE
                D2_CACHE_SIZE;
#else
                0x400000; // 2^22 bytes, ~= 4 million
#endif
            batch_size = std::max(size_t(expl2csz / opts.sketchsize_ / opts.fd_level_), size_t(1));
        } else {
            batch_size = opts.nthreads_;
        }
    }
    if(batch_size > std::max(opts.nthreads_, 1u)) {
        batch_size = opts.nthreads_;
    }
    return batch_size;
}

}
