#pragma once
#ifndef DASHING2_OPTIONS_H__
#define DASHING2_OPTIONS_H__
#include <enums.h>
#include <getopt.h>
#include "dedup_core.h"

namespace dashing2 {
#define LO_ARG(LONG, SHORT) {LONG, required_argument, 0, SHORT},
#define LO_NO(LONG, SHORT) {LONG, no_argument, 0, SHORT},
#define LO_FLAG(LONG, SHORT, VAR, VAL) {LONG, no_argument, (int *)&VAR, VAL},

using option_struct = struct option;

enum OptArg {
    OPTARG1 = 1000,
    OPTARG_BED,
    OPTARG_BIGWIG,
    OPTARG_LEAFCUTTER,
    OPTARG_OUTPREF,
    OPTARG_CMPOUT,
    OPTARG_BINARY_OUTPUT,
    OPTARG_FASTCMP,
    OPTARG_BBIT_SIGS,
    OPTARG_EXACT_KMER_DIST,
    OPTARG_ISZ,
    OPTARG_BED_NORMALIZE,
    OPTARG_PROTEIN,
    OPTARG_PROTEIN6,
    OPTARG_PROTEIN8,
    OPTARG_PROTEIN14,
    OPTARG_HPCOMPRESS,
    OPTARG_DOWNSAMPLE_FRACTION,
    OPTARG_REFINEEXACT,
    OPTARG_MASHDIST,
    OPTARG_SYMCONTAIN,
    OPTARG_CONTAIN,
    OPTARG_ASYMMETRIC_ALLPAIRS,
    OPTARG_SET,
    OPTARG_SPACING,
    OPTARG_RANDOM_SEED,
    OPTARG_FILTERSET,
    OPTARG_PARSEBYSEQ,
    OPTARG_HELP,
    OPTARG_CMP_BATCH_SIZE,
    OPTARG_GREEDY,
    OPTARG_NLSH,
    OPTARG_ENTROPYMIN,
    OPTARG_SIGRAMLIMIT,
    OPTARG_MAXCAND,
    OPTARG_SETSKETCH_AB,
    OPTARG_FASTCMPBYTES,
    OPTARG_FASTCMPSHORTS,
    OPTARG_FASTCMPWORDS,
    OPTARG_FASTCMPNIBBLES,
    OPTARG_FULL_SETSKETCH,
    OPTARG_PAIRLIST,
    OPTARG_USZ,
    OPTARG_DUMMY,
    OPTARG_SEQS_IN_RAM
};

#define SHARED_OPTS \
    LO_ARG("ffile", 'F')\
    LO_ARG("qfile", 'Q')\
    LO_ARG("threads", 'p')\
    LO_ARG("sketchsize", 'S')\
    LO_ARG("cmpout", OPTARG_CMPOUT)\
    LO_ARG("distout", OPTARG_CMPOUT)\
    LO_ARG("cmp-outfile", OPTARG_CMPOUT)\
    LO_ARG("outprefix", OPTARG_OUTPREF)\
    LO_ARG("prefix", OPTARG_OUTPREF)\
    LO_ARG("topk", 'K')\
    LO_ARG("top-k", 'K')\
    LO_ARG("similarity-threshold", 'T')\
    LO_ARG("fastcmp", OPTARG_FASTCMP)\
    LO_ARG("regsize", OPTARG_FASTCMP)\
    LO_ARG("countsketch-size", 'c')\
    LO_ARG("countmin-size", 'c')\
    LO_ARG("window-size", 'w')\
    LO_ARG("kmer-length", 'k')\
    LO_ARG("outfile", 'o')\
    LO_ARG("count-threshold", 'm')\
    LO_ARG("threshold", 'm')\
    LO_FLAG("binary-output", OPTARG_BINARY_OUTPUT, of, OutputFormat::MACHINE_READABLE) \
    LO_FLAG("emit-binary", OPTARG_BINARY_OUTPUT, of, OutputFormat::MACHINE_READABLE) \
    LO_FLAG("binary", OPTARG_BINARY_OUTPUT, of, OutputFormat::MACHINE_READABLE) \
    LO_FLAG("bagminhash", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("prob", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("probs", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("pminhash", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("pmh", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("PMH", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("probminhash", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("bed", OPTARG_BED, dt, DataType::BED)\
    LO_FLAG("bigwig", OPTARG_BIGWIG, dt, DataType::BIGWIG)\
    LO_FLAG("leafcutter", OPTARG_LEAFCUTTER, dt, DataType::LEAFCUTTER)\
    /*LO_FLAG("edit-distance", 'E', sketch_space, SPACE_EDIT_DISTANCE)*/\
    /*LO_FLAG("set", 'H', res, FULL_MMER_SET)*/\
    /*LO_FLAG("exact-kmer-dist", OPTARG_EXACT_KMER_DIST, exact_kmer_dist, true)*/\
    LO_FLAG("bbit-sigs", OPTARG_BBIT_SIGS, truncate_mode, 1)\
    LO_FLAG("intersection", OPTARG_ISZ, measure, INTERSECTION)\
    LO_FLAG("intersection-size", OPTARG_ISZ, measure, INTERSECTION)\
    LO_FLAG("union-size", OPTARG_USZ, measure, UNION_SIZE)\
    LO_FLAG("mash-distance", OPTARG_MASHDIST, measure, POISSON_LLR)\
    LO_FLAG("distance", OPTARG_MASHDIST, measure, POISSON_LLR)\
    LO_FLAG("symmetric-containment", OPTARG_SYMCONTAIN, measure, SYMMETRIC_CONTAINMENT)\
    LO_FLAG("containment", OPTARG_CONTAIN, measure, Measure::CONTAINMENT)\
    LO_FLAG("poisson-distance", OPTARG_MASHDIST, measure, POISSON_LLR)\
    LO_FLAG("compute-edit-distance", OPTARG_MASHDIST, measure, M_EDIT_DISTANCE)\
    LO_FLAG("multiset", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("bagminhash", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("bmh", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("BMH", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("countdict", 'J', res, FULL_MMER_COUNTDICT)\
    LO_FLAG("seq", 'G', res, FULL_MMER_SEQUENCE)\
    LO_FLAG("128bit", '2', use128, true)\
    LO_FLAG("long-kmers", '2', use128, true)\
    LO_FLAG("phylip", OPTARG_PHYLIP, ok, OutputKind::PHYLIP)\
    LO_FLAG("asymmetric-all-pairs", OPTARG_ASYMMETRIC_ALLPAIRS, ok, OutputKind::ASYMMETRIC_ALL_PAIRS)\
    LO_FLAG("asymmetric", OPTARG_ASYMMETRIC_ALLPAIRS, ok, OutputKind::ASYMMETRIC_ALL_PAIRS)\
    LO_FLAG("square", OPTARG_ASYMMETRIC_ALLPAIRS, ok, OutputKind::ASYMMETRIC_ALL_PAIRS)\
    LO_FLAG("seqs-in-ram", OPTARG_SEQS_IN_RAM, seqs_in_memory, 1) \
    LO_ARG("regbytes", OPTARG_FASTCMP)\
    /*LO_ARG("set", 'H')*/\
    /*{"fastcmp-nibbles", no_argument, 0, OPTARG_FASTCMPNIBBLES},*/\
    {"fastcmp-bytes", no_argument, 0, OPTARG_FASTCMPBYTES},\
    {"fastcmp-shorts", no_argument, 0, OPTARG_FASTCMPSHORTS},\
    {"fastcmp-words", no_argument, 0, OPTARG_FASTCMPWORDS},\
    {"save-kmers", no_argument, 0, 's'},\
    {"save-kmercounts", no_argument, 0, 'N'},\
    {"hp-compress", no_argument, 0, OPTARG_HPCOMPRESS},\
    {"refine-exact", no_argument, 0, OPTARG_REFINEEXACT},\
    {"edit-distance", no_argument, 0, 'E'},\
    {"oneperm-setsketch", no_argument, 0, 'Z'},\
    {"oneperm", no_argument, 0, 'Z'},\
    {"one-perm", no_argument, 0, 'Z'},\
    {"oph", no_argument, 0, 'Z'},\
    {"doph", no_argument, 0, 'Z'},\
    {"full", no_argument, 0, OPTARG_FULL_SETSKETCH},\
    {"full-setsketch", no_argument, 0, OPTARG_FULL_SETSKETCH},\
    {"normalize-intervals", no_argument, 0, OPTARG_BED_NORMALIZE},\
    {"protein", no_argument, 0, OPTARG_PROTEIN},\
    {"protein20", no_argument, 0, OPTARG_PROTEIN},\
    {"enable-protein", no_argument, 0, OPTARG_PROTEIN},\
    {"protein6", no_argument, 0, OPTARG_PROTEIN6},\
    {"protein8", no_argument, 0, OPTARG_PROTEIN8},\
    {"protein14", no_argument, 0, OPTARG_PROTEIN14},\
    {"downsample", required_argument, 0, OPTARG_DOWNSAMPLE_FRACTION},\
    {"cache", no_argument, 0, 'W'},\
    {"cache-sketches", no_argument, 0, 'W'},\
    {"no-canon", no_argument, 0, 'C'},\
    {"set", no_argument, 0, OPTARG_SET},\
    {"exact-kmer-dist", no_argument, 0, OPTARG_EXACT_KMER_DIST},\
    {"spacing", required_argument, 0, OPTARG_SPACING},\
    {"seed", required_argument, 0, OPTARG_RANDOM_SEED},\
    {"filterset", required_argument, 0, OPTARG_FILTERSET},\
    {"parse-by-seq", no_argument, 0, OPTARG_PARSEBYSEQ},\
    {"help", no_argument, 0, OPTARG_HELP},\
    {"batch-size", required_argument, 0, OPTARG_CMP_BATCH_SIZE},\
    {"greedy", required_argument, 0, OPTARG_GREEDY},\
    {"nlsh", required_argument, 0, OPTARG_NLSH},\
    {"nLSH", required_argument, 0, OPTARG_NLSH},\
    {"entmin", no_argument, 0, OPTARG_ENTROPYMIN},\
    {"by-chrom", no_argument, (int *)&by_chrom, 1},\
    {"sketch-size-l2", required_argument, 0, 'L'},\
    {"sig-ram-limit", required_argument, 0, OPTARG_SIGRAMLIMIT},\
    {"maxcand", required_argument, 0, OPTARG_MAXCAND},\
    {"setsketch-ab", required_argument, 0, OPTARG_SETSKETCH_AB},\
    {"pairlist", required_argument, 0, OPTARG_PAIRLIST},\
    {"verbose", no_argument, 0, 'v'}



const std::vector<std::string> VALID_LONG_OPTION_STRINGS {{
    "128bit",
    "BMH",
    "PMH",
    "asymmetric",
    "asymmetric-all-pairs",
    "bagminhash",
    "bagminhash",
    "batch-size",
    "bbit-sigs",
    "bed",
    "bigwig",
    "binary",
    "binary-output",
    "bmh",
    "by-chrom",
    "cache",
    "cache-sketches",
    "cmp-outfile",
    "cmpout",
    "compute-edit-distance",
    "containment",
    "count-threshold",
    "countdict",
    "countmin-size",
    "countsketch-size",
    "distance",
    "distout",
    "doph",
    "downsample",
    "edit-distance",
    "edit-distance",
    "emit-binary",
    "enable-protein",
    "entmin",
    "exact-kmer-dist",
    "exact-kmer-dist",
    "fastcmp",
    "fastcmp-bytes",
    "fastcmp-nibbles",
    "fastcmp-shorts",
    "fastcmp-words",
    "ffile",
    "filterset",
    "full",
    "full-setsketch",
    "greedy",
    "help",
    "hp-compress",
    "intersection",
    "intersection-size",
    "kmer-length",
    "leafcutter",
    "long-kmers",
    "mash-distance",
    "maxcand",
    "multiset",
    "nLSH",
    "nlsh",
    "no-canon",
    "normalize-intervals",
    "one-perm",
    "oneperm",
    "oneperm-setsketch",
    "oph",
    "outfile",
    "outprefix",
    "pairlist",
    "parse-by-seq",
    "phylip",
    "pmh",
    "pminhash",
    "poisson-distance",
    "prefix",
    "prob",
    "probminhash",
    "probs",
    "protein",
    "protein14",
    "protein20",
    "protein6",
    "protein8",
    "qfile",
    "refine-exact",
    "regbytes",
    "regsize",
    "save-kmercounts",
    "save-kmers",
    "seed",
    "seq",
    "set",
    "set",
    "set",
    "setsketch-ab",
    "sig-ram-limit",
    "similarity-threshold",
    "sketch-size-l2",
    "sketchsize",
    "spacing",
    "square",
    "symmetric-containment",
    "threads",
    "threshold",
    "top-k",
    "topk",
    "union-size",
    "verbose",
    "window-size",
    "seqs-in-ram",
    //"machine-readable", // added this line so that I can use --machine-readable as a terminal flag to enable outputting distance computations to binary files
}};

// This function takes a vector of additional strings,
// since subcommands may have special flags unique to them.
// Add those strings to the second argument when validating those options.
// For instance, --presketched exists on cmp but not sketch, and must be permitted for it but not sketch.
static inline void validate_options(char **const opts, const std::vector<std::string>& extras={}) {
    for(char **ptr = opts; *ptr; ++ptr) {
        const int32_t len = std::strlen(*ptr);
        if(len > 2 && std::memcmp(*ptr, "--", 2) == 0) {
            const std::string flag((*ptr) + 2, (*ptr) + len);
            if(std::find(std::cbegin(extras), std::cend(extras), flag) != std::cend(extras)) {
                continue;
            }
            if(std::find(std::cbegin(VALID_LONG_OPTION_STRINGS), std::cend(VALID_LONG_OPTION_STRINGS), flag) == std::cend(VALID_LONG_OPTION_STRINGS)) {
                std::fprintf(stderr, "flag %s not found in expected set. See usage.\n", flag.data());
                throw std::runtime_error(std::string("Flag ") + flag + " not found");
            }
        }
    }
}



#define TOPK_FIELD case 'K': {ok = OutputKind::KNN_GRAPH; topk_threshold = std::atoi(optarg); break;}
#define SIMTHRESH_FIELD case 'T': {ok = OutputKind::NN_GRAPH_THRESHOLD; similarity_threshold = std::atof(optarg); break;}
#define GREEDY_FIELD case OPTARG_GREEDY: {\
    char *eptr;\
    ok = OutputKind::DEDUP; similarity_threshold = std::strtod(optarg, &eptr);\
    for(;*eptr;++eptr) {\
        if((*eptr | 32) == 'e') {exhaustive_dedup = true;}\
        if((*eptr | 32) == 'f') {fasta_dedup = true;}\
    }\
    break;\
}

#define CMPOUT_FIELD case OPTARG_CMPOUT: {cmpout = optarg; break;}
#define FASTCMP_FIELD case OPTARG_FASTCMP: {nbytes_for_fastdists = std::atof(optarg);\
            if(nbytes_for_fastdists != 8. && nbytes_for_fastdists != 4. && nbytes_for_fastdists != 2. && nbytes_for_fastdists != 1.){\
                std::fprintf(stderr, "--fastcmp must have 8, 4, 2, or 1 as the argument. These are the only register sizes supported.\n");\
                throw std::runtime_error("See usage for --fastcmp instructions.");\
            }\
            break;\
            }
#define PROT_FIELD case OPTARG_PROTEIN: {rht = bns::PROTEIN20; canon = false; std::fprintf(stderr, "Parsing 20-character amino acid sequences\n"); break;} \
    case OPTARG_PROTEIN6: {rht = bns::PROTEIN_6; canon = false; std::fprintf(stderr, "Parsing amino acid sequences with 6-letter compressed alphabet.\n"); break;}\
    case OPTARG_PROTEIN14: {rht = bns::PROTEIN14; canon = false; std::fprintf(stderr, "Parsing amino acid sequences with 14-letter compressed alphabet.\n"); break;}\
    case OPTARG_PROTEIN8: {rht = bns::PROTEIN8; canon = false; std::fprintf(stderr, "Using 3-bit protein encoding\n"); break;}
#define REFINEEXACT_FIELD case OPTARG_REFINEEXACT: {refine_exact = true; break;}

#define SHARED_FIELDS TOPK_FIELD SIMTHRESH_FIELD \
        CMPOUT_FIELD FASTCMP_FIELD PROT_FIELD REFINEEXACT_FIELD \
        GREEDY_FIELD \
        case 'E': sketch_space = SPACE_EDIT_DISTANCE; break;\
        case 'C': canon = false; break;\
        case 'p': nt = std::atoi(optarg); break;\
        case OPTARG_EXACT_KMER_DIST: {\
            exact_kmer_dist = true;\
            std::fprintf(stderr, "Exact kmer distance is set to %d -- %d/%d\n", exact_kmer_dist, c, option_index);\
            break;\
        }\
        case 'S': sketchsize = std::atoi(optarg); break;\
        case 'N': save_kmers = save_kmercounts = true; break;\
        case 's': save_kmers = true; break;\
        case OPTARG_SET: case 'H': res = FULL_MMER_SET; break;\
        case 'J': res = FULL_MMER_COUNTDICT; break;\
        case 'G': res = FULL_MMER_SEQUENCE; break;\
        case '2': use128 = true; break;\
        case 'm': count_threshold = std::atoi(optarg); break;\
        case 'F': ffile = optarg; break;\
        case 'Q': qfile = optarg; ok = PANEL; break;\
        case OPTARG_BED_NORMALIZE: normalize_bed = true; break;\
        case 'o': outfile = optarg; break;\
        case 'c': cssize = std::strtoull(optarg, nullptr, 10); break;\
        case 'k': k = std::atoi(optarg); break;\
        case 'w': w = std::atoi(optarg); break;\
        case 'W': cache = true; break;\
        case 'B': sketch_space = SPACE_MULTISET; res = FULL_SETSKETCH; break;\
        case 'P': sketch_space = SPACE_PSET; res = FULL_SETSKETCH; break;\
        case 'Z': res = ONE_PERM; break;\
        case 'v': ++verbosity; break;\
        case OPTARG_ISZ: measure = INTERSECTION; break;\
        case OPTARG_OUTPREF: {\
            outprefix = optarg; DBG_ONLY(std::fprintf(stderr, "outprefix: %s\n", outprefix.data());) break;\
        } \
        case OPTARG_HPCOMPRESS: {\
            hpcompress = true; break;\
        }\
        case OPTARG_DOWNSAMPLE_FRACTION: {\
            downsample_frac = std::atof(optarg); break;\
        }\
        case OPTARG_SPACING: spacing = optarg; canon = false; break;\
        case OPTARG_RANDOM_SEED: {seedseed = std::strtoull(optarg, 0, 10);} break;\
        case OPTARG_FILTERSET: fsarg = optarg; break;\
        case OPTARG_PARSEBYSEQ: parse_by_seq = true; break;\
        case OPTARG_CMP_BATCH_SIZE: batch_size = std::strtoull(optarg, 0, 10); break;\
        case OPTARG_NLSH: nLSH = std::atoi(optarg); break;\
        case OPTARG_ENTROPYMIN: entmin = true; break;\
        case 'L': {\
            const int ssl2 = std::atoi(optarg);\
            sketchsize = size_t(1) << ssl2;\
            if(ssl2 <= 0 || ssl2 >= 64) {\
                std::fprintf(stderr, "Error: ssl2 is out of bounds. ssl2: %d. sketchsize: %zu. Did you mean to specify sketchsize \n", ssl2, sketchsize);\
                std::exit(1);\
            }\
            std::fprintf(stderr, "Using log_2 sketchsize = %d, yielding sketchsize = %zu\n", ssl2, sketchsize);\
            break;\
        }\
        case OPTARG_SIGRAMLIMIT: {\
            MEMSIGTHRESH = std::strtoull(optarg, nullptr, 10);\
        } break;\
        case OPTARG_MAXCAND: {\
            maxcand_global = std::atoi(optarg);\
            if(maxcand_global < 0) {std::fprintf(stderr, "Warning: maxcand_global < 0. This defaults to heuristics for selecting the number of candidates. This may be in error.\n");}\
        } break;\
        case OPTARG_FULL_SETSKETCH: \
            res = FULL_SETSKETCH;sketch_space = SPACE_SET; break;\
        case OPTARG_SETSKETCH_AB: {\
            char *e;\
            compressed_a = std::strtold(optarg, &e);\
            compressed_b = std::strtold(e + 1, nullptr);\
            if(std::min(compressed_a, compressed_b) <= 0.) THROW_EXCEPTION(std::invalid_argument("Compressed A and B parameters must be greater than 0."));\
            std::fprintf(stderr, "setsketch ab: %Lg, %Lg\n", compressed_a, compressed_b);\
        } break;\
        case OPTARG_FASTCMPSHORTS: {\
            compressed_a = sketch::setsketch::ShortSetS::DEFAULT_A;\
            compressed_b = sketch::setsketch::ShortSetS::DEFAULT_B;\
            nbytes_for_fastdists = 2.;\
        } break;\
        case OPTARG_FASTCMPBYTES: {\
            compressed_a = sketch::setsketch::ByteSetS::DEFAULT_A;\
            compressed_b = sketch::setsketch::ByteSetS::DEFAULT_B;\
            nbytes_for_fastdists = 1.;\
        } break;\
        case OPTARG_FASTCMPNIBBLES: {\
            compressed_a = sketch::setsketch::NibbleSetS::DEFAULT_A;\
            compressed_b = sketch::setsketch::NibbleSetS::DEFAULT_B;\
            nbytes_for_fastdists = .5;\
        } break;\
        case OPTARG_FASTCMPWORDS: {\
            compressed_a = sketch::setsketch::UintSetS::DEFAULT_A;\
            compressed_b = sketch::setsketch::UintSetS::DEFAULT_B;\
            nbytes_for_fastdists = 4.;\
        } break;\
        case OPTARG_PAIRLIST: {\
            if(paths.size()) throw std::runtime_error("Expected empty paths list for pairlist command. Either provide pairlist or paths, not both.");\
            flat_hash_map<std::string, uint32_t> pathids;\
            std::ifstream ifs(optarg);\
            for(std::string line;std::getline(ifs, line);) {\
                auto spacepos = std::find_if(line.data(), line.data() + line.size(), [](auto x) {return std::isspace(x);});\
                if(*spacepos == '\0') throw std::runtime_error("Expected two paths separated by a space in pairlist.");\
                auto lhs = std::string(line.data(), spacepos - line.data());\
                auto rhs = std::string(spacepos + 1, line.data() + line.size() - spacepos - 1);\
                assert(!std::isspace(lhs.back()));\
                assert(!std::isspace(lhs.front()));\
                assert(!std::isspace(rhs.back()));\
                assert(!std::isspace(rhs.front()));\
                auto lit = pathids.find(lhs);\
                if(lit == pathids.end()) pathids.emplace(lhs, pathids.size()).first;\
                auto rit = pathids.find(rhs);\
                if(rit == pathids.end()) pathids.emplace(rhs, pathids.size()).first;\
                compareids.emplace_back(lit->second, rit->second);\
            }\
            for(auto &pair: pathids) paths.emplace_back(std::move(pair.first));\
            break;\
        }



static constexpr const char *siglen =
        sizeof(RegT) == 2 ? "16":
        sizeof(RegT) == 1 ? "8":
        sizeof(RegT) == 4 ? "32":
        sizeof(RegT) == 8 ? "64": "128";
#define SHARED_DOC_LINES \
        "General usage:\n"\
        "dashing2 <subcommand> <flags> [optional: positional arguments]\n\n"\
        "dashing2 takes all positional arguments as sketch objects.\nAlternatively, -F/--ffile takes inputs from a file in addition to positional arguments."\
        "Threading parallelism is set with -p/--threads.\n"\
        "\nThere are two main subcommands: sketch and cmp\n"\
        "\nSketch only generates sketches, while cmp performs comparisons across inputs, but sketches if necessary.\n"\
        "\n\nFile Format Options --\n"\
        "Dashing2 can sketch 4 kinds of files:"\
        "Fastq/Fasta, which has specific encoding options (default)\n"\
        "--bed to sketch BED files for interval sets\n"\
        "--bigwig to sketch BigWig files for coverage vectors\n"\
        "and --leafcutter to sketch LeafCutter splicing output\n"\
        "\n\nSequence Parsing Options --\n"\
        "If parsing fasta or fastq data, you can set several options:\n"\
        "  1. Alphabet [Default: DNA]. See 'Sequence Alphabet Options' below for more details.\n"\
        "  2. K-mer length (-k/--kmer-length) [Default: Maximum expressible in uint64_t (32 for DNA)]\n"\
        "  3. Window size (-w/--window-size) [Default: k-mer length]. Selects minimum-hash item from window length. Larger windows yields fewer minimizers.\n"\
        "  4. K-mer spacing. (--spacing) [Default: unspaced]. Allows for some positions to be ignored and others used. See below for more detail.\n\n"\
        "  5. K-mer encoding size. Defaults to uint64_t. 128-bit integers (up to 64bp DNA) can be enabled with --long-kmers.\n"\
        "     Dashing2 supports unbounded k-mer length by switching to rolling hashing if k is greater than the maximum possible length for the alphabet chosen.\n"\
        "  6. Disabling Canonicalization (--no-canon). By default, DNA alphabet k-mers are canonicalized to abstract strand. --no-canon causes Dashing2 to be strand-specific.\n"\
        "  7. Seed (--seed). To draw samples from a new hash function, you can provide a seed to the analysis.\n"\
        "  8. K-mer filtering, such as from a set of sequences, downsampling, or minimum count. See Detailed Filtering Options below for details.\n"\
        "Detailed Sequence Parsing Options\n"\
        "-k/--kmer-length: set k. Defaults to the largest k expressible directly in uint64_t.\n"\
        "If k is greater than this limit (31 for DNA, 14 for --protein, 22 for --protein8, 24 for --protein6), then rolling hashes will be generated instead of exact k-mer encodings.\n"\
        "-w/--window-size: set window size for winnowing; by default, all k-mers are used. If w > k, then only the minimum-hash k-mer in each window is processed\n"\
        "This can be useful to speed up sketching or to reduce the number of items sketched in exact mode.\n"\
        "--entmin: If -w/--window-size is enabled, this option weights the hash value by the entropy of the k-mer itself.\nThis is only valid for k-mers short enough to be encoded exactly in 64-bit or 128-bit integers, depending on if --long-kmers is enabled.\n"\
        "--spacing: Set a spacing scheme for spaced minimizers\n"\
        "This must have 1 less integer than the k-mer length.\n"\
        "e.g., -k5 --spacing 0,1,1,0 specifies a match length of 5 with a match pattern of `KK$K$KK`, where $ positions are ignored and `K` positions are kept.\n"\
        "This can also be run-length compressed in <space, num> format.\n"\
        "For example, --spacing 0,1x2,0 is equivalent to 0,1,1,0.\n"\
        "-2/--128bit/long-kmers: Use 128-bit k-mer hashes instead of 64-bit\n"\
        "Detailed Sequence Alphabet Options\n\n"\
        "Dashing2 sketches DNA by default. This can be changed with the following flags; this will disable canonicalization.\n"\
        "--enable-protein: Use 20 character amino acid alphabet.\n"\
        "--protein14: Use 14 character amino acid alphabet.\n"\
        "--protein6: Use 6 character amino acid alphabet.\n"\
        "--protein8: Use 8 character (3-bit) amino acid alphabet.\n"\
        "--no-canon: If DNA is being encoded, this disables canonicalization. By default, DNA sequence is canonicalized with its reverse-complement.\n"\
        "            Otherwise, this is ignored\n"\
        "--seed: Set a seed for k-mer hashing; If 0, this disables k-mer XORing and k-mers are encoded directly if a k-mer type can represent it.\n"\
        "        Otherwise, this changes the hash function applied to k-mers when generated sorted hash sets. This makes it easy to decode quickly, but we can still get good bottom-k estimates using these hashes\n"\
        "        the xor value for u64 kmers (unless --long-kmers is enabled) is the Wang 64-bit hash of the seed.\n"\
        "        u128 kmers (--long-kmers) have the same lower 64 bits, but the upper 64 bits are the Wang 64-bit hash of the u64 xor value.\n"\
        "Detailed Filtering Options\n\n"\
        "--downsample\t Downsample minimizers at fraction <arg> . Default is 1: IE, all minimizers pass.\n"\
        "-m/--threshold/--count-threshold <arg>: Set a count threshold for inclusion. If set to > 1, this will only sketch k-mers with count >= <arg>\n"\
        "If there are a set of common k-mers or artefactual sequence, you can specify --filterset to skip k-mers in this file when sketching other files.\n"\
        "By default, this converts it into a sorted hash set and skips k-mers which are found in the set.\n"\
        "`--filterset [path]` yields this.\n"\
        "\n\nInput File Options --\n"\
        "By default, dashing2 reads positional arguments and sketches them. You may want to use flags instructing it\n"\
        "to read from paths in <file>. Additionally, you can put multiple files separated by spaces into a single line "\
        "to place them all into a single sketch instead.\n"\
        "-F/--ffile: read paths from file in addition to positional arguments.\n"\
        "-Q/--qfile: read query paths from file; this is used for asymmetric queries (e.g., containment).\n"\
        "If multiple files are space-delimited in a single line, they will be sketched jointly.\n\n"\
        "\nSketch options\n"\
        "-S/--sketchsize: Set sketchsize (1024)\n"\
        "Warning: this is a significant change from Dashing1, where sketch size was required to be a power of 2, and so the sketch size would be 2**<arg>.\n"\
        "This restriction is no longer in place, and so any positive integer can be selected.\n"\
        "e.g. - instead of -S12 in Dashing, you would specify -S 4096\n"\
        "For convenience, we offer another argument to specify the size in log2.\n"\
        "-L/--sketch-size-l2: Set sketchsize to 2^<arg>. Must be > 0 and < 64\n"\
        "--cache/--cache-sketches: Save sketches to disk instead of re-computing each time.\n"\
        "\tDependent option:\n"\
        "\t                 --outprefix: specifies directory in which to save sketches instead of adjacent to the input files.\n"\
        "\t                 aliases: --prefix.\n"\
        "\t                 Note: You must have permission to write in the specified folder.\n"\
        "\nSketching Mode Options -- \n\n"\
        "Inputs can be summarized into several structures, and flags determine which is chosen.\n"\
        "1. SetSketch (one-permutation). Treats inputs as sets, ignoring multiplicities. The fastest option.\n"\
        "   This is faster at sketching, but has a small probability of failure which grows with sketch size. Big one-permutation sketches may perform poorly.\n"\
        "   (Default setting.)\n"\
        "2. FullSetSketch. Also treats inputs as sets but has better behavior for larger sketches and small sets.\n"\
        "   FullSetSketch is slower than one-permutation at sketching, and equally fast at comparisons than One-Permutation.\n"\
        "   --full/--full-setsketch to enable.\n"\
        "\n"\
        " We provide to weighted sketching algorithms -- WeightedSetSketch and DiscreteProbabilitySetSketch\n"\
        " Both require counting for sequence data, but do not for input methods with counting already performed, e.g. BigWig.\n"\
        " Full k-mer counting is enabled by default, but memory requirements can be fixed by using a count-min sketch during sketching.\n"\
        " Enabled by --countmin-size [number-registers], this allows for weighted sketching with fixed memory usage at the expense of some approximation.\n"\
        " This is only relevant to WeightedSetSketch and DiscreteProbabilitySetSketch.\n"\
        "3. WeightedSetSketch: Weighted Sets sketched by BagMinHash\n"\
        "   Multiset sketching via BagMinHash is an LSH for the weighted Jaccard similarity, which treats k-mer counts as weighted sets.\n"\
        "   -B/--multiset/--bagminhash to enable.\n"\
        "4. DiscreteProbabilitySetSketch: Discrete Probability Distributions using ProbMinHash\n"\
        "   ProbMinHash is an LSH for the probability Jaccard index, which normalizes observations by total counts, yielding a discrete probability distribution for each collection.\n"\
        "   ProbMinHash is most applicable for datasets where sampling fractions are important, such as expression or splicing counts.\n"\
        "   --prob/--pminhash\n"\
        "   ProbMinHash is 2-10x as fast as BagMinHash at sketching.\n"\
        " We also provide some exact comparison and sketching modes.\n"\
        " These include full k-mer sets (--set), a k-mer count dictionary (--countdict), or a sequence of minimizers (--seq)\n"\
        " 5. K-mer Sets.\n"\
        "   This generates a sorted hash set for k-mers in the data. If the parser is windowed (-w is fairly large), this could even be rather small.\n"\
        "   -H/--set to enable\n"\
        "   If an LSH table is generated, then weighted bottom-k hashes are used to build an LSH table\n"\
        " 6. Full k-mer countdict. \n"\
        "    This generates a sorted hash set for k-mers in the data, and additionally saves the associated counts for these k-mers.\n"\
        "    If an LSH table is generated, then weighted bottom-k hashes as in Cohen, E. \"Summarizing Data using Bottom-K Sketches\"\n"\
        "   -J/--countdict to enable \n"\
        " 7. Full k-mer (or minimizer) sequence. This faster than building the hash set, and can be used to build a minimizer index afterwards\n"\
        "          If you use --parse-by-seq with this and an output path is provided, then the stacked minimizer sequences will be written to it.\n"\
        "          The format is the similar to the standard stacked sketches, except that the cardinality fields instead represent minimizer sequence lengths (in 64-bit registers).\n"\
        "          Specifically, it consists of a header: [uint64_t nitems, uint32_t k, uint32_t w], followed by `nitems` [double], specifying sequence lengths of 64-bit registers\n"\
        "   -G/--seq to enable.\n"\
        "    Dependent option:\n"\
        "          --hp-compress:\n"\
        "              Minimizer sequence will be homopolymer-compressed before emission. \n"\
        "              This makes the sequences ignore the lengths of minimizer stretches.\n"\
        "\n\nOther Sketching Options -- \n"\
        "--parse-by-seq: Parse each sequence in each file as a separate entity. For workloads using edit distance, or for the --greedy mode, this will store all sequences in a temporary file in $TMPDIR.\n"\
        "                Previous versions of Dashing2 stored all sequences in memory with high memory usage for --parse-by-seq. This is reduced in v2.1.18.\n"\
        "                For faster use but more memory (restoring previous behavior), add --seqs-in-ram to avoid spilling to disk.\n"\
        "-s/--save-kmers: Save k-mers. This puts the k-mers saved into .kmer files to correspond with the minhash samples. \n"\
        "  If an output path is specified for dashing2 and --save-kmers is enabled, stacked k-mers will be written to <arg>.kmer64, and names will be written to <arg>.kmer.names.txt\n"\
        "  This has a 16-byte header containing a 32-bit integer describing the alphabet used, 32 bits describing sketch size, one 32-bit integer for k, and one 32-bit integer for window-length.\n"\
        "  This database can be used for dashing2 contain.\n"\
        "-N/--save-kmercounts: Save k-mer counts for sketches. This puts the k-mer counts saved into .kmercounts.f64 files to correspond with the k-mers.\n"\
        "-o/--outfile: sketches are stacked into a single file and written to <arg>\n"\
        "  This is the path for the stacked sketches; to set output location, use --cmpout instead. (This is the distance matrix betweek sketches).\n"\
        "  This can also reduce memory requirements, as the destination file is memory-mapped instead of held in memory.\n"\
        "\n\nComparison Options -- \n"\
        "We provide exhaustive (all-vs-all) comparisons, top-k nearest neighbor (KNN) graph generation and clustering. \n"\
        "  Within Exhaustive Comparisons, we have dense and sparse. The first 3 are dense:\n"\
        "    1. Upper Triangular PHYLIP (default)\n"\
        "    2. Square distance matrix (--asymmetric-all-pairs/--square)\n"\
        "    3. Rectangular distance matrix (--qfile/-Q)\n"\
        "  Instead of performing pairwise comparisons across one set, you can instead compare all in set X against all in set Y.\n"\
        "  We call this rectangular. When enabled with a path, entries on each line of the file at that path are treated as entities.\n"\
        "  positional arguments and -F paths are treated as a reference set;\n"\
        "  Paths provided in -Q/--qfile are are treated as a query set.\n"\
        "  Performing --asymmetric-all-pairs with the same input for -F and -Q should yield equivalent results.\n"\
        "  The output shape then has |F| rows and |Q| columns, (F, Q) in row-major format\n\n"\
        " We also support Sparse Exhaustive Comparisons, where the results are limited to more important entries.\n"\
        "  These include: \n"\
        "  1. Top-K, (--topk) only the top-k nearest neighbors are emitted for each item\n"\
        "  2. Thresholded, (--similarity-threshold), which emits entries only if the similarity exceeds the given threshold.\n"\
        "Both change the output format from a full matrix into compressed-sparse row (CSR) format listing only the filtered entries;\n"\
        "All of these are powered by the use of an LSH table built over the sketches, with the exception of exact mode (--countdict or --set), which use an LSH index built over their bottom-k hashes.\n"\
        "For details on LSH table parameters, see `LSH Options` below.\n"\
        "Top-K (K-Nearest-Neighbor) mode -- \n"\
        "--topk/--top-k <arg>\tMaximum number of nearest neighbors to list. If <arg> is greater than N - 1, pairwise distances are instead emitted.\n"\
        "\nThresholded Mode -- \n"\
        "--similarity-threshold <arg>\tMinimum fraction similarity for inclusion.\n\tIf this is enabled, only pairwise similarities over <arg> will be emitted.\n"\
        "\n\n"\
        "Greedy HIT Clustering Options --\n"\
        "In addition to exhaustive comparisons, we also perform greedy clustering using the CD High-Identity with Tolerance (CD-HIT) algorithm.\n"\
        "--greedy <float (0-1]> For greedy clustering by a given similarity threshold.\n"\
        "This uses an LSH index by default. \n"\
        "    To compare all points to all clusters, add E to the end of the flag. (e.g., '--greedy 0.8E')"\
        "    By default, this emits the names of the entities (sequence names, if --parse-by-seq, and filenames otherwise).\n"\
        "    You may want to emit fasta-formatted output. You can do this by adding F to the end of the --greedy argument.\n"\
        "    Example: '--greedy 0.8F' or '--greedy 0.8FE'.\n"\
        "    This is only allowed for --parse-by-seq.\n"\
        "  As this number approaches 1, the number and uniformity of clusters grows.\n"\
        "  For human-readable output, this emits one line per cluster listing its constituents, ordered by similarity\n"\
        "  For machine-readable output, this file consists of 2 64-bit integers (nclusters, nsets), followed by (nclusters + 1) 64-bit integers, followed by nsets 64-bit integers, identifying which sets belonged to which clusters.\n"\
        "  This is a vector in Compressed-Sparse notation.\n"\
        "  Python code for parsing the binary representation is available at https://github.com/dnbaker/dashing2/blob/main/python/parse.py.\n"\
        "\nDistance Register Size Options --\n"\
        "By default, we compare items with full hash function registers; to trade accuracy for speed, these sketches can be compressed before comparisons.\n"\
        "To truncate for faster comparisons, you can either select a register size to which to truncate to after generating full floating-point values, which will allow Dashing2 to use as much resolution as possible in a fixed number of bits.\n"\
        "On the other hand, this means that initial sketching needs more memory.\n"\
        "If you provide register values ahead of time, you can accumulate in smaller registers directly to reduce memory requirements.\n"\
        "--fastcmp/--regsize <arg>\tEnable faster comparisons using n-byte signatures rather than full registers. By default, these are logarithmically-compressed\n"\
        "  You can use this first approach with --fastcmp/--regsize. These are logarithmically compressed.\n"\
        "  For example, --fastcmp 1 uses byte-sized sketches, with a and b parameters inferred by the data to minimize information loss\n"\
        "  <arg> may be 8 (64 bits), 4 (32 bits), 2 (16 bits), or 1 (8 bits)\n"\
        "  Results may even be somewhat stabilized by the use of smaller registers.\n"\
        "\tDependent Options:\n"\
        "\t          --setsketch-ab <float1>,<float2>\n"\
        "\t            Example: --setsketch-ab 0.4,1.005\n"\
        "\t            If you specify a, b before using this method, then SetSketches of --fastcmp size bytes will be generated directly, rather than truncating after sketching all items.\n"\
        "\t            However, this is only supported for the SetSketch. (e.g., not for --bagminhash or --probminhash).\n"\
        "\t            This can allow you to scale to larger collections.\n"\
        "\t            We also provide several pre-set parameter values.\n"\
        "\t          --fastcmp-bytes sets a and b to 20 and 1.2, and sets --fastcmp to 1\n"\
        "\t          --fastcmp-shorts sets a and b to .06 and 1.0005, and sets --fastcmp to 2.\n"\
        "\t          --fastcmp-words sets a and b to 19.77 and 1.0000000109723500835 and sets --fastcmp to 4.\n"\
        /*"\t          --fastcmp-nibbles sets a and b to .0005 and 2.71828, and sets --fastcmp to .5\n"*/\
        "\n"\
        "If you instead want to truncate to the bottom-b bits of the signature --\n"\
        "\t          --bbit-sigs: truncate to bottom-<arg> bytes of signatures instead of logarithmically-compressed.\n"\
        "The runtime is effectively equivalent to the setsketch.\n"\
        "\n\nComparison Function Options --\n"\
        "The default comparison emitted is similarity. For MinHash/HLL/SetSketch sketches, this is the fraction of shared registers.\n"\
        "This can be changed to a distance (--mash-distance) for k-mer similarity, where it can be used for hierarchical clustering.\n"\
        "--mash-distance/--poisson-distance/--distance\t Emit distances, as estimated by the Poisson model for k-mer distances.\n"\
        "--symmetric-containment\t Use symmetric containment as the distance. e.g., (|A & B| / min(|A|, |B|))\n"\
        "--containment\t Use containment as the distance. e.g., (|A & B| / |A|). This is asymmetric, so you must consider that when deciding the output shape.\n"\
        "--intersection-size/--intersection\t Emit the cardinality of the intersection between entities. IE, the number of k-mers shared between the two.\n"\
        "--union-size\t Emit the cardinality of the union between entities. IE, the number of k-mers in the union of the two.\n"\
        "--compute-edit-distance\t For edit distance, perform actual edit distance calculations rather than returning the distance in LSH space.\n"\
        "                       \t This means that the LSH index eliminates the quadratic barrier in candidate generation, but they are refined using actual edit distance.\n"\
        "\n"\
        "Distance Output Options --\n"\
        "--binary-output\tEmit binary output rather than human-readable.\n\n"\
        "\t For symmetric pairwise, this emits condensed distance matrix in f32\n"\
        "\t For -Q/--qfile usage, this emits a full rectangular distance matrix in of shape (|F|, |Q|)\n"\
        "\t For asymmetric pairwise, this emits a full distance matrix in f32 in row-major storage.\n"\
        "\t For asymmetric pairwise, this emits a flat distance matrix in f32\n"\
        "\t For top-k filtered, this emits a matrix of min(k, |N|) x |N| of IDs and distances\n"\
        "In `dashing2 sketch`, distances are not automatically computed; If set, they are written to <arg>\n"\
        "In `dashing2 cmp`, this defaults to stdout.\n"\
        "--cmpout/--distout/--cmp-outfile\tCompute distances and emit them to <arg>.\n"\
        "\t For similarity-thresholded distances, this emits a compressed-sparse row (CSR) formatted matrix with 64-bit indptr, 32-bit indices, and 32-bit floats for distances\n"\
        "\n\nLSH Options --\n"\
        "There are a variety of heuristics in the LSH tables; however, the most important besides sketch size is the number of hash tables used.\n"\
        "--nLSH <int=2>\t\n"\
        "   Aliases: --nlsh\n"\
        "This sets the number of LSH tables. The first 3 tables use register grouping sizes of powers of 2, and subsequent tables use 2 times the index.\n"\
        "If 2 is used (default), these will be of sizes (1, 2), but 4 yields (1, 2, 4, 6) and 5 yields (1, 2, 4, 6, 8) registers per composite key.\n"\
        "Increase this number to pay more memory/time for higher accuracy.\n"\
        "Decrease this number for higher speed and lower accuracy.\n"\
        "This is ignored for exact sketching (--countdict or --set), where a single permutation is generated and a single hash table is used.\n"\
        "--maxcand <int>\t Set the maximum number of candidates to fetch from the LSH index before evaluating distances against them.\n"\
        "                  This is always used in --greedy mode.\n"\
        "                  By default, this number is heuristically selected by the number of items in the index.\n"\
        "                  If N <= 100000, uses max(N / 50, max(cbrt(N), 3)). Otherwise, uses (log(N)^3).\n"\
        "                  Note: this is a maximum; if fewer items have matches, fewer will be reported.\n"\
        "                  This option is ignored in --topk mode, as the number of samples is ceil(3.5 * <topk>).\n"\
        "                  If set in --similarity-threshold mode, the number of items compared will be truncated to <maxcand> even if further samples are above the similarity threshold.\n"\
        "                  This can prevent quadratic complexity for the (rare) case that all items are within threshold distaance of each other.\n"\


extern size_t MEMSIGTHRESH;

}

#endif /* DASHING2_OPTIONS_H__ */
