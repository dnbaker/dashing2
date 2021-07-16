#pragma once
#ifndef DASHING2_OPTIONS_H__
#define DASHING2_OPTIONS_H__
#include <getopt.h>

namespace dashing2 {
#define LO_ARG(LONG, SHORT) {LONG, required_argument, 0, SHORT},
#define LO_NO(LONG, SHORT) {LONG, no_argument, 0, SHORT},
#define LO_FLAG(LONG, SHORT, VAR, VAL) {LONG, no_argument, (int *)&VAR, VAL},

using option_struct = struct option;

enum OptArg{
    OPTARG1 = 1000,
    OPTARG_BED,
    OPTARG_BIGWIG,
    OPTARG_LEAFCUTTER,
    OPTARG_OUTPREF,
    OPTARG_CMPOUT,
    OPTARG_BINARY_OUTPUT,
    OPTARG_FASTCMP,
    OPTARG_REGBYTES,
    OPTARG_BBIT_SIGS,
    OPTARG_EXACT_KMER_DIST,
    OPTARG_ISZ,
    OPTARG_BED_NORMALIZE,
    OPTARG_PROTEIN,
    OPTARG_HPCOMPRESS,
    OPTARG_DOWNSAMPLE_FRACTION,
    OPTARG_REFINEEXACT,
    OPTARG_MASHDIST,
    OPTARG_SYMCONTAIN,
    OPTARG_CONTAIN,
    OPTARG_ASYMMETRIC_ALLPAIRS,
    OPTARG_SET,
    OPTARG_DUMMY
};

#define SHARED_OPTS \
    LO_ARG("save-kmercounts", 'N')\
    LO_ARG("ffile", 'F')\
    LO_ARG("qfile", 'Q')\
    LO_ARG("threads", 'p')\
    LO_ARG("sketchsize", 'S')\
    LO_ARG("cmpout", OPTARG_CMPOUT)\
    LO_ARG("distout", OPTARG_CMPOUT)\
    LO_ARG("cmp-outfile", OPTARG_CMPOUT)\
    LO_ARG("outprefix", OPTARG_OUTPREF)\
    LO_ARG("topk", 'K')\
    LO_ARG("similarity-threshold", 'T')\
    LO_ARG("fastcmp", OPTARG_FASTCMP)\
    LO_ARG("countsketch-size", 'c')\
    LO_ARG("window-size", 'w')\
    LO_ARG("kmer-length", 'k')\
    LO_ARG("outfile", 'o')\
    LO_ARG("count-threshold", 'm')\
    LO_ARG("threshold", 'm')\
    LO_FLAG("binary-output", OPTARG_BINARY_OUTPUT, of, OutputFormat::MACHINE_READABLE) \
    LO_FLAG("parse-by-seq", OPTARG_PARSEBYSEQ, parse_by_seq, true)\
    LO_FLAG("bagminhash", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("prob", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("bed", OPTARG_BED, dt, DataType::BED)\
    LO_FLAG("bigwig", OPTARG_BIGWIG, dt, DataType::BIGWIG)\
    LO_FLAG("leafcutter", OPTARG_LEAFCUTTER, dt, DataType::LEAFCUTTER)\
    /*LO_FLAG("edit-distance", 'E', sketch_space, SPACE_EDIT_DISTANCE)*/\
    /*LO_FLAG("set", 'H', res, FULL_MMER_SET)*/\
    /*LO_FLAG("exact-kmer-dist", OPTARG_EXACT_KMER_DIST, exact_kmer_dist, true)*/\
    LO_FLAG("bbit-sigs", OPTARG_BBIT_SIGS, truncate_mode, 1)\
    LO_FLAG("intersection", OPTARG_ISZ, measure, INTERSECTION)\
    LO_FLAG("mash-distance", OPTARG_MASHDIST, measure, POISSON_LLR)\
    LO_FLAG("distance", OPTARG_MASHDIST, measure, POISSON_LLR)\
    LO_FLAG("symmetric-containment", OPTARG_SYMCONTAIN, measure, SYMMETRIC_CONTAINMENT)\
    LO_FLAG("containment", OPTARG_CONTAIN, measure, SYMMETRIC_CONTAINMENT)\
    LO_FLAG("poisson-distance", OPTARG_MASHDIST, measure, POISSON_LLR)\
    LO_FLAG("compute-edit-distance", OPTARG_MASHDIST, measure, M_EDIT_DISTANCE)\
    LO_FLAG("multiset", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("countdict", 'J', res, FULL_MMER_COUNTDICT)\
    LO_FLAG("seq", 'G', res, FULL_MMER_SEQUENCE)\
    LO_FLAG("128bit", '2', use128, true)\
    LO_FLAG("long-kmers", '2', use128, true)\
    LO_FLAG("save-kmers", 's', save_kmers, true)\
    LO_FLAG("asymmetric-all-pairs", OPTARG_ASYMMETRIC_ALLPAIRS, ok, OutputKind::ASYMMETRIC_ALL_PAIRS)\
    LO_ARG("regbytes", OPTARG_REGBYTES)\
    /*LO_ARG("set", 'H')*/\
    {"hp-compress", no_argument, 0, OPTARG_HPCOMPRESS},\
    {"refine-exact", no_argument, 0, OPTARG_REFINEEXACT},\
    {"edit-distance", no_argument, 0, 'E'},\
    {"full-setsketch", no_argument, 0, 'Z'},\
    {"normalize-intervals", no_argument, 0, OPTARG_BED_NORMALIZE},\
    {"enable-protein", no_argument, 0, OPTARG_PROTEIN},\
    {"downsample", required_argument, 0, OPTARG_DOWNSAMPLE_FRACTION},\
    {"cache", no_argument, 0, 'W'},\
    {"no-canon", no_argument, 0, 'C'},\
    {"set", no_argument, 0, OPTARG_SET},\
    {"exact-kmer-dist", no_argument, 0, OPTARG_EXACT_KMER_DIST}


#define TOPK_FIELD case 'K': {ok = OutputKind::KNN_GRAPH; topk_threshold = std::atoi(optarg); break;}
#define SIMTHRESH_FIELD case 'T': {ok = OutputKind::NN_GRAPH_THRESHOLD; similarity_threshold = std::atof(optarg); break;}
#define CMPOUT_FIELD case OPTARG_CMPOUT: {cmpout = optarg; break;}
#define FASTCMP_FIELD case OPTARG_FASTCMP: {nbytes_for_fastdists = std::atof(optarg); break;}
#define PROT_FIELD case OPTARG_PROTEIN: {rht = bns::PROTEIN; break;}
#define REFINEEXACT_FIELD case OPTARG_REFINEEXACT: {refine_exact = true; break;}

#define SHARED_FIELDS TOPK_FIELD SIMTHRESH_FIELD CMPOUT_FIELD FASTCMP_FIELD PROT_FIELD REFINEEXACT_FIELD \
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
        case OPTARG_SET: case 'H': res = FULL_MMER_SET; std::fprintf(stderr, "Result should now be FULL_MMER_SET %s\n", to_string(res).data()); break;\
        case 'J': res = FULL_MMER_COUNTDICT; break;\
        case 'G': res = FULL_MMER_SEQUENCE; break;\
        case '2': use128 = true; break;\
        case 'm': count_threshold = std::atof(optarg); break;\
        case 'F': ffile = optarg; break;\
        case 'Q': qfile = optarg; ok = PANEL; break;\
        case OPTARG_BED_NORMALIZE: normalize_bed = true; break;\
        case OPTARG_REGBYTES: nbytes_for_fastdists = std::atof(optarg); break;\
        case 'o': outfile = optarg; break;\
        case 'c': cssize = std::strtoull(optarg, nullptr, 10); break;\
        case 'k': k = std::atoi(optarg); break;\
        case 'w': w = std::atoi(optarg); break;\
        case 'W': cache = true; break;\
        case 'B': sketch_space = SPACE_MULTISET; res = FULL_SETSKETCH; break;\
        case 'P': sketch_space = SPACE_PSET; res = FULL_SETSKETCH; break;\
        case 'Z': res = FULL_SETSKETCH; break;\
        case OPTARG_ISZ: measure = INTERSECTION; break;\
        case OPTARG_OUTPREF: {\
            outprefix = optarg; break;\
        } \
        case OPTARG_HPCOMPRESS: {\
            hpcompress = true; break;\
        }\
        case OPTARG_DOWNSAMPLE_FRACTION: {\
            downsample_frac = std::atof(optarg); break;\
        }


static constexpr const char *siglen =
        sizeof(RegT) == 2 ? "16":
        sizeof(RegT) == 1 ? "8":
        sizeof(RegT) == 4 ? "32":
        sizeof(RegT) == 8 ? "64": "128";
#define SHARED_DOC_LINES \
        "\n\nDistance Options (shared)\n"\
        "--cmpout/--distout/--cmp-outfile\tCompute distances and emit them to [arg].\n"\
        "--similarity-threshold [arg]\tMinimum fraction similarity for inclusion.\n\tIf this is enabled, only pairwise similarities over [arg] will be emitted.\n"\
        "This changes the output format from a full matrix into compressed-sparse row (CSR) notation.\n"\
        "--topk [arg]\tMaximum number of nearest neighbors to list. If [arg] is greater than N - 1, pairwise distances are instead emitted.\n"\
        "--binary-output\tEmit binary output rather than human-readable.\n\t For symmetric pairwise, this emits condensed distance matrix in f32\n"\
        "-Q/--qfile\tThis causes rectangular output between -F filenames ane -Q filenames. This option is listed twice, as it changes the format of the data emitted.\n"\
        "Use this if you want to compare a set of queries to a set of references rather than complete all-pairs. Note: -F must be provided, or reference files should be added as positional arguments\n"\
        "\t For asymmetric pairwise, this emits a flat distance matrix in f32\n"\
        "\t For top-k filtered, this emits a matrix of min(k, |N|) x |N| of IDs and distances\n"\
        "\t For similarity-thresholded distances, this emits a compressed-sparse row (CSR) formatted matrix\n"\
        "--fastcmp [arg]\tEnable faster comparisons using n-byte signatures rather than full registers. By default, these are set-sketch compressed\n"\
        "For example, --fastcmp 1 uses byte-sized sketches, with a and b parameters inferred by the data to minimize information loss\n"\
        "\t If --bbit-sigs is enabled, this random signatures truncated to [arg] bytes will be replaced.\n"\
        "\t The tradeoff is that you may get better accuracy in set space comparisons at the expense of information regarding the sizes of the sets\n"\
        "--exact-kmer-dist\tThis uses exact k-mer distances instead of approximate methods. In edit distance space, this means calculating exact edit distance.\n"\
        "For minhash sketches (setsketch, probminhash, and bagminhash), this uses full registers instead of compressed.\n"\
        "--refine-exact\tThis causes the candidate KNN graph to be refined to a final KNN graph using full distances.\tIf using sketches, then full hash registers are used.\nOtherwise, exact k-mer comparison functions are used.\n"\
        "--downsample\t Downsample minimizers at fraction [arg] . Default is 1: IE, all minimizers pass.\n"\
        "\n\nDistance specifications\n"\
        "The default value emitted is similarity. For MinHash/HLL/SetSketch sketches, this is the fraction of shared registers.\n"\
        "This can be changed to a distance (--mash-distance) for k-mer similarity, where it can be used for hierarchical clustering.\n"\
        "--mash-distance/--poisson-distance\t Emit distances, as estimated by the Poisson model for k-mer distances.\n"\
        "--symmetric-containment\t Use symmetric containment as the distance. e.g., (|A & B| / min(|A|, |B|))\n"\
        "--containment\t Use containment as the distance. e.g., (|A & B| / |A|). This is asymmetric, so you must consider that when deciding the output shape.\n"\
        "--compute-edit-distance\t For edit distance, perform actual edit distance calculations rather than returning the distance in LSH space.\n"\
        "                       \t This means that the LSH index eliminates the quadratic barrier in candidate generation, but they are refined using actual edit distance.\n"



}

#endif /* DASHING2_OPTIONS_H__ */
