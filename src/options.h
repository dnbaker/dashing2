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
    LO_FLAG("canon", 'C', canon, true)\
    LO_FLAG("cache", 'W', cache, true)\
    LO_FLAG("binary-output", OPTARG_BINARY_OUTPUT, of, OutputFormat::MACHINE_READABLE) \
    LO_FLAG("parse-by-seq", OPTARG_PARSEBYSEQ, parse_by_seq, true)\
    LO_FLAG("bagminhash", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("prob", 'P', sketch_space, SPACE_PSET)\
    LO_FLAG("bed", OPTARG_BED, dt, DataType::BED)\
    LO_FLAG("bigwig", OPTARG_BIGWIG, dt, DataType::BIGWIG)\
    LO_FLAG("leafcutter", OPTARG_LEAFCUTTER, dt, DataType::LEAFCUTTER)\
    /*LO_FLAG("edit-distance", 'E', sketch_space, SPACE_EDIT_DISTANCE)*/\
    LO_FLAG("set", 'H', res, FULL_MMER_SET)\
    LO_FLAG("exact-kmer-dist", OPTARG_EXACT_KMER_DIST, exact_kmer_dist, true)\
    LO_FLAG("bbits-sigs", OPTARG_BBIT_SIGS, truncate_mode, 1)\
    LO_FLAG("intersection", OPTARG_ISZ, measure, INTERSECTION)\
    LO_FLAG("multiset", OPTARG_DUMMY, sketch_space, SPACE_MULTISET)\
    LO_FLAG("countdict", 'J', res, FULL_MMER_COUNTDICT)\
    LO_FLAG("seq", 'G', res, FULL_MMER_SEQUENCE)\
    LO_FLAG("128bit", '2', use128, true)\
    LO_FLAG("long-kmers", '2', use128, true)\
    LO_FLAG("save-kmers", 's', save_kmers, true)\
    LO_ARG("regbytes", OPTARG_REGBYTES)\
    {"edit-distance", no_argument, 0, 'E'},\
    {"full-setsketch", no_argument, 0, 'Z'},\
    {"normalize-intervals", no_argument, 0, OPTARG_BED_NORMALIZE},\
    {"enable-protein", no_argument, 0, OPTARG_PROTEIN},\



#define TOPK_FIELD case 'K': {ok = OutputKind::KNN_GRAPH; topk_threshold = std::atoi(optarg); break;}
#define SIMTHRESH_FIELD case 'T': {ok = OutputKind::NN_GRAPH_THRESHOLD; similarity_threshold = std::atof(optarg); break;}
#define CMPOUT_FIELD case OPTARG_CMPOUT: {cmpout = optarg; break;}
#define FASTCMP_FIELD case OPTARG_FASTCMP: {nbytes_for_fastdists = std::atoi(optarg); break;}
#define PROT_FIELD case OPTARG_PROTEIN: {rht = bns::PROTEIN; break;}

#define SHARED_FIELDS TOPK_FIELD SIMTHRESH_FIELD CMPOUT_FIELD FASTCMP_FIELD PROT_FIELD \
        case 'E': sketch_space = SPACE_EDIT_DISTANCE; break;\
        case 'C': canon = true; break;\
        case 'p': nt = std::atoi(optarg); break;\
        case 'S': sketchsize = std::atoi(optarg); break;\
        case 'N': save_kmers = save_kmercounts = true; break;\
        case 's': save_kmers = true; break;\
        case 'H': res = FULL_MMER_SET; break;\
        case 'J': res = FULL_MMER_COUNTDICT; break;\
        case 'G': res = FULL_MMER_SEQUENCE; break;\
        case '2': use128 = true; break;\
        case 'm': count_threshold = std::atof(optarg); break;\
        case 'F': ffile = optarg; break;\
        case 'Q': qfile = optarg; break;\
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
        "\t For asymmetric pairwise, this emits a flat distance matrix in f32\n"\
        "\t For top-k filtered, this emits a matrix of min(k, |N|) x |N| of IDs and distances\n"\
        "\t For similarity-thresholded distances, this emits a compressed-sparse row (CSR) formatted matrix\n"\
        "--fastcmp [arg]\tEnable faster comparisons using n-byte signatures rather than full registers. By default, these are set-sketch compressed\n"\
        "For example, --fastcmp 1 uses byte-sized sketches, with a and b parameters inferred by the data to minimize information loss\n"\
        "\t If --bbit-sigs is enabled, this random signatures truncated to [arg] bytes will be replaced.\n"\
        "\t The tradeoff is that you may get better accuracy in set space comparisons at the expense of information regarding the sizes of the sets\n"


}

#endif /* DASHING2_OPTIONS_H__ */
