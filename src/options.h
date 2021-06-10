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
    LO_ARG("cmpout", OPTARG_CMPOUT)\
    LO_ARG("distout", OPTARG_CMPOUT)\
    LO_ARG("cmp-outfile", OPTARG_CMPOUT)\
    LO_ARG("topk", 'K')\
    LO_ARG("similarity-threshold", 'T')\
    LO_ARG("fastcmp", OPTARG_FASTCMP)\
    LO_FLAG("bbit-sigs", OPTARG_BBIT_SIGS, truncate_mode, 1)\
    LO_FLAG("binary-output", OPTARG_BINARY_OUTPUT, of, OutputFormat::MACHINE_READABLE) \
    LO_FLAG("parse-by-seq", OPTARG_PARSEBYSEQ, parse_by_seq, true)\
    LO_FLAG("bagminhash", OPTARG_DUMMY, s, SPACE_MULTISET)\
    LO_FLAG("prob", 'P', s, SPACE_PSET)\
    LO_FLAG("edit-distance", 'E', s, SPACE_EDIT_DISTANCE)\
    LO_FLAG("set", 'H', res, FULL_MMER_SET)\
    LO_FLAG("exact-kmer-dist", OPTARG_EXACT_KMER_DIST, exact_kmer_dist, true)\
    {"full-setsketch", no_argument, 0, 'Z'},\
    {"normalize-intervals", no_argument, 0, OPTARG_BED_NORMALIZE},\
    {"enable-protein", no_argument, 0, OPTARG_PROTEIN},\
    LO_FLAG("intersection", OPTARG_ISZ, measure, INTERSECTION)


#define TOPK_FIELD case 'K': {ok = OutputKind::KNN_GRAPH; topk_threshold = std::atoi(optarg); break;}
#define SIMTHRESH_FIELD case 'T': {ok = OutputKind::NN_GRAPH_THRESHOLD; similarity_threshold = std::atof(optarg); break;}
#define CMPOUT_FIELD case OPTARG_CMPOUT: {cmpout = optarg; break;}
#define FASTCMP_FIELD case OPTARG_FASTCMP: {nbytes_for_fastdists = std::atoi(optarg); break;}
#define PROT_FIELD case OPTARG_PROTEIN: {rht = bns::PROTEIN; std::fprintf(stderr, "Enabled protein\n"); break;}

#define SHARED_FIELDS TOPK_FIELD SIMTHRESH_FIELD CMPOUT_FIELD FASTCMP_FIELD PROT_FIELD


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
