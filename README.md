## Introduction

Dashing2 is the second version of the Dashing sequence sketching and comparison system.

There have been several major changes, but you can still can get a quick start to compare a group of sequence collections [here](#quickstart).

[Installation instructions](#installation) can be found below.

### New Features

tl;dr --
1. Faster and more accurate distances.
2. Weight-informed sketching.
3. Near-linear K-nearest neighbor and thresholded neighbor graphs via LSH tables.
4. Minimizer sequence transduction.
5. Faster/easier installation than Dashing1
  1. Compilation in 45s in parallel, 2 min serial

Input Formats -- See [inputs](#inputs) below for details.
 1. Fastq/Fasta -- enhanced alphabet support
   1. Default -- ACGT
   2. Protein encoding - 20 characters (`--protein`)
   3. Reduced alphabets
     1. Protein - 14, 8, and 6-character alphabets for long-distance homology.
   4. Optional -- generating 128-bit k-mers (`--long-kmers`)
   5. All k-mer parsing can be winnowed by setting `--window-size` to be > k
   6. Seeds can be spaced by providing a `--spacing` option, which provides the number of ignored characters in between used characters for seeds.
   7. Exact multiset comparisons (`--countdict`)
   8. Weight-aware sketching -- multiset and probability distribution
     1. We use BagMinHash for weighted sets (`--bagminhash` or `--multiset`), and ProbMinHash for discrete probability distributions (`--prob`)
   9. Minimizer sequence transduction
     1. By enabling `--seq`, a sequence of minimizer values are emitted as a string.
     2. This can be used for simple minimizer generation, or these minimizer sequences's edit distances can be compared in downstream analysis.
 2. Splicing data
   1. LeafCutter splicing output files `--leafcutter`
 3. Interval Sets
   1. BigWig files (`--bigwig`) [Weighted Interval Set]
   2. BED files (`--bed`)
 4. Binary dumps
   1. Binary files made from 64-bit identifiers
   2. Optional: paired weights for each feature
   3. Can be sketched as sets, weighted sets, and discrete probability distributions.


There's also a two-step process; we can use the LSH index to generate probable candidates, but switch to exact distance comparisons after refinement.
This can be enabled with the `--refine-exact` flag when using filtered exact-set (`--set`) or exact-multiset (`--countdict`).

For sketching algorithms that are not OrderMinHash (whether the default SetSketch, BagMinHash or ProbMinHash), this means that candidates are generated with fast, compressed sketches,
but final results are produced using the full hash values.

If `--edit-distance` has been enabled, `--refine-exact` and `--compute-edit-distance` cause candidates to be generated by LSH table querying, but whose final distances are computed by exact edit ditsance.
This allows one to prune the edit distances with LSH but still get exact results out.

## Outputs

Outputs
   1. All-pairs symmetric (default), resulting in a compressed distance matrix of size N-choose-2 = `(N (N - 1)) / 2`.
   2. Rectangular comparisons (`--qfile/-Q` against `--ffile/-F`) by providing query/reference sets
   3. Filtered output
      1. Jaccard thresholded-results (`--similarity-threshold`) and top-k results (`--topk <k>`)
        2. Dashing2 builds and queries an LSH index to avoid all-pairs comparisons when generating these neighbor-graphs.
        3. This affords us near-linear time- and space- comparisons, which we can use for clustering, indexing, and summarization.
   4. Asymmetric all-pairs (`--asymmetric-all-pairs`), which performs c(A, B) for all-pairs, resulting in an NxN matrix.

All of these can be emitted in binary format `--binary-output` to avoid parsing/print formatting costs.

1. The matrices are emitted in flat, row-major format in float32.
   1. If human-readable, the all-pairs symmetric is emitted as PHYLIP.
   2. All-pairs asymmetric is emitted as a flat tsv.
   3. Panel (rectangular) results are emitted as a flat tsv.
2. Thresholded results are emitted in CSR-format; see `dashing2 cmp --help` for more details.
   1. Human-readable thresholded results are emitted as a tsv, with potentially varying numbers of items per line.
   2. Individual entries consist of the distance/similarity and the corresponding entity.


## Inputs

Dashing2 sketches input datasets, and then compares resulting sketches.

The formats supported are:

• FastA and FastQ
  - Uncompressed or gzip-compressed is transparently processed
  - Note: this supports both nucleotide and protein sketching.
  - Canonicalization is off by default.
  -       For k <= 64, DNA/RNA uses exact encoding, and for higher k values, a rolling hash is used.
  -       If protein is enabled (--enable-protein), then canonicalization does not apply and is not performed.
  - There is currently no quality filtering; to filter by quality scores, mask bases below desired quality as Ns.
• BigWig
  - Is not stratified by chromosome
• BED files
  -  By default, BED files are sketched as sets of reference base/contig pairs.
  -  These rows can be normalized (--normalize-intervals) to treat each *interval* as having the same weight,
     or can be treated as-is to treat each contig/reference base pair as an item to sketch.
  - Note: BED7 files have strand information, but Dashing2 does not attempt to stratify by strand
  -     Convert these files into artificial BED files with the strand information included to work around this for now.
• LeafCutter output files for splicing datasets

# Sketching Formats
## Sketching Algorithms
By default, we do set sketching (One-Permutation MinHash or SetSketch Minhash), but we can do multiset (BagMinHash) or discrete prob set (ProbMinHash) sketching.
BagMinhash is probably most appropriate for genome assemblies (--multiset), and ProbMinhash (--prob) is most appropriate for splicing and expression datasets, due to normalization.

We also have untested support of edit distance LSH. It seems to work, but hasn't been carefully vetted.

## Counting
Counting is necessary for both BagMinHash and ProbMinHash. Default behavior, which works with sets and discards quantities, does not require it.

For counting, we either use a (single-row) CountSketch to approximate weighted set comparisons or we compute exact counts using a hash table.
Counting is exact by default.

To enable count-sketch approximated counting, use the flag --countsketch-size [size]. The larger this parameter, the closer to exact the weighted sketching will be.

## Exact sketching
In addition to creating sketches for m-mer sets, Dashing2 can perform full m-mer sets and m-mer count dictionaries.
This can be enabled with --set or --countdict, respectively. This will be slower, but exact.


#### BED sketching
##### Normalization
Default normalization treats each base pair in an interval as equally weighted. For SPACE\_SET, this doesn't matter, but for weighted sketching, it does.

That is, an interval from 0:100-200 will have twice the weight as 0:100-150. If, instead, you want to treat each interval as weighing the same as another interval, enable --normalize-intervals.


#### LeafCutter sketching
Similar to BED, except that the normalize flag causes the splicing event weight to be the fraction of reads supporting the junction rather than the absolute number  of reads.

IE, 3/5 would have weight 0.6 if --normalize-intervals is enabled and weight 3 otherwise.


Sketches can be:
    SetSketch (One-permutation minhash or SetSketch)
    BagMinHash
    ProbMinHash
    OrderMinHash

However, OrderMinHash is only available for sequences.

When handling multiplicities, you can use exact counting, which may be slower, or you can approximate the count vectors with a single-row count-sketch
by setting opts.cssize\_ > 0.

ProbMinHash is usually significantly (2-20+x) faster than BagMinHash, although multiset jaccard may be more appropriate for some problems.
For instance, expression and chromatic accessibility might be better considered discrete probability distributions and therefore fit ProbMinHash, whereas genomic sequences better match the multiset concept and benefit from BagMinHash, which is a MinHash algorithm for weighted sets.


#### QuickStart

We expect most usage to involve the `sketch` subcommand; to sketch and perform distances, add `--cmpout <file>`, where <file> is the destination file or '-' to represent stdout.

**Use 1 -- Sketch \+ All-pairs Compare Sequence Collections**

```
dashing2 sketch [options] --cmpout <outfile> genome1.fa genome2.fa <...>
# Alternate, if filenames are in F.txt
# dashing2 sketch [options] --cmpout <outfile> -F F.txt
```

Full usage is found via `dashing2 --help` and `dashing2 <subcommand> --help`, where  <subcommand> is one of the dashing2 subcommands.

`dashing2 sketch` performs sketching/summarization of a set of input files, or sequence-by-sequence processing of one or more sequence files.
It also optionally performs comparisons and emits results to `--cmpout`.

Adding `--cache` causes Dashing2 to cache sketches to disk adjacent to the input files;
     this location can be changed with `--outprefix`.


We support a variety of alphabets -- DNA, Protein, and reduced amino acid alphabets for long-range homology (--protein14, --protein8, --protein6).

Alternate file-types supported include BigWigs (`--bigwig`), BED (`--bed`), and LeafCutter outputs (`--leafcutter`).


Clustering: scipy.cluster.hierarchy and fastcluster.hierarchy yield fast, concise clusterings, if distances are emitted.

To perform query-set vs reference-set comparison, see `-Q/--qfile` usage -- this yields a full rectangular matrix.

This is particularly useful for asymmetric similarities, such as containment.


**Use 2 -- Sketch \+ Top-k NN graphs**

For many applications, the quadratic space and time complexity of pairwise comparisons is prohibitive;
we use locality-sensitive hashing (LSH) table to only perform comparisons between candidates highly probable to be nearest neighbors.

To generate a K-NN graph between a collection of sequences, with k = 250:

```
dashing2 sketch <comparison options...> --cmpout <outfile> --topk 250 -F F.txt
```

The output file will be a table with the names and similarities for the nearest neighbors so generated; the corresponding binary output is Compressed Sparse Row-notation results.

Clustering: This can be used for spectral clustering community detection algorithms such as Louvain and Leiden.

**Use 3 -- Sketch \+ Jaccard-thresholded similarity graphs**

Alternative to `--topk [k]`, once can select a Jaccard similarity threshold below which the algorithm can ignore.

```
dashing2 sketch <comparison options...> --cmpout <outfile> --similarity-threshold 0.7 -F F.txt
```

The output file will be a table with the names and similarities for the nearest neighbors so generated; the corresponding binary output is Compressed Sparse Row-notation results.

Clustering: This can be used for spectral clustering community detection algorithms such as Louvain and Leiden.



**Use 4 -- Protein sequence similarity similarity**

The feature `--parse-by-seq` allows us to sketch and compute similarities between collections of sequences in a single file;
in particular, this is useful for sequence files of protein sequences.

The default similarity/distance is k-mer set comparisons in sketch space; however, edit distance is perhaps more useful for protein sequences.
Also, sketch sizes for proteins should be substantially smaller than for a full genome sequence.

```
dashing2 sketch -S256 --cmpout prot.k5.table --parse-by-seq -k5 {--protein,--protein14,--protein8,--protein6} uniref50.fa
```

If `--edit-distance` is enabled, OrderMinHash (Guillaume Marcais, et al.) is used to build an LSH table, and comparisons are generated between OrderMinHash sketches;
this is only allowed in parse-by-seq mode, as it is undefined on a collection of sequences as opposed to a single sequence.

```
dashing2 sketch --topk 25 --edit-distance --compute-edit-distance -S256 --cmpout prot.top25.k5.table --parse-by-seq -k5 {--protein,--protein14,--protein8,--protein6} uniref50.fa

```

If either `--refine-exact` and `--compute-edit-distance` is enabled, then the final distances will be computed via edit distance.

This is particularly useful if the result is Jaccard or top-k thresholded, allowing linear-time top-k nearest neighbor lists.


**Use 5: weighted sketching of feature sets**

There's also `dashing2 wsketch`, which can be used for hashing weighted sets for comparisons.
This is for the case where there are a set of integral identifiers for sketching.
See `dashing2 wsketch --help` for usage and examples.


### Installation

Easiest installation is `git clone --recursive https://github.com/dnbaker/dashing2 && make -j4`.

Dashing2 is written in C++17, and therefore needs a relatively recent compiler.


## Versions + Configuration

1. More than 2^32 items -
`dashing2` uses 32-bit hash identifiers in LSH tables for speed and memory efficiency.
To use more than 4.3 billion, use `dashing2-64`, which switched to 64-bit identifiers and hashes.

The default version of Dashing2 is dashing2, which uses 32-bit LSH keys and ID types in its NN tables;
this is faster and more memory-efficient, but less specific;

2. Hardware cache size
When comparing sketches, computations are grouped for better cache efficiency.
A group size is selected to fit as many sketches as possible in cache;
The default cache size estimate is 4MB. To change this, set the `D2_CACHE_SIZE` environment variable.

```sh
# set cache size to 64 MB
export D2_CACHE_SIZE=67108864
```

As an aside, these sketches are stored contiguously to reduce fragmentation compared to Dashing1.
