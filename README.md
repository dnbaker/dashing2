## Input formats

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
BagMinhash is probably most appropriate for genome assemblies (--multiset), and ProbMinhash (--probs) is most appropriate for splicing and expression datasets, due to normalization.

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

However, OrderMinHash is only available for sequences

However, OrderMinHash is only available for sequences

When handling multiplicities, you can use exact counting, which may be slower, or you can approximate the count vectors with a single-row count-sketch
by setting opts.cssize\_ > 0.

ProbMinHash is usually significantly (20x++) faster than BagMinHash, although multiset jaccard may be more appropriate for some problems.

