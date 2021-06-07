Draft outline for Dashing2

# Encoding kinds, data formats, and options
By default, we do set sketching, but we can do multiset (BagMinHash) or discrete prob set (ProbMinHash) sketching.

We support fastx files (fastq and fasta), several k-mer/m-mer encoding schemes, and minimizers.
For those, we do set-based sketching (default) unless counting is enabled.

For counting, we either use a (single-row) CountSketch to approximate weighted set comparisons or we compute exact counts using a hash table.

### Completed
#### BigWig sketching
BigWigs done.

If opts.sspace_ set to SPACE\_SET yields SetSketches (one-permutation or full, depending on opts.one_perm_).
This flattens all genomic locations to single count (presence/absence).
Otherwise, it yields BagMinHash for SPACE\_MULTISET and ProbMinHash for SPACE\_PSET.

#### BED sketching
For BED files, we only support sets currently.
By default, BED files are sketched as sets. This can be changed to MULTISET or PSET, but this is more expensive, as it requires coverage counting first.
##### Normalization
Default normalization treats each base pair in an interval as equally weighted. For SPACE\_SET, this doesn't matter, but for weighted sketching, it does.

That is, an interval from 0:100-200 will have twice the weight as 0:100-150. If, instead, you want to treat each interval as weighing the same as another interval, enable ParseOptions.bed_parse_normalize_intervals_.


#### LeafCutter sketching
Similar to BED, except that the normalize flag causes the splicing event weight to be the fraction of reads supporting the junction rather than the absolute number  of reads.

IE, 3/5 would have weight 0.6 if bed_parse_normalize_intervals_ is enabled and weight 3 otherwise.

#### Sequence sketching

This can be done per file or per sequence. Per sequence is not yet supported, but we plan to.

Sketches can be:
    SetSketch (OPC or C)
    BagMinHash
    ProbMinHash
    OrderMinHash

However, OrderMinHash is only available for sequences

When handling multiplicities, you can use exact counting, which may be slower, or you can approximate the count vectors with a single-row count-sketch
by setting opts.cssize_ > 0.

ProbMinHash is usually significantly (20x++) faster than BagMinHash, although multiset jaccard may be more appropriate for some problems.
