Drafting for Dashing2

# Encoding kinds, data formats, and options
By default, we do set sketching, but we can do multiset (BagMinHash) or discrete prob set (ProbMinHash) sketching.

We support fastx files (fastq and fasta), several k-mer encoding schemes, and minimizers.
For those, we do set-based sketching (default) unless counting is enabled.

For counting, we either use a (single-row) CountSketch to approximate weighted set comparisons or we compute exact counts using a hash table.

For BigWigs, this is already done (as it is coverage summarization). This can be flattened to binary presence/absence if desired. (Select SET space.)

For BED files, we only support sets currently.

# After sketching, what can we do?

1. Build LSH index for set space, multiset space, or edit distance
2. In the process, build up top-k/nearest-neighbor lists.
3. Create candidate top-k/nn-lists
4. Perform sketch comparisons for final list of neighbors.
5. Emit (distance matrix or neighbor list or list of final genomes after de-duplication)
