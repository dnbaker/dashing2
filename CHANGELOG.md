# Changelog

## v2.1.21 (2026-10-20)

### Fixed
- `sketch --set --parse-by-seq` segfaulted: parse-by-seq built only SetSketch signatures while
  `--set` sends comparison to the exact on-disk path (whose file list this mode never fills).
  Each sequence's sorted, deduplicated k-mer set is now kept in memory and compared exactly;
  values match file-mode `--set`.
- setsketch trim loops (erase while iterating) skipped elements and, with a vector-backed map,
  could read past `end()`.
- `flat_hash_set` alias dropped its `Hash` parameter.
- bwsketch: the parallel contig loop only looks up, never inserts.

### Changed
- `flat_hash_map`/`flat_hash_set` use `ankerl::unordered_dense` (robin_hood removed).
- bonsai v0.6.0, sketch v0.20.0 (densification, sparse Jaccard and counting b-bit MinHash
  fixes; explicit `add_threadsafe`). Sketch outputs on test.fa (one-permutation and full
  SetSketch, BagMinHash, ProbMinHash, both with count sketch) are byte-identical to v2.1.20.
