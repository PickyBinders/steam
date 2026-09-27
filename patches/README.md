# STEAM patches for bundled MMseqs2

CMake applies these patches in filename order to the pinned MMseqs2 submodule.
Already-applied patches are detected by a reverse check; incompatible source
versions fail configuration. Apply the complete set before compiling.

| Patch | Purpose |
| --- | --- |
| `0001-steam-bucket-filter.patch` | Per-database frequent-seed thresholds and exclusion in QueryMatcher; fixed 10% default. |
| `0002-steam-cluster-prefilter-submat.patch` | Keep TEA prefilter matrices separate from AA alignment matrices in existing clustering workflows. |
| `0003-steam-linclust-bucket-filter.patch` | Frequent-seed exclusion in linclust, which has its own matcher; 10% default. |
| `0004-combined-tea-aa-ungapped.patch` | Companion-AA database and matrix loading, plus combined TEA+AA diagonal scoring before candidate truncation. |
| `0005-preserve-spaced-pattern-index.patch` | Write nonempty spaced patterns into cached-index metadata, preserving k5/k6 seed selection. |
| `0006-split-invariant-prefilter.patch` | Complete-target seed counts across splits, global candidate caps after split merging, and deterministic tie ordering. |

Patch `0006` depends on the bucket-filter interfaces in `0001` and the
companion-AA interfaces in `0004`. Apply the complete series in filename order.

`STEAM_BUCKET_COVERAGE_PCT` sets frequent-seed exclusion by cumulative
seed-occurrence coverage (default 10%). STEAM searches supply
`STEAM_SPLIT_INVARIANT=1` and the companion-AA scoring inputs internally.
The alignment-supported score correction and MinHash E-value model live in
STEAM itself rather than these dependency patches.
