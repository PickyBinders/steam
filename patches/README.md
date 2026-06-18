# STEAM patches for the bundled mmseqs2 submodule

Each `*.patch` in this directory is applied to `lib/mmseqs/` at CMake
configure time (see the patch-apply block in the top-level `CMakeLists.txt`).

## Existing patches

### `0001-steam-bucket-filter.patch`

**Target:** `src/prefiltering/QueryMatcher.cpp::matchQuery`.

Adds a hot-k-mer skip filter to the prefilter. On the TEA alphabet a small
number of structural-motif k-mers (e.g. `HWHACW`) accumulate huge index
buckets that dominate prefilter compute without contributing discriminative
signal. The filter skips query lookups for k-mers whose target bucket falls
within the top *X*% of the index by cumulative entry coverage.

Configurable via one environment variable (no new CLI flag so the patch stays
self-contained): `STEAM_BUCKET_COVERAGE_PCT` (default `20`; set `0` to disable).
The threshold is derived once per index from the sorted bucket-size
distribution.

### `0002-steam-cluster-prefilter-submat.patch`

**Target:** `src/workflow/Cluster.cpp` (cascaded cluster) and
`src/workflow/Linclust.cpp` (linclust).

Lets the cluster/linclust workflows use different substitution matrices for
prefilter and alignment. STEAM's prefilter runs on TEA (matrix `matcha.out`)
while gapped alignment runs on the AA alphabet (BLOSUM62).
The patch wraps the prefilter and kmermatcher parameter-string generation in
`STEAM_SWAP_PREFILTER_SUBMAT` / `STEAM_RESTORE_PREFILTER_SUBMAT` macros that
swap `par.scoringMatrixFile` to the value of `MMSEQS_PREFILTER_SUBMAT` and
back. When the env var is unset the macros are no-ops.

Swap points: 4 in `Cluster.cpp` (`PREFILTER0_PAR`, the cascaded `PREFILTERi_PAR`
loop, `PREFILTER_REASSIGN_PAR`, and single-step `PREFILTER_PAR`) and 1 in
`Linclust.cpp` (`KMERMATCHER_PAR`).

Driven from STEAM by the `teacluster` / `teaLinclust` wrappers in
`src/workflow/TeaClusterWrappers.cpp`, which `setenv("MMSEQS_PREFILTER_SUBMAT",
"matcha.out", 1)` before delegating to upstream `clusteringworkflow()` /
`linclust()`.

If `STEAM_CLUSTER_DEFAULTS=1` (default), `setWorkflowDefaults` in `Cluster.cpp` switches to
foldseek-style defaults (`--max-seqs 1000` and `-e 0.01`).

### `0003-steam-linclust-bucket-filter.patch`

**Target:** `src/linclust/kmermatcher.cpp::assignGroup`.

Linclust analogue of `0001`. Linclust does its own k-mer matching outside of
`QueryMatcher`, so the prefilter patch doesn't apply. This patch pre-scans the
sorted `hashSeqPair` array to build a bucket-size distribution, picks the
threshold that yields the top *X*% cumulative coverage of k-mer entries in the
split, then refuses to emit pairs from any bucket larger than the threshold
(entries marked `SIZE_T_MAX`).

Same environment variable as `0001`: `STEAM_BUCKET_COVERAGE_PCT` (default
`20`; set `0` to disable). A floor of 10 entries per bucket prevents the
filter from triggering when the bucket-size distribution is too uniform for
the top-X% concept to be meaningful.

If `STEAM_CLUSTER_DEFAULTS=1` (default), `setLinclustWorkflowDefaults` in `Linclust.cpp` switches to
foldseek-style defaults (`-e 0.01` and `kmersPerSequence 300`).
