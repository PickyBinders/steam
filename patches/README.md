# STEAM patches for the bundled mmseqs2 submodule

Each `*.patch` in this directory is applied to `lib/mmseqs/` at CMake
configure time (see the patch-apply block in the top-level `CMakeLists.txt`).

## Existing patches

### `0001-steam-bucket-filter.patch`

Adds a hot-k-mer skip filter in `QueryMatcher::matchQuery`. On the TEA
alphabet a small number of structural-motif k-mers (e.g. `HWHACW`) accumulate
huge index buckets that dominate prefilter compute without contributing
discriminative signal. The filter skips query lookups for k-mers whose target
bucket falls within the top *X*% of the index by cumulative entry coverage.

Configurable via one environment variable (no new CLI flag so the patch stays
self-contained): `STEAM_BUCKET_COVERAGE_PCT` (default `20`). The threshold is
derived once per index from the sorted bucket-size distribution.
