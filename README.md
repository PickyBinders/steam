# STEAM — Search with TEA against Many

STEAM is heavily adapted from [Foldseek](https://github.com/steineggerlab/foldseek) (Van Kempen et al., Nature Biotechnology 2024), replacing Foldseek's 3Di structural alphabet with [TEA](https://github.com/PickyBinders/tea). This means STEAM can be applied to any protein sequence, no 3D structure required. Like Foldseek, STEAM is built on the [MMseqs2](https://github.com/soedinglab/MMseqs2) framework.

If you have just a few proteins, consider using the [STEAM web-server](https://pickybinders.org/tea/steam) instead.

See [the preprint](https://doi.org/10.1101/2025.11.27.690975) to learn more.

## Requirements

- CMake >= 3.15
- GCC >= 7 or Clang
- For TEA sequence generation: [TEA](https://github.com/PickyBinders/tea) (`pip install git+https://github.com/PickyBinders/tea.git`)

## Installation

```bash
# Install build dependencies (if needed)
mamba install -c conda-forge cmake gxx_linux-64

# Build
git clone --branch dev --recursive https://github.com/PickyBinders/steam.git
cd steam
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
```

The binary will be at `build/src/steam`.

## Quick start

### 1. Generate TEA sequences

Convert amino acid sequences to the TEA structural alphabet using the `tea_convert` tool:

```bash
tea_convert -f proteins.fasta -o proteins_tea.fasta
```

This requires a GPU and the [TEA](https://github.com/PickyBinders/tea) package. The output is a FASTA file with TEA sequences in the same order as the input.

### 2. Search

```bash
steam easy-search query_tea.fasta query_aa.fasta \
                   target_tea.fasta target_aa.fasta \
                   result.m8 tmp
```

Useful flags:

| Flag | Default | Notes |
|-----------|---------|-------------|
| `-e` | 100 | E-value threshold |
| `--max-seqs` | 2000 | Maximum results per query from prefiltering |
| `--min-seq-id` | 0 | Minimum **amino acid** sequence identity |

### 3. Cluster

`easy-cluster` runs cascaded clustering (sensitive) and `easy-linclust` runs linear-time clustering (faster, less sensitive). Both take paired TEA/AA FASTAs:

```bash
# Cascaded clustering
steam easy-cluster proteins_tea.fasta proteins_aa.fasta clusterResult tmp

# Linear-time clustering (large datasets)
steam easy-linclust proteins_tea.fasta proteins_aa.fasta clusterResult tmp
```

Outputs three files alongside `clusterResult`:

- `clusterResult_cluster.tsv` — `<representative>	<member>` adjacency list
- `clusterResult_rep_seq.fasta` — one AA sequence per cluster representative
- `clusterResult_all_seqs.fasta` — FASTA grouped by cluster

Useful flags:

| Flag | Default | Notes |
|---|---|---|
| `--min-seq-id` | 0 | Minimum **amino acid** sequence identity for cluster members |
| `-c` | 0.8 | Minimum coverage |
| `-e` | 0.01 | E-value threshold |
| `--cov-mode` | 0 | 0=bidirectional, 1=target, 2=query |
| `--cluster-reassign` | off | Cascaded only: corrects criteria-violations from cascaded merging |
| `--single-step-cluster` | off | Cascaded only: skip cascading, single pass |

## Commands

| Command | Description |
|---------|-------------|
| `easy-search` | Search FASTA pairs against FASTA pairs or a pre-built database |
| `easy-cluster` | Cluster paired TEA/AA FASTAs (cascaded, sensitive) |
| `easy-linclust` | Cluster paired TEA/AA FASTAs (linear-time, faster) |
| `createdb` | Create a STEAM database from paired TEA/AA FASTA files |
| `computediversity` | Compute or refresh fixed MinHash metadata for an existing paired database |
| `search` | Search pre-built databases (faster for repeated searches) |
| `prefilter` | Generate native candidate rows directly |
| `cluster` | Cluster a pre-built database (cascaded) |
| `linclust` | Cluster a pre-built database (linear-time) |
| `convertalis` | Convert alignment results to various output formats |
| `createsubdb` | Subset a STEAM database (keeps `_aa` companion in sync) |

## Database workflow

For searching or clustering the same database multiple times, pre-build it:

```bash
# Create database (one time)
steam createdb target_tea.fasta target_aa.fasta targetDB

# Search against pre-built database (fast, repeatable)
steam easy-search query_tea.fasta query_aa.fasta targetDB result.m8 tmp

# Cluster the pre-built database
steam cluster targetDB clusterDB tmp
```

## Output format

Default BLAST-tab format (same as MMseqs2/BLAST -outfmt 6):

```text
query  target  fident  alnlen  mismatch  gapopen  qstart  qend  tstart  tend  evalue  bits
```

Custom output with `--format-output` adds TEA-specific output columns:

| Column | Description |
|--------|-------------|
| `tfident` | TEA fractional identity |
| `tpident` | TEA percent identity |
| `qteaseq` | Query TEA full sequence |
| `tteaseq` | Target TEA full sequence |
| `qteaaln` | Query TEA aligned sequence |
| `tteaaln` | Target TEA aligned sequence |

Standard MMseqs2 output columns (`fident`, `alnlen`, `qcov`, `tcov`, `evalue`, `raw`, `bits`, etc.) are also available.

## Scoring

The alignment score at each position is the sum of:

- **MATCHA score**: substitution score from the TEA alphabet matrix, scaled by `--tea-scale` (default 1; use 2 for half-bit matrices)
- **AA score**: BLOSUM62 substitution score, weighted by `--aa-weight` (default 3)

By default, prefilter candidates are ranked using the combined MATCHA +
weighted AA ungapped diagonal score before the `--max-seqs` cutoff. Pass
`--ungapped-tea-aa 0` for TEA-only ungapped ranking. K-mer matching still uses
TEA sequences.

For an exported TEA matrix, verify the native loaded integers before searching:

```bash
cmake --build build --target steam_dump_matrix
build/src/steam_dump_matrix tea.out 2
```
After alignment traceback, the score becomes $floor(S*(1+e)+0.5)$, where `e` is the mean Markov rarity weight of distinct query TEA words supported by gap-free seven-column seed spans.

## E-value computation

STEAM uses a continuous piecewise-loglinear E-value model fitted to the
empirical false-positive score tail. E-values are computed as:

```
E(s) = D_target * 10^(c + m_low*s + (m_high-m_low)*max(0, s-b))
```

where `s` is the final corrected score, `D_target` is the complete target
database's AA3+TEA5 MinHash effective target count, and `b` is the fitted
score breakpoint.

The selected coefficients are:

```text
m_low = -0.005628286904365784
m_high = -0.0013916072449292318
c = 0.24858538629087792
b = 578
```

### Legacy TEA

Pass the legacy matrix and gap, and disable the seed correction. For example:

```bash
steam easy-search query_tea.fasta query_aa.fasta targetDB result.m8 tmp \\
  --matcha legacy_matrix.out --seed-correction 0 \\
  --tea-scale 1 --aa-weight 1.4 --gap-open 14 --gap-extend 2 \\
  --comp-bias-corr 0 --ungapped-tea-aa 0 -e inf
```
