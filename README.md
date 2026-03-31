# STEAM — Search with TEA against Many

STEAM is heavily adapted from [Foldseek](https://github.com/steineggerlab/foldseek) (Van Kempen et al., Nature Biotechnology 2024), replacing Foldseek's 3Di structural alphabet (derived from backbone coordinates) with [TEA](https://github.com/PickyBinders/tea). This means STEAM can be applied to any protein sequence, no 3D structure required. Like Foldseek, STEAM is built on the [MMseqs2](https://github.com/soedinglab/MMseqs2) framework.

## How it works

STEAM uses paired representations for each protein: a TEA sequence encoding structural features learned from ESM2, and the standard amino acid sequence. Both are scored simultaneously during alignment using the MATCHA substitution matrix (for TEA) and BLOSUM62 (for amino acids), providing better sensitivity than either signal alone.

### Search pipeline

1. **Prefiltering**: Fast k-mer matching on TEA sequences (k=6, spaced pattern) to find candidate pairs
2. **Ungapped rescoring**: TEA+AA ungapped diagonal alignment filters low-scoring hits (coverage-weighted score threshold, default 10). Disable with `--min-ungapped-score 0`.
3. **Gapped alignment**: Smith-Waterman alignment with combined TEA+AA scoring
4. **E-value estimation**: Log-linear statistical model calibrated on SCOP structural classification

### Clustering pipeline

Clustering follows the MMseqs2/Foldseek cascaded approach:

1. **Linclust** (linear-time initial clustering):
   - k-mer matching with auto-selected k and reduced alphabet
   - TEA+AA ungapped rescoring
   - Pre-clustering from ungapped scores
   - Gapped SW alignment on cluster representative-member pairs only
   - Merge pre-clusters with gapped clusters
2. **Cascaded refinement** (optional, iterative):
   - Prefilter remaining sequences with increasing sensitivity
   - Gapped alignment and clustering at each step
   - Default 3 cascade steps

## Requirements

- CMake >= 3.15
- GCC >= 7 or Clang
- For TEA sequence generation: [TEA](https://github.com/PickyBinders/tea) (`pip install git+https://github.com/PickyBinders/tea.git`)

## Installation

```bash
git clone --recursive https://github.com/PickyBinders/steam.git
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

### 3. Cluster

```bash
steam easy-cluster sequences_tea.fasta sequences_aa.fasta \
                    result tmp \
                    -c 0.9 -e 0.01
```

### 4. Interpret results

Results are sorted by E-value (best first). An E-value of 10 means approximately 10 false positive hits are expected per query at that score threshold. Hits with E < 1 are generally significant.

## Commands

| Command | Description |
|---------|-------------|
| `easy-search` | Search FASTA pairs against FASTA pairs or a pre-built database |
| `easy-cluster` | Cluster paired TEA+AA FASTA files |
| `createdb` | Create a STEAM database from paired TEA/AA FASTA files |
| `search` | Search pre-built databases (faster for repeated searches) |
| `cluster` | Cluster a pre-built database |
| `convertalis` | Convert alignment results to various output formats |
| `tearescorediagonal` | Ungapped TEA+AA diagonal rescoring |

## Database workflow

For searching the same target database multiple times, pre-build it:

```bash
# Create database (one time)
steam createdb target_tea.fasta target_aa.fasta targetDB

# Search against pre-built database (fast, repeatable)
steam easy-search query_tea.fasta query_aa.fasta targetDB result.m8 tmp
```

The `createdb` command reads TEA and AA FASTA files in lockstep -- headers must match in order. TEA headers with entropy suffixes from `tea_convert` (e.g., `>seq1|H=0.123`) are stripped automatically for matching.

## Key parameters

### Search

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--aa-weight` | 1.4 | Weight for amino acid scores in combined scoring (0 = TEA only) |
| `--matcha` | matcha.out (bundled) | MATCHA substitution matrix for TEA scoring |
| `-e` | 100 | E-value threshold |
| `--max-seqs` | 2000 | Maximum results per query from prefiltering |
| `--gap-open` | 14 | Gap open penalty |
| `--gap-extend` | 2 | Gap extension penalty |
| `-k` | 6 | k-mer length for prefiltering |
| `--exact-kmer-matching` | 1 | Use exact k-mer matching |
| `--min-ungapped-score` | 10 | Minimum ungapped covScore to pass to gapped alignment (0 = disable) |

### Clustering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-c` | 0.8 | Coverage threshold |
| `-e` | 0.01 | E-value threshold |
| `--min-seq-id` | 0 | Minimum sequence identity (also controls auto-sensitivity) |
| `--cluster-steps` | 3 | Number of cascaded clustering iterations |
| `--cluster-mode` | auto | Clustering algorithm (auto-selected: set-cover or greedy) |
| `--cov-mode` | 0 | Coverage mode (0=bidirectional, 1=target, 2=query) |
| `--single-step-clustering` | false | Skip cascaded steps, use linclust only |

## Output format

Default BLAST-tab format (same as MMseqs2/BLAST -outfmt 6):

```text
query  target  fident  alnlen  mismatch  gapopen  qstart  qend  tstart  tend  evalue  bits
```

Custom output with `--format-output`:

```bash
steam easy-search ... --format-output query,target,evalue,fident,alnlen,qcov,tcov
```

### TEA-specific output columns

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

- **MATCHA score**: substitution score from the TEA alphabet matrix
- **AA score**: BLOSUM62 substitution score, weighted by `--aa-weight` (default 1.4)

Results are ranked by a coverage-weighted score: `raw_score * sqrt(min(query_coverage, target_coverage))`. This penalizes partial alignments that only match a small fragment while preserving sensitivity for full-domain matches.

## E-value computation

STEAM uses a log-linear E-value model following [Edgar & Sahakyan (2025)](https://doi.org/10.1101/2025.07.17.665375), which accounts for the heavy-tailed false positive score distribution observed in structural alphabet searches. E-values are computed as:

```
E(s) = (H/Q) * 10^(m*s + c)
```

where `s` is the coverage-weighted alignment score, `H/Q` is the average number of reported hits per query (computed at runtime from prefilter results), and `m` and `c` are parameters fitted on SCOP40c. Unlike traditional Karlin-Altschul E-values, this model does not assume a Gumbel distribution.
