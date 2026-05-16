# STEAM — Search with TEA against Many

STEAM is heavily adapted from [Foldseek](https://github.com/steineggerlab/foldseek) (Van Kempen et al., Nature Biotechnology 2024), replacing Foldseek's 3Di structural alphabet (derived from backbone coordinates) with [TEA](https://github.com/PickyBinders/tea). This means STEAM can be applied to any protein sequence, no 3D structure required. Like Foldseek, STEAM is built on the [MMseqs2](https://github.com/soedinglab/MMseqs2) framework.

## Requirements

- CMake >= 3.15
- GCC >= 7 or Clang
- For TEA sequence generation: [TEA](https://github.com/PickyBinders/tea) (`pip install git+https://github.com/PickyBinders/tea.git`)

## Installation

```bash
# Install build dependencies (if needed)
mamba install -c conda-forge cmake gxx_linux-64

# Build
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

## Commands

| Command | Description |
|---------|-------------|
| `easy-search` | Search FASTA pairs against FASTA pairs or a pre-built database |
| `createdb` | Create a STEAM database from paired TEA/AA FASTA files |
| `search` | Search pre-built databases (faster for repeated searches) |
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

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-e` | 100 | E-value threshold |
| `--max-seqs` | 2000 | Maximum results per query from prefiltering |

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

- **MATCHA score**: substitution score from the TEA alphabet matrix
- **AA score**: BLOSUM62 substitution score, weighted by `--aa-weight` (default 1.4)

## E-value computation

STEAM uses a log-linear E-value model following [Edgar & Sahakyan (2025)](https://doi.org/10.1101/2025.07.17.665375). E-values are computed as:

```
E(s) = (H/Q) * 10^(m*s + c)
```

where `s` is the raw alignment score, `H/Q` is the average number of reported hits per query (computed at runtime from prefilter results), and `m` and `c` are parameters fitted on SCOP40c.
