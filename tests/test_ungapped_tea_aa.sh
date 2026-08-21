#!/bin/bash

set -euo pipefail

steam_bin="${1:-build/src/steam}"
test_dir="${2:-$(mktemp -d /tmp/steam_ungapped_tea_aa.XXXXXX)}"
fixture_dir="tests/data/ungapped_tea_aa"

mkdir -p "${test_dir}/flag_off" "${test_dir}/flag_on"

"${steam_bin}" easy-search \
    "${fixture_dir}/query_tea.fasta" "${fixture_dir}/query_aa.fasta" \
    "${fixture_dir}/target_tea.fasta" "${fixture_dir}/target_aa.fasta" \
    "${test_dir}/flag_off.m8" "${test_dir}/flag_off" \
    --max-seqs 1 --threads 1 -v 1

"${steam_bin}" easy-search \
    "${fixture_dir}/query_tea.fasta" "${fixture_dir}/query_aa.fasta" \
    "${fixture_dir}/target_tea.fasta" "${fixture_dir}/target_aa.fasta" \
    "${test_dir}/flag_on.m8" "${test_dir}/flag_on" \
    --ungapped-tea-aa 1 --max-seqs 1 --threads 1 -v 3

test -s "${test_dir}/flag_off.m8"
test -s "${test_dir}/flag_on.m8"
grep -q $'^query\ttea_best\t' "${test_dir}/flag_off.m8"
grep -q $'^query\tcombined_best\t' "${test_dir}/flag_on.m8"

printf 'combined ungapped smoke test passed\n'
printf 'without flag: '
cut -f1-2 "${test_dir}/flag_off.m8"
printf 'with flag:    '
cut -f1-2 "${test_dir}/flag_on.m8"
