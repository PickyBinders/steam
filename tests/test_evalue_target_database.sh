#!/bin/bash

set -euo pipefail

steam_bin="${1:-build/src/steam}"
test_dir="${2:-$(mktemp -d /tmp/steam_evalue_target_database.XXXXXX)}"
fixture_dir="tests/data/ungapped_tea_aa"

mkdir -p "${test_dir}"

cat "${fixture_dir}/query_tea.fasta" > "${test_dir}/query_batch_tea.fasta"
sed 's/^>query$/>query_copy/' "${fixture_dir}/query_tea.fasta" >> "${test_dir}/query_batch_tea.fasta"
cat "${fixture_dir}/query_aa.fasta" > "${test_dir}/query_batch_aa.fasta"
sed 's/^>query$/>query_copy/' "${fixture_dir}/query_aa.fasta" >> "${test_dir}/query_batch_aa.fasta"

run_case() {
    name="$1"
    query_tea="$2"
    query_aa="$3"
    max_seqs="$4"
    coverage="$5"
    pattern="$6"
    mkdir -p "${test_dir}/${name}"
    STEAM_BUCKET_COVERAGE_PCT="${coverage}" "${steam_bin}" easy-search \
        "${query_tea}" "${query_aa}" \
        "${fixture_dir}/target_tea.fasta" "${fixture_dir}/target_aa.fasta" \
        "${test_dir}/${name}.m8" "${test_dir}/${name}" \
        --ungapped-tea-aa 0 --max-seqs "${max_seqs}" --threads 1 -e 1e30 \
        --loglinear-m -0.01 --loglinear-m-high -0.002 \
        --loglinear-breakpoint 50 --loglinear-c 0 --p-fp 1 \
        --spaced-kmer-pattern "${pattern}" \
        --format-output query,target,evalue,raw -v 1
    awk -F '\t' '$1 == "query" && $2 == "tea_best" {print $3, $4}' \
        "${test_dir}/${name}.m8" > "${test_dir}/${name}.common"
    test -s "${test_dir}/${name}.common"
}

"${steam_bin}" createdb \
    "${fixture_dir}/target_tea.fasta" "${fixture_dir}/target_aa.fasta" \
    "${test_dir}/expected_target" --threads 1 -v 1
effective_targets="$(awk -F '\t' '$1 == "effective_targets" {print $2}' \
    "${test_dir}/expected_target.steam-diversity")"
test -n "${effective_targets}"

run_case baseline "${fixture_dir}/query_tea.fasta" "${fixture_dir}/query_aa.fasta" 2 0 1101101
run_case max_seqs_1 "${fixture_dir}/query_tea.fasta" "${fixture_dir}/query_aa.fasta" 1 0 1101101
run_case query_batch "${test_dir}/query_batch_tea.fasta" "${test_dir}/query_batch_aa.fasta" 2 0 1101101
run_case alternate_prefilter "${fixture_dir}/query_tea.fasta" "${fixture_dir}/query_aa.fasta" 2 20 1101101

for name in max_seqs_1 query_batch alternate_prefilter; do
    cmp "${test_dir}/baseline.common" "${test_dir}/${name}.common"
done

awk -v effective_targets="${effective_targets}" '
    {
        raw = $2
        exponent = -0.01 * raw
        if (raw > 50) exponent += (-0.002 - -0.01) * (raw - 50)
        expected = effective_targets * exp(log(10) * exponent)
        relative = ($1 > expected ? $1 - expected : expected - $1) / expected
        if (relative > 0.001) exit 1
    }
' "${test_dir}/baseline.common"
printf 'target-database diversity E-value invariance test passed\n'
