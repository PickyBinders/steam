#!/bin/bash

set -euo pipefail

steam_bin="${1:-build/src/steam}"
test_dir="${2:-$(mktemp -d /tmp/steam_evalue_minhash.XXXXXX)}"
fixture_dir="tests/data/ungapped_tea_aa"
mkdir -p "${test_dir}"

"${steam_bin}" createdb \
    "${fixture_dir}/target_tea.fasta" "${fixture_dir}/target_aa.fasta" \
    "${test_dir}/target" --threads 2 -v 1
test -s "${test_dir}/target.steam-diversity"
grep -q $'^method\taa3-tea5-doph64-hll14-v1$' \
    "${test_dir}/target.steam-diversity"
grep -q $'^sequences\t2$' "${test_dir}/target.steam-diversity"

# Exact duplication changes the number of target rows but not database
# diversity.  Rename headers so the native database still has unique IDs.
for copy in $(seq 1 10); do
    awk -v copy="${copy}" '/^>/{print $0 "_copy" copy; next} {print}' \
        "${fixture_dir}/target_tea.fasta" >> "${test_dir}/duplicate_tea.fasta"
    awk -v copy="${copy}" '/^>/{print $0 "_copy" copy; next} {print}' \
        "${fixture_dir}/target_aa.fasta" >> "${test_dir}/duplicate_aa.fasta"
done
"${steam_bin}" createdb \
    "${test_dir}/duplicate_tea.fasta" "${test_dir}/duplicate_aa.fasta" \
    "${test_dir}/duplicate" --threads 2 -v 1
grep -q $'^sequences\t20$' "${test_dir}/duplicate.steam-diversity"

base_diversity="$(awk -F '\t' '$1 == "effective_targets" {print $2}' \
    "${test_dir}/target.steam-diversity")"
duplicate_diversity="$(awk -F '\t' '$1 == "effective_targets" {print $2}' \
    "${test_dir}/duplicate.steam-diversity")"
test "${base_diversity}" = "${duplicate_diversity}"

# The explicit command must reproduce createdb's automatic metadata exactly.
cp "${test_dir}/target.steam-diversity" "${test_dir}/target.expected"
rm "${test_dir}/target.steam-diversity"
"${steam_bin}" computediversity "${test_dir}/target" --threads 2 -v 1
cmp "${test_dir}/target.expected" "${test_dir}/target.steam-diversity"

# Subsetting changes database diversity, so the paired createsubdb wrapper must
# recompute metadata rather than copy the parent's value.
awk -F '\t' 'NR == 1 {print $1}' "${test_dir}/target.lookup" \
    > "${test_dir}/subset.list"
"${steam_bin}" createsubdb "${test_dir}/subset.list" \
    "${test_dir}/target" "${test_dir}/subset" --threads 2 -v 1
test -s "${test_dir}/subset.steam-diversity"
grep -q $'^sequences\t1$' "${test_dir}/subset.steam-diversity"

run_case() {
    name="$1"
    target="$2"
    mkdir -p "${test_dir}/${name}"
    STEAM_BUCKET_COVERAGE_PCT=0 "${steam_bin}" easy-search \
        "${fixture_dir}/query_tea.fasta" "${fixture_dir}/query_aa.fasta" \
        "${target}" "${test_dir}/${name}.m8" "${test_dir}/${name}" \
        --threads 2 --max-seqs 20 --ungapped-tea-aa 0 -e 1e30 \
        --loglinear-m -0.008074326697592253 \
        --loglinear-m-high -0.0028050667758289924 \
        --loglinear-breakpoint 431 --loglinear-c 0.1632259424343494 \
        --format-output query,target,evalue,raw -v 1
}

run_case minhash_base "${test_dir}/target"
run_case minhash_duplicate "${test_dir}/duplicate"
awk -F '\t' '$1 == "query" && $2 == "tea_best" {print $3, $4}' \
    "${test_dir}/minhash_base.m8" > "${test_dir}/minhash_base.common"
awk -F '\t' '$1 == "query" && $2 == "tea_best_copy1" {print $3, $4}' \
    "${test_dir}/minhash_duplicate.m8" > "${test_dir}/minhash_duplicate.common"
test -s "${test_dir}/minhash_base.common"
cmp "${test_dir}/minhash_base.common" "${test_dir}/minhash_duplicate.common"

# Calibrated E-values must fail closed rather than silently reverting to target count.
rm "${test_dir}/duplicate.steam-diversity"
if run_case missing_metadata "${test_dir}/duplicate" \
        >"${test_dir}/missing_metadata.log" 2>&1; then
    echo "MinHash search unexpectedly accepted missing metadata" >&2
    exit 1
fi
grep -q 'Cannot compute diversity-adjusted E-values' \
    "${test_dir}/missing_metadata.log"

printf 'MinHash database-diversity E-value tests passed\n'
