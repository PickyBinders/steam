#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
}

# Steam's createdb takes <tea_fasta> <aa_fasta> <out_db>, not a variadic
# fasta list, so we cannot reuse mmseqs's easycluster.sh as-is.
if notExists "${TMP_PATH}/input.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${INPUT_TEA}" "${INPUT_AA}" "${TMP_PATH}/input" ${CREATEDB_PAR} \
        || fail "createdb died"
fi

if notExists "${TMP_PATH}/clu.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" "${CLUSTER_MODULE}" "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/clu_tmp" ${CLUSTER_PAR} \
        || fail "Clustering died"
fi

if notExists "${TMP_PATH}/cluster.tsv"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtsv "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/cluster.tsv" ${THREADS_PAR} \
        || fail "createtsv died"
fi

if notExists "${TMP_PATH}/rep_seq.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" result2repseq "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/clu_rep" ${RESULT2REPSEQ_PAR} \
        || fail "result2repseq died"
    # shellcheck disable=SC2086
    "$MMSEQS" result2flat "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/clu_rep" "${TMP_PATH}/rep_seq.fasta" --use-fasta-header ${VERBOSITY_PAR} \
        || fail "result2flat died"
fi

if notExists "${TMP_PATH}/all_seqs.fasta"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createseqfiledb "${TMP_PATH}/input" "${TMP_PATH}/clu" "${TMP_PATH}/clu_seqs" ${THREADS_PAR} \
        || fail "createseqfiledb died"
    # shellcheck disable=SC2086
    "$MMSEQS" result2flat "${TMP_PATH}/input" "${TMP_PATH}/input" "${TMP_PATH}/clu_seqs" "${TMP_PATH}/all_seqs.fasta" ${VERBOSITY_PAR} \
        || fail "result2flat died"
fi

mv "${TMP_PATH}/all_seqs.fasta" "${RESULTS}_all_seqs.fasta"
mv "${TMP_PATH}/rep_seq.fasta"  "${RESULTS}_rep_seq.fasta"
mv "${TMP_PATH}/cluster.tsv"    "${RESULTS}_cluster.tsv"

if [ -n "${REMOVE_TMP}" ]; then
    "$MMSEQS" rmdb "${TMP_PATH}/input"    ${VERBOSITY_PAR}
    "$MMSEQS" rmdb "${TMP_PATH}/input_h"  ${VERBOSITY_PAR}
    "$MMSEQS" rmdb "${TMP_PATH}/input_aa" ${VERBOSITY_PAR}
    "$MMSEQS" rmdb "${TMP_PATH}/clu_seqs" ${VERBOSITY_PAR}
    "$MMSEQS" rmdb "${TMP_PATH}/clu_rep"  ${VERBOSITY_PAR}
    "$MMSEQS" rmdb "${TMP_PATH}/clu"      ${VERBOSITY_PAR}
    rm -rf "${TMP_PATH}/clu_tmp"
    rm -f  "${TMP_PATH}/easyteacluster.sh"
fi
