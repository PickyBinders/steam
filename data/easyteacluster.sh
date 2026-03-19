#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# Create TEA DB if input is FASTA pair
if notExists "${TMP_PATH}/input.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${INPUT_TEA}" "${INPUT_AA}" "${TMP_PATH}/input" ${CREATEDB_PAR} \
        || fail "createteadb died"
fi
INPUT="${TMP_PATH}/input"

INTERMEDIATE="${TMP_PATH}/clu"
if notExists "${INTERMEDIATE}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" cluster "${INPUT}" "${INTERMEDIATE}" "${TMP_PATH}/cluster_tmp" ${CLUSTER_PAR} \
        || fail "Cluster died"
fi

# shellcheck disable=SC2086
"$MMSEQS" createtsv "${INPUT}" "${INPUT}_h" "${INTERMEDIATE}" "${RESULTS}_cluster.tsv" ${THREADS_PAR} \
    || fail "createtsv died"

# shellcheck disable=SC2086
"$MMSEQS" result2repseq "${INPUT}" "${INTERMEDIATE}" "${TMP_PATH}/rep" ${RESULT2REPSEQ_PAR} \
    || fail "result2repseq died"

# shellcheck disable=SC2086
"$MMSEQS" result2flat "${INPUT}" "${INPUT}_h" "${TMP_PATH}/rep" "${RESULTS}_rep_seq.fasta" --use-fasta-header ${VERBOSITY} \
    || fail "result2flat died"

# shellcheck disable=SC2086
"$MMSEQS" createseqfiledb "${INPUT}" "${INTERMEDIATE}" "${TMP_PATH}/clu_seqs" ${THREADS_PAR} \
    || fail "createseqfiledb died"

# shellcheck disable=SC2086
"$MMSEQS" result2flat "${INPUT}" "${INPUT}_h" "${TMP_PATH}/clu_seqs" "${RESULTS}_all_seqs.fasta" ${VERBOSITY} \
    || fail "result2flat died"

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/rep" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/clu_seqs" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${INTERMEDIATE}" ${VERBOSITY}
    if [ -f "${TMP_PATH}/input.dbtype" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/input" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/input_h" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/input_aa" ${VERBOSITY}
    fi
    rm -rf "${TMP_PATH}/cluster_tmp"
    rm -f "${TMP_PATH}/easyteacluster.sh"
fi
