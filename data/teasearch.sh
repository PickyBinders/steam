#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# 1. Prefilter (k-mer matching on TEA sequences)
if notExists "${TMP_PATH}/pref.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefilter "${QUERY}" "${TARGET}${INDEXEXT}" "${TMP_PATH}/pref" ${PREFILTER_PAR} \
        || fail "Prefilter step died"
fi

# 2. TEA+AA combined alignment
if notExists "${TMP_PATH}/aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" align "${QUERY}" "${TARGET}${INDEXEXT}" "${TMP_PATH}/pref" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
        || fail "TEA alignment step died"
fi

# shellcheck disable=SC2086
"$MMSEQS" mvdb "${TMP_PATH}/aln" "${RESULTS}" ${VERBOSITY}

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
fi
