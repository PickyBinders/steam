#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# 1. Prefilter (k-mer matching on TEA sequences) or exhaustive all-pairs
if notExists "${TMP_PATH}/pref.dbtype"; then
    if [ -n "${EXHAUSTIVE}" ]; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" prefilter "${QUERY}" "${TARGET}${INDEXEXT}" "${TMP_PATH}/pref" ${EXHAUSTIVE_PREFILTER_PAR} \
            || fail "Exhaustive prefilter step died"
    else
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" prefilter "${QUERY}" "${TARGET}${INDEXEXT}" "${TMP_PATH}/pref" ${PREFILTER_PAR} \
            || fail "Prefilter step died"
    fi
fi

ALIGN_INPUT="${TMP_PATH}/pref"

# 2. Ungapped TEA+AA rescoring (optional speed filter, disabled by default)
if [ -n "${RESCORE_PAR}" ]; then
    if notExists "${TMP_PATH}/pref_rescore.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$MMSEQS" tearescorediagonal "${QUERY}" "${TARGET}${INDEXEXT}" "${TMP_PATH}/pref" "${TMP_PATH}/pref_rescore" ${RESCORE_PAR} \
            || fail "TEA ungapped rescoring step died"
    fi
    ALIGN_INPUT="${TMP_PATH}/pref_rescore"
fi

# 3. TEA+AA gapped alignment
if notExists "${TMP_PATH}/aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" align "${QUERY}" "${TARGET}${INDEXEXT}" "${ALIGN_INPUT}" "${TMP_PATH}/aln" ${ALIGNMENT_PAR} \
        || fail "TEA alignment step died"
fi

# shellcheck disable=SC2086
"$MMSEQS" mvdb "${TMP_PATH}/aln" "${RESULTS}" ${VERBOSITY}

if [ -n "$REMOVE_TMP" ]; then
    echo "Removing temporary files"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
    if [ -n "${RESCORE_PAR}" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_rescore" ${VERBOSITY}
    fi
fi
