#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

# Build a synthetic prefilter result that pairs every query with every target.
# Borrowed from mmseqs's blastp.sh; the "data" file is a symlink to the
# target's .index so every query row references the same all-targets payload.
fake_pref() {
    QDB="$1"
    TDB="$2"
    RES="$3"
    ln -s "$(abspath "${TDB}.index")" "${RES}"
    INDEX_SIZE="$(wc -c < "${TDB}.index")"
    awk -v size="$INDEX_SIZE" '{ print $1"\t0\t"size; }' "${QDB}.index" > "${RES}.index"
    awk 'BEGIN { printf("%c%c%c%c",7,0,0,0); exit; }' > "${RES}.dbtype"
}

# 1. Prefilter (k-mer matching on TEA sequences) or true exhaustive all-pairs
if notExists "${TMP_PATH}/pref.dbtype"; then
    if [ -n "${EXHAUSTIVE}" ]; then
        # True all-vs-all: synthetic prefilter DB pairing every (query, target).
        fake_pref "${QUERY}" "${TARGET}${INDEXEXT}" "${TMP_PATH}/pref" \
            || fail "Exhaustive fake_pref step died"
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
