#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

# Create query TEA DB if input is FASTA pair
if notExists "${TMP_PATH}/query.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "${QUERY_TEA}" "${QUERY_AA}" "${TMP_PATH}/query" ${CREATEDB_PAR} \
        || fail "query createteadb died"
fi
QUERY="${TMP_PATH}/query"

# Create target TEA DB if needed
if notExists "${TARGET}.dbtype"; then
    if notExists "${TMP_PATH}/target.dbtype"; then
        # shellcheck disable=SC2086
        "$MMSEQS" createdb "${TARGET_TEA}" "${TARGET_AA}" "${TMP_PATH}/target" ${CREATEDB_PAR} \
            || fail "target createteadb died"
    fi
    TARGET="${TMP_PATH}/target"
fi

INTERMEDIATE="${TMP_PATH}/result"
if notExists "${INTERMEDIATE}.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" search "${QUERY}" "${TARGET}" "${INTERMEDIATE}" "${TMP_PATH}/search_tmp" ${SEARCH_PAR} \
        || fail "Search died"
fi

if notExists "${TMP_PATH}/alis.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" convertalis "${QUERY}" "${TARGET}${INDEXEXT}" "${INTERMEDIATE}" "${RESULTS}" ${CONVERT_PAR} \
        || fail "Convert Alignments died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    if [ -f "${TMP_PATH}/target.dbtype" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/target" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/target_h" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/target_aa" ${VERBOSITY}
    fi
    if [ -f "${TMP_PATH}/query.dbtype" ]; then
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_h" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/query_aa" ${VERBOSITY}
    fi
    rm -rf "${TMP_PATH}/search_tmp"
    rm -f "${TMP_PATH}/easyteasearch.sh"
fi
