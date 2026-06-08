#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

if [ -e "${IN}.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${LIST}" "${IN}" "${OUT}" ${CREATESUBDB1_PAR} \
        || fail "createsubdb died"
    sort -nk2 "${OUT}.index" > "${OUT}.indextmp"
fi

if [ -e "${IN}_aa.dbtype" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" base:createsubdb "${OUT}.indextmp" "${IN}_aa" "${OUT}_aa" ${CREATESUBDB2_PAR} \
        || fail "createsubdb on _aa died"
fi

if [ -e "${OUT}.sh" ]; then
    rm -f -- "${OUT}.sh"
    rm -f -- "${OUT}.indextmp"
fi
