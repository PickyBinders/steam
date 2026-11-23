#!/bin/bash
# Prove the ordered patch set recreates the modified dependency sources.
set -euo pipefail
patch_work=$(mktemp -d "${TMPDIR:-/tmp}/steam-patches.XXXXXX")
mapfile -t patch_files < <(awk '/^diff --git / {sub(/^a\//, "", $3); print $3}' patches/*.patch | sort -u)
git -C lib/mmseqs archive HEAD -- "${patch_files[@]}" | tar -x -C "$patch_work"
for patch in patches/*.patch; do
    git -C "$patch_work" apply "$PWD/$patch"
done
for source in "${patch_files[@]}"; do
    cmp "$patch_work/$source" "lib/mmseqs/$source"
done
printf 'The ordered patch series reproduces the dependency sources from pinned HEAD\n'
