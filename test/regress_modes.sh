#!/usr/bin/env bash
# dashing2 regression: fixed sketch + compare modes on test/data/seqs.fa
# (100 sequences), outputs compared byte for byte with test/golden/. The
# sketch release gate (sketch scripts/release_check.sh) runs this.
#   test/regress_modes.sh [./dashing2] [--update]
# --update rewrites the goldens: only for an intended change in results, and
# say why in the commit.
set -euo pipefail
cd "$(dirname "$0")/.."
D2=${1:-./dashing2}
UPDATE=${2:-}
T=$(mktemp -d "${TMPDIR:-/tmp}/d2_regress_XXXXXX")
trap 'rm -rf "$T"' EXIT
i=0; bad=0
while read -r args; do
    i=$((i + 1))
    "$D2" sketch --parse-by-seq -k 9 -S 128 -p 4 $args --cmpout "$T/m$i.txt" test/data/seqs.fa > "$T/log$i.txt" 2>&1 \
        || { cat "$T/log$i.txt"; echo "FAIL: mode $i ($args) exited nonzero" >&2; exit 1; }
    if [ "$UPDATE" = --update ]; then
        cp "$T/m$i.txt" "test/golden/m$i.txt"
    elif ! cmp -s "$T/m$i.txt" "test/golden/m$i.txt"; then
        echo "FAIL: mode $i ($args) differs from test/golden/m$i.txt" >&2; bad=1
    fi
done <<'MODES'
--set
-B
--prob
--prob --countsketch-size 256
-B --countsketch-size 256
--oph
MODES
[ $bad -eq 0 ] || exit 1
if [ "$UPDATE" = --update ]; then echo "regress_modes: $i goldens updated"; else echo "regress_modes: $i modes match the goldens"; fi
