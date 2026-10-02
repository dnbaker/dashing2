#!/usr/bin/env bash
# With --edit-distance (OrderMinHash), a record and an identical copy of it
# have identical sketches, so their similarity must be exactly 1. The
# registers are 64-bit hashes stored as doubles; about 1 in 2048 is a NaN bit
# pattern, and a floating-point comparison counts it as unequal. Records are
# 100 bp, and four of them with -S 4096 make a NaN register all but certain.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random
r = random.Random(3)
with open('dup.fa', 'w') as f:
    for i in range(4):
        s = ''.join(r.choice('ACGT') for _ in range(100))
        f.write(f'>s{i}\n{s}\n>s{i}copy\n{s}\n')
PY
"$D2" sketch -p1 -S4096 -k7 --parse-by-seq --edit-distance --cmpout out.tsv dup.fa 2>/dev/null; rc=$?
for i in 0 1 2 3; do
    got=$(awk -F'\t' -v a="s$i" -v b="s${i}copy" '/^#Sources/ {for (c = 2; c <= NF; c++) if ($c == b) col = c; next}
                                                /^#/ {next} {sub(/ +$/, "", $1)} $1 == a {print $col}' out.tsv 2>/dev/null)
    if [ $rc -eq 0 ] && [ "$got" = 1 ]; then echo "PASS s$i vs identical copy: got $got expected 1"
    else echo "FAIL s$i vs identical copy: status $rc, got '$got' expected 1"; fail=1; fi
done
exit $fail
