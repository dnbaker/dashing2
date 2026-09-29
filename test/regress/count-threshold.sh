#!/usr/bin/env bash
# -m/--count-threshold N keeps k-mers with count >= N (see --help). Checks the
# exact modes (--set, --countdict) against a Python oracle, -B (BagMinHash)
# within sampling error, and that a threshold leaving no k-mers does not crash.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
check() { # name got expected tolerance
    if python3 -c "import sys; sys.exit(not abs(float('$2') - float('$3')) <= $4)" 2>/dev/null; then
        echo "PASS $1: got $2 expected $3"
    else
        echo "FAIL $1: got $2 expected $3"; fail=1
    fi
}
# Four random 1.5 kbp units; each file repeats units a given number of times,
# so every k-mer of unit j has count reps[j] in that file. Inputs are listed in
# decreasing size order, the order in which they are sketched.
python3 - <<'PY'
import random
r = random.Random(3)
U = [''.join(r.choice('ACGT') for _ in range(1500)) for _ in range(4)]
for name, reps in [('m3', [10, 2, 4, 0]), ('m1', [5, 1, 2, 0]), ('m2', [1, 4, 0, 2]), ('u1', [1, 0, 0, 0]), ('u2', [0, 1, 0, 0])]:
    with open(name + '.fa', 'w') as f:
        f.writelines(f'>{j}.{i}\n{U[j]}\n' for j, c in enumerate(reps) for i in range(c))
PY
# Exact oracle over canonical 21-mers: prints "m pair setJ intersection weightedJ".
python3 - >oracle.txt <<'PY'
from collections import Counter
comp = str.maketrans('ACGT', 'TGCA')
def counts(fn):
    c = Counter()
    for line in open(fn):
        if line.startswith('>'): continue
        x = line.strip()
        for i in range(len(x) - 20):
            k = x[i:i + 21]; c[min(k, k.translate(comp)[::-1])] += 1
    return c
C = {n: counts(n + '.fa') for n in ('m3', 'm1', 'm2')}
for m in (1, 2, 3):
    F = {n: Counter({k: v for k, v in c.items() if v >= m}) for n, c in C.items()}
    for a, b in (('m3', 'm1'), ('m3', 'm2'), ('m1', 'm2')):
        A, B = F[a], F[b]; ks = set(A) | set(B); i = len(set(A) & set(B))
        wj = sum(min(A[k], B[k]) for k in ks) / sum(max(A[k], B[k]) for k in ks)
        print(m, f'{a}~{b}', i / len(ks), i, wj)
PY
# Upper-triangle entry for pair "row~col" from a 3x3 matrix in input order m3 m1 m2.
entry() { awk -v p="$2" 'BEGIN {split("m3~m1 1 3|m3~m2 1 4|m1~m2 2 4", t, "|"); for (x in t) {split(t[x], f, " "); R[f[1]] = f[2]; K[f[1]] = f[3]}}
    !/^#/ {n++; if (n == R[p]) print $K[p]}' "$1"; }
run() { "$D2" cmp -k 21 "$@" m3.fa m1.fa m2.fa 2>/dev/null; rm -f ./*.kmer*; }
for m in 1 2 3; do
    run --set -m $m >set.txt
    run --set --intersection -m $m >int.txt
    run --countdict -m $m >cd.txt
    run -B -S 4096 -m $m >bmh.txt
    while read -r om pair sj inter wj; do
        [ "$om" = "$m" ] || continue
        check "-m $m --set J $pair" "$(entry set.txt "$pair")" "$sj" 1e-4
        check "-m $m --set intersection $pair" "$(entry int.txt "$pair")" "$inter" 0
        check "-m $m --countdict weighted J $pair" "$(entry cd.txt "$pair")" "$wj" 1e-4
        check "-m $m -B weighted J $pair (+-0.03)" "$(entry bmh.txt "$pair")" "$wj" 0.03
    done <oracle.txt
done
# Every k-mer occurs once in u1 and u2, so -m 2 leaves both sets empty.
"$D2" cmp -k 21 --set --intersection -m 2 u1.fa u2.fa >/dev/null 2>&1; rc=$?
rm -f ./*.kmer*
check "--set -m 2 with no qualifying k-mers exit status" "$rc" 0 0
exit $fail
