#!/usr/bin/env bash
# --greedy T clusters items whose best match is within T. For similarities that
# means similarity >= T; for distances (--mash-distance) it means distance <= T.
# Fixtures: a.fa (random 5 kbp), b.fa (a with 10 substitutions, Mash distance
# about 0.002), c.fa and d.fa unrelated (distance inf). Expected: 3 clusters,
# one of them {a, b}, for the LSH and exhaustive (E) variants, -p 1 and -p 2,
# and with --set as well as --full.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(5)
a = ''.join(r.choice('ACGT') for _ in range(5000))
b = list(a)
for i in range(0, 5000, 500):
    b[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[b[i]]
open('a.fa', 'w').write('>a\n' + a + '\n')
open('b.fa', 'w').write('>b\n' + ''.join(b) + '\n')
for n in 'cd':
    open(n + '.fa', 'w').write('>' + n + '\n' + ''.join(r.choice('ACGT') for _ in range(5000)) + '\n')
PY
fail=0
for mode in --full --set; do
for g in "--mash-distance --greedy 0.05" "--mash-distance --greedy 0.05E" "--greedy 0.5" "--greedy 0.5E"; do
    for p in 1 2; do
        out=$("$D" cmp $mode -k 21 -p $p $g a.fa b.fa c.fa d.fa 2>/dev/null)
        n=$(echo "$out" | grep -c '^Cluster-')
        ab=$(echo "$out" | grep '^Cluster-' | grep 'a.fa' | grep -c 'b.fa')
        if [ "$n" -eq 3 ] && [ "$ab" -eq 1 ]; then
            echo "PASS $mode $g -p $p: got $n clusters, a+b together=$ab expected 3, 1"
        else
            echo "FAIL $mode $g -p $p: got $n clusters, a+b together=$ab expected 3, 1"; fail=1
        fi
        rm -f ./*.kmerset* ./*.kmercounts*
    done
done
done
exit $fail
