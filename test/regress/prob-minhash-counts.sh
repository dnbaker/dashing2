#!/usr/bin/env bash
# --prob (ProbMinHash) estimates the probability Jaccard index of k-mer count
# distributions, so it must work when k-mer counts differ between inputs.
# Fixtures are built from four random 1.5 kbp units u0..u3, each copy a
# separate record so every k-mer of a unit has exactly that unit's count:
#   m1 = u0 x5, u1 x1, u2 x2;  m2 = u0 x1, u1 x4, u3 x2;  m3 = 2 x m1.
# m1 and m3 are the same distribution (expected 1). The exact probability
# Jaccard of m1 and m2 is computed in Python (Ertl's J_P formula).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
JP=$(python3 - <<'EOF'
import random
from collections import Counter
R = random.Random(21)
u = [''.join(R.choice('ACGT') for _ in range(1500)) for _ in range(4)]
spec = {'m1': (5, 1, 2, 0), 'm2': (1, 4, 0, 2), 'm3': (10, 2, 4, 0)}
def canon_counts(recs, k=21):
    c = Counter()
    for s in recs:
        rc = s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        for i in range(len(s) - k + 1):
            c[min(s[i:i+k], rc[len(s)-i-k:len(s)-i])] += 1
    return c
C = {}
for name, counts in spec.items():
    recs = [u[j] for j, n in enumerate(counts) for _ in range(n)]
    with open(name + '.fa', 'w') as f:
        for r, s in enumerate(recs): f.write(f'>{name}_{r}\n{s}\n')
    C[name] = canon_counts(recs)
A, B = C['m1'], C['m2']
sa, sb = sum(A.values()), sum(B.values())
pairs = Counter((A.get(x, 0) / sa, B.get(x, 0) / sb) for x in set(A) | set(B))
print(sum(1 / sum(n * max(p / (A[x] / sa), q / (B[x] / sb)) for (p, q), n in pairs.items())
          for x in set(A) & set(B)))
EOF
)
fail=0
check() { # name got expected
    if python3 -c "import sys; sys.exit(abs($2 - $3) > 0.05)"; then
        echo "PASS $1: got $2 expected $3"
    else
        echo "FAIL $1: got $2 expected $3"; fail=1
    fi
}
out=$("$D2" cmp -k 21 --prob m1.fa m2.fa m3.fa 2>/dev/null)
check m1_vs_m2 "$(echo "$out" | awk '$1 == "m1.fa" {print $3}')" "$JP"
check m1_vs_m3 "$(echo "$out" | awk '$1 == "m1.fa" {print $4}')" 1
rm -f ./*.kmerset* ./*.pmh
exit $fail
