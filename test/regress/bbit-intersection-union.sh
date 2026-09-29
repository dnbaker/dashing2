#!/usr/bin/env bash
# With --bbit-sigs and compressed registers (--fastcmp 1, 2 or 4),
# --intersection must estimate |A & B| and --union-size |A | B|.
# Fixtures: a = 200 kbp random DNA, b = a[:100000] + 100 kbp of new sequence,
# so |A & B| is about 100k and |A | B| about 300k. Exact values over canonical
# 21-mers are computed in Python.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'EOF' > exact.txt
import random
R = random.Random(9)
rnd = lambda n: ''.join(R.choice('ACGT') for _ in range(n))
a = rnd(200000)
b = a[:100000] + rnd(100000)
open('a.fa', 'w').write('>a\n' + a + '\n')
open('b.fa', 'w').write('>b\n' + b + '\n')
def kmers(x, k=21):
    rc = x.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    n = len(x)
    return {min(x[i:i+k], rc[n-i-k:n-i]) for i in range(n - k + 1)}
A, B = kmers(a), kmers(b)
print(len(A & B), len(A | B))
EOF
read I U < exact.txt
fail=0
check() { # name got expected; sketches are estimates, so allow 10% relative error
    if python3 -c "import sys; sys.exit(abs($2 - $3) > 0.1 * $3)"; then
        echo "PASS $1: got $2 expected $3"
    else
        echo "FAIL $1: got $2 expected $3"; fail=1
    fi
}
val() { "$D2" cmp -k 21 --bbit-sigs "$@" a.fa b.fa 2>/dev/null | awk '$1 == "a.fa" {print $3}'; }
for n in 1 2 4; do
    check "intersection_fastcmp$n" "$(val --fastcmp $n --intersection)" "$I"
    check "union_fastcmp$n" "$(val --fastcmp $n --union-size)" "$U"
done
rm -f ./*.kmerset* ./*.opss
exit $fail
