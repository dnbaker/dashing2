#!/usr/bin/env bash
# The default one-permutation SetSketch must estimate Jaccard correctly when
# the sketches being compared have different numbers of empty buckets (inputs
# with fewer distinct k-mers than a few times the sketch size).
# Fixtures: x = 5 kbp random DNA, y = x with a substitution every 97 bp,
# h = the first half of x, x2 = identical copy of x. The exact Jaccard of
# canonical 21-mers is computed in Python and compared with dashing2 cmp.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'EOF' > exact.txt
import random
R = random.Random(7)
x = ''.join(R.choice('ACGT') for _ in range(5000))
y = list(x)
for i in range(0, len(y), 97): y[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[y[i]]
y = ''.join(y)
h = x[:2500]
for name, s in (('x', x), ('y', y), ('h', h), ('x2', x)):
    open(name + '.fa', 'w').write('>' + name + '\n' + s + '\n')
def kmers(s, k=21):
    rc = s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    n = len(s)
    return {min(s[i:i+k], rc[n-i-k:n-i]) for i in range(n - k + 1)}
X, Y, H = kmers(x), kmers(y), kmers(h)
print(len(X & Y) / len(X | Y), len(X & H) / len(X | H))
EOF
read JXY JXH < exact.txt
fail=0
check() { # name got expected tolerance
    if python3 -c "import sys; sys.exit(abs($2 - $3) > $4)"; then
        echo "PASS $1: got $2 expected $3 (tolerance $4)"
    else
        echo "FAIL $1: got $2 expected $3 (tolerance $4)"; fail=1
    fi
}
sim() { # sketch size, file1, file2
    "$D2" cmp -k 21 -S "$1" "$2" "$3" 2>/dev/null | awk -v f="$2" '$1 == f {print $3}'
}
for S in 1024 2048 4096 8192; do
    check "x_vs_y_S$S" "$(sim $S x.fa y.fa)" "$JXY" 0.05
    check "x_vs_half_S$S" "$(sim $S x.fa h.fa)" "$JXH" 0.05
    check "x_vs_copy_S$S" "$(sim $S x.fa x2.fa)" 1 0
done
rm -f ./*.kmerset* ./*.opss
exit $fail
