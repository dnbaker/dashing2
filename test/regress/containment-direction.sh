#!/usr/bin/env bash
# --containment is documented as |A & B| / |A| with A the row entity. Every
# sketch type must use the same direction. Fixtures: a = 400 kbp random DNA,
# s = a[100000:200000], so row a / column s is about 0.25 and row s / column a
# is 1.
# The exact values over canonical 21-mers are computed in Python.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'EOF' > exact.txt
import random
R = random.Random(5)
a = ''.join(R.choice('ACGT') for _ in range(400000))
s = a[100000:200000]
open('a.fa', 'w').write('>a\n' + a + '\n')
open('s.fa', 'w').write('>s\n' + s + '\n')
def kmers(x, k=21):
    rc = x.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    n = len(x)
    return {min(x[i:i+k], rc[n-i-k:n-i]) for i in range(n - k + 1)}
A, S = kmers(a), kmers(s)
print(len(A & S) / len(A), len(A & S) / len(S))
EOF
read CAS CSA < exact.txt
fail=0
check() { # name got expected; sketches are estimates, so allow 10% relative error
    if python3 -c "import sys; sys.exit(abs($2 - $3) > 0.1 * $3)"; then
        echo "PASS $1: got $2 expected $3"
    else
        echo "FAIL $1: got $2 expected $3"; fail=1
    fi
}
for mode in "" --full --set; do
    out=$("$D2" cmp -k 21 -S 4096 --containment --square $mode a.fa s.fa 2>/dev/null)
    check "row_a_col_s${mode:+_$mode}" "$(echo "$out" | awk '$1 == "a.fa" {print $3}')" "$CAS"
    check "row_s_col_a${mode:+_$mode}" "$(echo "$out" | awk '$1 == "s.fa" {print $2}')" "$CSA"
done
rm -f ./*.kmerset* ./*.opss
exit $fail
