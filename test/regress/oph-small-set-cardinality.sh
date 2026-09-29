#!/usr/bin/env bash
# The default one-permutation SetSketch must estimate the number of distinct
# k-mers without bias for sets smaller than a few times the sketch size, and
# must report 0 for an empty input. Fixtures are random DNA sequences; the
# exact number of distinct canonical 21-mers is computed in Python and
# compared with the cardinality dashing2 sketch writes to <out>.names.txt.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'EOF' > exact.txt
import random
R = random.Random(11)
open('empty.fa', 'w').write('')
print('empty.fa', 0)
for L in (120, 520, 1520, 20020):
    s = ''.join(R.choice('ACGT') for _ in range(L))
    open(f'r{L}.fa', 'w').write(f'>r{L}\n{s}\n')
    rc = s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    print(f'r{L}.fa', len({min(s[i:i+21], rc[L-i-21:L-i]) for i in range(L - 20)}))
EOF
"$D2" sketch -k 21 -S 1024 -o out empty.fa r120.fa r520.fa r1520.fa r20020.fa 2>/dev/null
fail=0
while read -r name exact; do
    got=$(awk -v f="$name" '$1 == f {print $2}' out.names.txt)
    # Allow 10% relative error (the sketch is an estimate); the empty set must give 0.
    if python3 -c "import sys; sys.exit(abs($got - $exact) > 0.1 * $exact)"; then
        echo "PASS card_$name: got $got expected $exact"
    else
        echo "FAIL card_$name: got $got expected $exact"; fail=1
    fi
done < exact.txt
rm -f ./*.kmerset* ./*.opss
exit $fail
