#!/usr/bin/env bash
# A joint -F entry "a.fa b.fa" is the union of both files. It must not share a
# cache or k-mer file with a.fa. Checks: joint~ab.fa (ab.fa is a.fa and b.fa
# concatenated) must be 1; joint~a.fa must equal the exact Jaccard index of the
# canonical 21-mer sets, computed below.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
expected=$(python3 - <<'PY'
import random
r = random.Random(1)
a = ''.join(r.choice('ACGT') for _ in range(5000))
b = ''.join(r.choice('ACGT') for _ in range(5000))
open('a.fa', 'w').write('>a\n' + a + '\n')
open('b.fa', 'w').write('>b\n' + b + '\n')
open('ab.fa', 'w').write('>a\n' + a + '\n>b\n' + b + '\n')
open('list', 'w').write('a.fa b.fa\na.fa\nab.fa\n')
rc = str.maketrans('ACGT', 'TGCA')
def kmers(s, k=21):
    return {min(s[i:i+k], s[i:i+k].translate(rc)[::-1]) for i in range(len(s) - k + 1)}
A, U = kmers(a), kmers(a) | kmers(b)
print(len(A & U) / len(A | U))
PY
)
fail=0
for mode in --set "--full --cache"; do
    for run in 1 2; do   # the second run exercises files left by the first
        out=$("$D" cmp $mode -k 21 -F list 2>/dev/null)
        # Row for the joint entry: columns are joint, a.fa, ab.fa
        read -r j_a j_ab <<<"$(echo "$out" | awk 'NR==4{print $(NF-1), $NF}')"
        ok_a=$(python3 -c "print(abs(float('${j_a:-nan}') - $expected) < 0.02)")
        if [ "$ok_a" = True ]; then echo "PASS $mode run $run joint~a: got $j_a expected $expected"
        else echo "FAIL $mode run $run joint~a: got ${j_a:-none} expected $expected"; fail=1; fi
        if [ "${j_ab:-none}" = 1 ]; then echo "PASS $mode run $run joint~ab: got $j_ab expected 1"
        else echo "FAIL $mode run $run joint~ab: got ${j_ab:-none} expected 1"; fail=1; fi
    done
    rm -f a.fa.* b.fa.* ab.fa.*
done
exit $fail
