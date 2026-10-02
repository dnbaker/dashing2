#!/usr/bin/env bash
# contain must report each reference's coverage and depth in that reference's
# column. x86 builds format the columns in blocks of 8 (AVX2) or 16 (AVX-512)
# and the rest one at a time, so this uses 21 references to reach both.
# Fixtures (Python, random DNA): 21 references of 2000 bp; the query holds the
# first 300 + 80 * i bp of reference i, repeated i % 3 + 1 times, so every
# column has its own coverage and depth. Expected value for column i: contain
# against a database holding only reference i (one column, so no block path).
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(5)
refs = [''.join(r.choice('ACGT') for _ in range(2000)) for _ in range(21)]
for i, s in enumerate(refs):
    open('r%02d.fa' % i, 'w').write('>r%02d\n%s\n' % (i, s))
open('q.fa', 'w').write(''.join('>q%d\n%s\n' % (i, s[:300 + 80 * i] * (i % 3 + 1)) for i, s in enumerate(refs)))
PY
"$D" sketch --full -k 21 --save-kmers -o db r*.fa >/dev/null 2>&1
fail=0
got=($("$D" contain db.kmer64 q.fa 2>/dev/null | awk -F'\t' '$1=="q.fa"{for(i=2;i<=NF;i++) print $i}'))
for i in $(seq 0 20); do
    f=$(printf 'r%02d.fa' $i)
    "$D" sketch --full -k 21 --save-kmers -o one $f >/dev/null 2>&1
    exp=$("$D" contain one.kmer64 q.fa 2>/dev/null | awk -F'\t' '$1=="q.fa"{print $2}')
    if [ "${got[$i]:-missing}" = "$exp" ]; then echo "PASS column $i: got ${got[$i]} expected $exp"
    else echo "FAIL column $i: got ${got[$i]:-missing} expected $exp"; fail=1; fi
done
rm -f ./*.kmerset* ./*.kmer64* db* one*
exit $fail
