#!/usr/bin/env bash
# --set --cache when only some inputs are cached: the result must match an
# uncached run of the same comparison.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
a = ''.join(r.choice('ACGT') for _ in range(5000))
b = list(a)
for i in range(0, len(b), 97):
    b[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[b[i]]
open('a.fa', 'w').write('>a\n' + a + '\n')
open('b.fa', 'w').write('>b\n' + ''.join(b) + '\n')
PY
val() { "$D" cmp "$@" a.fa b.fa 2>/dev/null | awk 'NR==4{print $3}'; }
fail=0
for mode in --set; do
    fresh=$(val $mode -k 21)
    rm -f a.fa.* b.fa.*
    "$D" cmp $mode -k 21 --cache a.fa >/dev/null 2>&1   # cache a.fa only
    partial=$(val $mode -k 21 --cache)
    rm -f a.fa.* b.fa.*
    if [ -n "$fresh" ] && [ "$partial" = "$fresh" ]; then
        echo "PASS $mode with only a.fa cached: got $partial expected $fresh"
    else
        echo "FAIL $mode with only a.fa cached: got '$partial' expected $fresh"; fail=1
    fi
done
exit $fail
