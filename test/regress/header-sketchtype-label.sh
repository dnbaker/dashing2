#!/usr/bin/env bash
# The #Dashing2Options header names the sketch type. Every alias of a sketch type
# must give the label of its short form, and --countdict must be named.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
a = ''.join(r.choice('ACGT') for _ in range(3000))
b = a[:2000] + ''.join(r.choice('ACGT') for _ in range(1000))
open('a.fa', 'w').write('>a\n' + a + '\n')
open('b.fa', 'w').write('>b\n' + b + '\n')
PY
fail=0
label() { "$D" cmp -k 15 -S 256 "$@" a.fa b.fa 2>/dev/null | sed -n 's/.*;sketchtype:\([^;]*\);.*/\1/p'; }
check() {  # check <flag> <expected label>
    local got; got=$(label "$1")
    if [ "$got" = "$2" ]; then echo "PASS $1: got $got expected $2"
    else echo "FAIL $1: got ${got:-nothing} expected $2"; fail=1; fi
}
for f in -B --multiset --bagminhash --bmh --BMH; do check "$f" bagminhash; done
for f in -P --prob --probs --pminhash --pmh --PMH --probminhash; do check "$f" probminhash; done
check -J mmerset64,kmercountsf64
check --countdict mmerset64,kmercountsf64
# The long forms give the same values as the short forms (the label is the only difference).
v() { "$D" cmp -k 15 -S 256 "$@" a.fa b.fa 2>/dev/null | awk 'NR==4{print $3}'; }
for p in "-B --multiset" "-P --prob"; do
    set -- $p
    s=$(v "$1"); l=$(v "$2")
    if [ -n "$s" ] && [ "$s" = "$l" ]; then echo "PASS $2 value: got $l expected $s (as $1)"
    else echo "FAIL $2 value: got ${l:-nothing} expected ${s:-nothing} (as $1)"; fail=1; fi
done
exit $fail
