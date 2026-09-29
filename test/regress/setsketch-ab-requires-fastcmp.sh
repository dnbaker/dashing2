#!/usr/bin/env bash
# --setsketch-ab A,B parameterizes compressed SetSketch registers, so it needs a
# register size of 1, 2 or 4 bytes. Without --fastcmp it must be rejected up
# front with a message naming the option, rather than failing after sketching.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random
r = random.Random(1)
a = ''.join(r.choice('ACGT') for _ in range(3000))
b = a[:1500] + ''.join(r.choice('ACGT') for _ in range(1500))
open('a.fa', 'w').write(f'>a\n{a}\n')
open('b.fa', 'w').write(f'>b\n{b}\n')
PY
for sub in "cmp" "sketch -o out"; do
    msg=$("$D2" $sub -k 21 --full --setsketch-ab 20,1.2 a.fa b.fa 2>&1 >/dev/null)
    rc=$?
    rm -f ./*.ss* out*
    if [ $rc -ne 0 ] && grep -q -- '--setsketch-ab requires --fastcmp' <<< "$msg"; then
        echo "PASS $sub --full --setsketch-ab without --fastcmp: got rejection (exit $rc) expected rejection"
    else
        echo "FAIL $sub --full --setsketch-ab without --fastcmp: got exit $rc, '$(tail -1 <<< "$msg")' expected '--setsketch-ab requires --fastcmp'"; fail=1
    fi
done
# With a register size the option is accepted (by-seq mode, which sketches compressed registers directly).
cat a.fa b.fa > ab.fa
got=$("$D2" cmp -k 21 --parse-by-seq --full --setsketch-ab 20,1.2 --fastcmp 1 ab.fa 2>/dev/null | grep -v '^#' | awk 'NR == 1 {print $3}')
if [ -n "$got" ]; then echo "PASS --setsketch-ab 20,1.2 --fastcmp 1 accepted: got $got expected a similarity"
else echo "FAIL --setsketch-ab 20,1.2 --fastcmp 1 accepted: got nothing expected a similarity"; fail=1; fi
exit $fail
