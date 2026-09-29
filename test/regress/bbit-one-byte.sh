#!/usr/bin/env bash
# --bbit-sigs --fastcmp 1 keeps one byte per register and corrects the match
# rate for the 1/256 chance that two different registers agree on that byte.
# Two unrelated random sequences share no 21-mers, so the corrected similarity
# must be about 0 (the standard error at -S 16384 is about 0.0005), and a file
# compared with an identical copy must give 1.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'EOF'
import random
R = random.Random(13)
for name in ('a', 'd'):
    open(name + '.fa', 'w').write('>' + name + '\n' + ''.join(R.choice('ACGT') for _ in range(200000)) + '\n')
open('a2.fa', 'w').write(open('a.fa').read())
EOF
fail=0
check() { # name got expected tolerance
    if python3 -c "import sys; sys.exit(abs($2 - $3) > $4)"; then
        echo "PASS $1: got $2 expected $3 (tolerance $4)"
    else
        echo "FAIL $1: got $2 expected $3 (tolerance $4)"; fail=1
    fi
}
out=$("$D2" cmp -k 21 -S 16384 --bbit-sigs --fastcmp 1 a.fa d.fa a2.fa 2>/dev/null)
check unrelated "$(echo "$out" | awk '$1 == "a.fa" {print $3}')" 0 0.003
check identical "$(echo "$out" | awk '$1 == "a.fa" {print $4}')" 1 0
rm -f ./*.kmerset* ./*.opss
exit $fail
