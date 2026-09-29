#!/usr/bin/env bash
# `cmp --presketched` on k-mer sequence files (.mmerseq64/.mmerseq128 from
# `sketch -G`) reports the number of positions whose k-mers are equal. Positions
# past the end of the shorter sequence are not matches, the result must not
# depend on argument order, and 128-bit k-mers must be compared as 128-bit items.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
# a: 5 kbp random; b: first 3 kbp of a plus 1 kbp new sequence.
python3 - <<'PY'
import random
r = random.Random(3)
a = ''.join(r.choice('ACGT') for _ in range(5000))
b = a[:3000] + ''.join(r.choice('ACGT') for _ in range(1000))
open('a.fa', 'w').write(f'>a\n{a}\n')
open('b.fa', 'w').write(f'>b\n{b}\n')
PY
# Oracle: count equal items over the common prefix of the two emitted streams.
equal_positions() {
    python3 - "$1" "$2" "$3" <<'PY'
import sys
w = int(sys.argv[3])
def items(fn):
    d = open(fn, 'rb').read()
    return [d[i:i + w] for i in range(0, len(d), w)]
print(sum(x == y for x, y in zip(items(sys.argv[1]), items(sys.argv[2]))))
PY
}
check() { # name got expected
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
for bits in 64 128; do
    flag=""; [ $bits = 128 ] && flag="-2"
    "$D2" sketch $flag -G -k 21 -o stacked a.fa b.fa 2>/dev/null
    A=$(ls a.fa*.mmerseq$bits); B=$(ls b.fa*.mmerseq$bits)
    cp "$A" acopy.mmerseq$bits
    exp=$(equal_positions "$A" "$B" $((bits / 8)))
    ab=$("$D2" cmp --presketched "$A" "$B" 2>/dev/null | grep -v '^#' | awk 'NR == 1 {print $3}')
    ba=$("$D2" cmp --presketched "$B" "$A" 2>/dev/null | grep -v '^#' | awk 'NR == 1 {print $3}')
    self=$("$D2" cmp --presketched "$A" acopy.mmerseq$bits 2>/dev/null | grep -v '^#' | awk 'NR == 1 {print $3}')
    check "mmerseq$bits (a,b)" "$ab" "$exp"
    check "mmerseq$bits (b,a)" "$ba" "$exp"
    # A file against an identical copy matches at every position.
    check "mmerseq$bits (a,copy of a)" "$self" "$(( $(wc -c < "$A") / (bits / 8) ))"
    rm -f ./*.mmerseq* stacked*
done
exit $fail
