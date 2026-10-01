#!/usr/bin/env bash
# Directly sketched compressed SetSketches (--full with --fastcmp-bytes,
# --fastcmp-shorts or --fastcmp-words) must not depend on which thread sketches
# which file. Each thread reuses one sketch object for the files it processes,
# so sketching a file after another one must give the same sketch as sketching
# it alone, and with a fixed seed -p 8 must print exactly what -p 1 prints.
# Fixtures (Python, random DNA): big = 200 kbp, small = 2 kbp; m1 = a 20 kbp
# sequence s plus a second record s[:5000], m2 = s alone (with -m 2, m2 has no
# k-mer seen twice, so its cardinality must stay at its value alone however it is
# ordered after m1); g00..g15 = 16 sequences of 2 kbp.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
REPS="${REPS:-10}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
R = random.Random(8)
rnd = lambda n: ''.join(R.choice('ACGT') for _ in range(n))
def w(name, *seqs):
    open(name, 'w').write(''.join('>r%d\n%s\n' % (i, s) for i, s in enumerate(seqs)))
w('big.fa', rnd(200000)); w('small.fa', rnd(2000))
s = rnd(20000); w('m1.fa', s, s[:5000]); w('m2.fa', s)
for i in range(16):
    w('g%02d.fa' % i, rnd(2000))
PY
fail=0
check() { # name got expected
    if [ -n "$2" ] && [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got '$2' expected '$3'"; fail=1; fi
}
# Cardinality of the last input, read from the diagonal of --union-size --square.
card() { "$D2" cmp -k 21 -S 1024 --full -p 1 --union-size --square "$@" 2>/dev/null | awk '!/^#/ {v = $NF} END {print v}'; }
for w in bytes shorts words; do
    check "small_after_big_$w" "$(card --fastcmp-$w big.fa small.fa)" "$(card --fastcmp-$w small.fa)"
    check "m2_after_m1_count2_$w" "$(card --fastcmp-$w -m 2 m1.fa m2.fa)" "$(card --fastcmp-$w -m 2 m2.fa)"
    ref=$("$D2" cmp -k 21 -S 1024 --full --fastcmp-$w -p 1 --union-size --square g*.fa 2>/dev/null | grep -v '^#' | cksum)
    bad=0
    for i in $(seq "$REPS"); do
        out=$("$D2" cmp -k 21 -S 1024 --full --fastcmp-$w -p 8 --union-size --square g*.fa 2>/dev/null | grep -v '^#' | cksum)
        [ "$out" = "$ref" ] || bad=$((bad + 1))
    done
    check "p8_equals_p1_$w" "$bad of $REPS runs differ" "0 of $REPS runs differ"
done
rm -f ./*.kmerset* ./*.ss ./*.opss
exit $fail
