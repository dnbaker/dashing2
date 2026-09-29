#!/usr/bin/env bash
# `contain` must find the k-mers of a database sketched with -2 (--long-kmers).
# Database: x.fa and y.fa (unrelated random 20 kbp). Invariants: x.fa against
# its own sketch is 100%, y.fa against x's sketch is 0%, and h.fa (the first
# half of x.fa) covers about half of x's sampled k-mers (1024 samples, so the
# standard error is about 1.6%; accept 40 to 60%). Checked on the direct
# encoder (k 21, k 40) and the rolling hasher (k 70), with -p 1 and -p 2. The
# k 40 case uses --no-canon so that it does not depend on 128-bit
# canonicalization in bonsai.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
x = ''.join(r.choice('ACGT') for _ in range(20000))
y = ''.join(r.choice('ACGT') for _ in range(20000))
open('x.fa', 'w').write('>x\n' + x + '\n')
open('h.fa', 'w').write('>h\n' + x[:10000] + '\n')
open('y.fa', 'w').write('>y\n' + y + '\n')
PY
fail=0
check() { # name got lo hi
    if awk -v g="$2" -v l="$3" -v h="$4" 'BEGIN{exit !(g >= l && g <= h)}'; then
        echo "PASS $1: got $2% expected $3..$4%"
    else
        echo "FAIL $1: got $2% expected $3..$4%"; fail=1
    fi
}
for k in 21 40 70; do
    nc=; [ $k = 40 ] && nc=--no-canon
    "$D" sketch --full -2 -k $k $nc --save-kmers -o db x.fa y.fa 2>/dev/null
    for p in 1 2; do
        out=$("$D" contain -p $p db.kmer64 x.fa h.fa y.fa 2>/dev/null)
        cov() { echo "$out" | awk -v q="$1" '$1==q{split($2,a,"%"); print a[1]}'; }
        check "-2 k=$k -p $p x.fa in x" "$(cov x.fa)" 100 100
        check "-2 k=$k -p $p h.fa in x" "$(cov h.fa)" 40 60
        check "-2 k=$k -p $p y.fa in x" "$(cov y.fa)" 0 0
    done
    rm -f db db.*
done
exit $fail
