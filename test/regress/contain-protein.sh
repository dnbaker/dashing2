#!/usr/bin/env bash
# `contain` must parse queries with the database's alphabet. Database: p.fa and
# q.fa (unrelated random 20,000-residue proteins). Invariants: p.fa against its
# own sketch is 100%, q.fa is 0%, and h.fa (the first half of p.fa) covers about
# half of p's 1024 sampled k-mers (standard error about 1.6%; accept 40 to 60%).
# The k values cover the direct encoder, the rolling hasher (k above the 64-bit
# capacity of the alphabet) and -2.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
aa = 'ACDEFGHIKLMNPQRSTVWY'
p = ''.join(r.choice(aa) for _ in range(20000))
q = ''.join(r.choice(aa) for _ in range(20000))
open('p.fa', 'w').write('>p\n' + p + '\n')
open('h.fa', 'w').write('>h\n' + p[:10000] + '\n')
open('q.fa', 'w').write('>q\n' + q + '\n')
PY
fail=0
check() { # name got lo hi
    if awk -v g="$2" -v l="$3" -v h="$4" 'BEGIN{exit !(g >= l && g <= h)}'; then
        echo "PASS $1: got $2% expected $3..$4%"
    else
        echo "FAIL $1: got $2% expected $3..$4%"; fail=1
    fi
}
for opts in "--protein -k 8" "--protein -k 20" "--protein -2 -k 20" "--protein14 -k 10" "--protein8 -k 12" "--protein6 -k 14"; do
    "$D" sketch --full $opts --save-kmers -o db p.fa q.fa 2>/dev/null
    out=$("$D" contain db.kmer64 p.fa h.fa q.fa 2>/dev/null)
    cov() { echo "$out" | awk -v q="$1" '$1==q{split($2,a,"%"); print a[1]}'; }
    check "$opts p.fa in p" "$(cov p.fa)" 100 100
    check "$opts h.fa in p" "$(cov h.fa)" 40 60
    check "$opts q.fa in p" "$(cov q.fa)" 0 0
    rm -f db db.*
done
exit $fail
