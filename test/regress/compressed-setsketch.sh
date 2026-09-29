#!/usr/bin/env bash
# SetSketches generated directly with compressed registers (--full with
# --fastcmp-bytes, --fastcmp-shorts, --fastcmp-words, or --setsketch-ab A,B
# --fastcmp N) must run and estimate Jaccard. Fixtures: a = 200 kbp random DNA,
# b = a[:100000] + 100 kbp of new sequence (exact Jaccard of canonical 21-mers
# computed in Python, about 1/3), and a2 = identical copy of a. A second run
# with --cache must reload the saved sketches and give the same value.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
J=$(python3 - <<'EOF'
import random
R = random.Random(9)
rnd = lambda n: ''.join(R.choice('ACGT') for _ in range(n))
a = rnd(200000)
b = a[:100000] + rnd(100000)
for name, s in (('a', a), ('b', b), ('a2', a)):
    open(name + '.fa', 'w').write('>' + name + '\n' + s + '\n')
def kmers(x, k=21):
    rc = x.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    n = len(x)
    return {min(x[i:i+k], rc[n-i-k:n-i]) for i in range(n - k + 1)}
A, B = kmers(a), kmers(b)
print(len(A & B) / len(A | B))
EOF
)
fail=0
check() { # name got expected tolerance
    if [ -n "$2" ] && python3 -c "import sys; sys.exit(abs($2 - $3) > $4)"; then
        echo "PASS $1: got $2 expected $3 (tolerance $4)"
    else
        echo "FAIL $1: got '$2' expected $3 (tolerance $4)"; fail=1
    fi
}
for opt in --fastcmp-bytes --fastcmp-shorts --fastcmp-words "--setsketch-ab 20,1.2 --fastcmp 1"; do
    out=$("$D2" cmp -k 21 --full $opt a.fa b.fa a2.fa 2>/dev/null)
    name="${opt// /_}"
    check "related$name" "$(echo "$out" | awk '$1 == "a.fa" {print $3}')" "$J" 0.05
    check "identical$name" "$(echo "$out" | awk '$1 == "a.fa" {print $4}')" 1 0
done
first=$("$D2" cmp -k 21 --full --fastcmp-bytes --cache a.fa b.fa 2>/dev/null | awk '$1 == "a.fa" {print $3}')
second=$("$D2" cmp -k 21 --full --fastcmp-bytes --cache a.fa b.fa 2>/dev/null | awk '$1 == "a.fa" {print $3}')
check cache_reload "$second" "${first:-nan}" 0
rm -f ./*.kmerset* ./*.ss ./*.opss
exit $fail
