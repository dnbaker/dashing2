#!/usr/bin/env bash
# The cardinality of a directly sketched compressed SetSketch (--full with
# --fastcmp-bytes, --fastcmp-shorts or --fastcmp-words, or --setsketch-ab A,B
# --fastcmp N) must be 0 for an input without k-mers and close to the true
# count for small inputs. --fastcmp-shorts uses a = 0.06, so many of its
# registers stay 0 for inputs of fewer than about 100 k-mers. Fixtures (Python):
# empty.fa = 10 bases (shorter than k = 21); n10.fa = 30 random bases, whose
# distinct canonical 21-mers are counted exactly; n20k.fa = 20 kbp, also
# counted exactly. Cardinalities are read from --union-size --square (-S 1024).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
COUNTS=$(python3 - <<'PY'
import random
R = random.Random(4)
rnd = lambda n: ''.join(R.choice('ACGT') for _ in range(n))
def kmers(x, k=21):
    rc = x.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    n = len(x)
    return {min(x[i:i+k], rc[n-i-k:n-i]) for i in range(n - k + 1)}
seqs = {'empty': rnd(10), 'n10': rnd(30), 'n20k': rnd(20000)}
for name, s in seqs.items():
    open(name + '.fa', 'w').write('>' + name + '\n' + s + '\n')
print(len(kmers(seqs['n10'])), len(kmers(seqs['n20k'])))
PY
)
read -r N10 N20K <<< "$COUNTS"
fail=0
check() { # name got expected tolerance
    if [ -n "$2" ] && python3 -c "import sys; sys.exit(abs($2 - $3) > $4)"; then
        echo "PASS $1: got $2 expected $3 (tolerance $4)"
    else
        echo "FAIL $1: got '$2' expected $3 (tolerance $4)"; fail=1
    fi
}
card() { "$D2" cmp -k 21 -S 1024 -p 1 --full --union-size --square "$@" 2>/dev/null | awk '!/^#/ {print $NF; exit}'; }
for opt in --fastcmp-bytes --fastcmp-shorts --fastcmp-words "--setsketch-ab 20,1.2 --fastcmp 1"; do
    name="${opt// /_}"
    check "empty$name" "$(card $opt empty.fa)" 0 0
    # Relative standard error with 1024 registers is about 3%, so 30% is a loose bound.
    check "n10$name" "$(card $opt n10.fa)" "$N10" 3
    check "n20k$name" "$(card $opt n20k.fa)" "$N20K" "$((N20K * 15 / 100))"
done
rm -f ./*.kmerset* ./*.ss ./*.opss
exit $fail
