#!/usr/bin/env bash
# --topk and --similarity-threshold neighbor lists built with -p 4 must equal
# the -p 1 lists and never name the same neighbor twice in a row. Fixtures: 30
# random 5 kbp genomes in 5 clusters of 6 (member j of cluster c has
# 0.2% * j * (c + 1) substitutions from member 0). The race is timing
# dependent, so each mode is run 10 times.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(7)
sub = {'A': 'CGT', 'C': 'AGT', 'G': 'ACT', 'T': 'ACG'}
for c in range(5):
    base = ''.join(r.choice('ACGT') for _ in range(5000))
    for j in range(6):
        s = list(base)
        for p in r.sample(range(5000), int(5000 * 0.002 * j * (c + 1))):
            s[p] = r.choice(sub[s[p]])
        open('c%d_%d.fa' % (c, j), 'w').write('>c%d_%d\n' % (c, j) + ''.join(s) + '\n')
PY
dups() { awk -F'\t' '!/^#/ && NF>1{delete s; for(i=2;i<=NF;i++){split($i,a,":"); if(a[1] in s) d++; s[a[1]]=1}} END{print d+0}'; }
fail=0
for mode in "--topk 3" "--similarity-threshold 0.1"; do
    "$D" cmp -p 1 --full -k 21 $mode c*.fa 2>/dev/null | grep -v '^#' > ref.txt
    nd=0; ndiff=0
    for i in 1 2 3 4 5 6 7 8 9 10; do
        "$D" cmp -p 4 --full -k 21 $mode c*.fa 2>/dev/null | grep -v '^#' > out.txt
        nd=$((nd + $(dups < out.txt)))
        cmp -s out.txt ref.txt || ndiff=$((ndiff + 1))
    done
    if [ $nd -eq 0 ] && [ $ndiff -eq 0 ]; then
        echo "PASS $mode -p 4: got $nd duplicate neighbors, $ndiff runs differing from -p 1 expected 0, 0"
    else
        echo "FAIL $mode -p 4: got $nd duplicate neighbors, $ndiff runs differing from -p 1 expected 0, 0"; fail=1
    fi
done
exit $fail
