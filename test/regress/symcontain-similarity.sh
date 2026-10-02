#!/usr/bin/env bash
# --symmetric-containment, |A & B| / min(|A|, |B|), is a similarity in [0, 1]:
# --topk must list the most contained neighbours first, --similarity-threshold T
# must keep pairs with value >= T and --greedy T must join items whose best match
# is >= T. Fixtures: two clusters of three related 4 kbp genomes (a base, a copy
# with 1% substitutions, a 3 kbp prefix with 2% substitutions); the clusters share
# no k-mers. Expected lists are read off the --square matrix of the same run:
# each genome's cluster mates score above 0.6 and the other cluster scores 0.
# Mash distance is a true distance and must still rank ascending.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
def mut(s, rate):
    return ''.join(c if r.random() > rate else r.choice('ACGT'.replace(c, '')) for c in s)
for c in range(2):
    a = ''.join(r.choice('ACGT') for _ in range(4000))
    for i, s in enumerate([a, mut(a, 0.01), mut(a[:3000], 0.02)]):
        open(f'c{c}_{i}.fa', 'w').write(f'>c{c}_{i}\n{s}\n')
PY
F="c0_0.fa c0_1.fa c0_2.fa c1_0.fa c1_1.fa c1_2.fa"
O="-k 17 --full"
"$D" cmp $O --symmetric-containment --square $F > square.txt 2>/dev/null
fail=0
check() { # name, expected, got
    if [ "$2" = "$3" ]; then echo "PASS $1: got $3 expected $2"; else echo "FAIL $1: got $3 expected $2"; fail=1; fi
}
# Neighbour lists as "query:n1,n2;..." in printed order, and the same lists
# derived from the square matrix (sorted by descending value for similarities).
lists() { python3 -c '
import sys
out = []
for l in open(sys.argv[1]):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    out.append(f[0] + ":" + ",".join(x.rsplit(":", 1)[0] for x in f[1:]))
print(";".join(out))' "$1"; }
expect() { python3 -c '
import sys
mode, t = sys.argv[2], float(sys.argv[3])
rows = [l.split() for l in open(sys.argv[1]) if not l.startswith("#")]
names = [r[0] for r in rows]
out = []
for i, r in enumerate(rows):
    v = [(float(x), names[j]) for j, x in enumerate(r[1:]) if j != i]
    v.sort(key=lambda x: -x[0])
    keep = v[:int(t)] if mode == "topk" else [x for x in v if x[0] >= t]
    out.append(names[i] + ":" + ",".join(n for _, n in keep))
print(";".join(out))' square.txt "$1" "$2"; }
for p in 1 2; do
    "$D" cmp $O -p $p --symmetric-containment --similarity-threshold 0.5 $F > thr.txt 2>/dev/null
    check "threshold 0.5 -p $p" "$(expect thr 0.5)" "$(lists thr.txt)"
    "$D" cmp $O -p $p --symmetric-containment --topk 2 $F > topk.txt 2>/dev/null
    check "topk 2 -p $p" "$(expect topk 2)" "$(lists topk.txt)"
    for g in 0.5 0.5E; do
        c=$("$D" cmp $O -p $p --symmetric-containment --greedy $g $F 2>/dev/null | grep '^Cluster-' \
            | python3 -c 'import sys; print(sorted(",".join(sorted(x.split(":")[0] for x in l.split()[1:])) for l in sys.stdin))')
        check "greedy $g -p $p" "['c0_0.fa,c0_1.fa,c0_2.fa', 'c1_0.fa,c1_1.fa,c1_2.fa']" "$c"
    done
done
# Mash distance ranks ascending: c0_1 (1% divergent) is closer to c0_0 than c0_2.
"$D" cmp $O --mash-distance --topk 2 $F > mash.txt 2>/dev/null
check "mash-distance topk 2 order for c0_0" "c0_0.fa:c0_1.fa,c0_2.fa" "$(lists mash.txt | cut -d';' -f1)"
rm -f ./*.kmerset* ./*.kmercounts*
exit $fail
