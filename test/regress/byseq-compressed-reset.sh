#!/usr/bin/env bash
# --parse-by-seq sketches each record separately. With compressed --full
# sketches (--fastcmp-bytes) the sketcher must be cleared between records, or
# record b is sketched as a u b and the similarity drifts toward |a|/|a u b|.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
# a: 3 kbp random; b: first half of a plus 1.5 kbp new sequence.
python3 - <<'PY'
import random
r = random.Random(1)
a = ''.join(r.choice('ACGT') for _ in range(3000))
b = a[:1500] + ''.join(r.choice('ACGT') for _ in range(1500))
open('ab.fa', 'w').write(f'>a\n{a}\n>b\n{b}\n')
PY
# Exact Jaccard of the canonical 21-mer sets of the two records.
exp=$(python3 - <<'PY'
k = 21; comp = str.maketrans('ACGT', 'TGCA')
recs = [l.strip() for l in open('ab.fa') if not l.startswith('>')]
a, b = ({min(s[i:i + k], s[i:i + k].translate(comp)[::-1]) for i in range(len(s) - k + 1)} for s in recs)
print(f'{len(a & b) / len(a | b):.4f}')
PY
)
for opts in "--full" "--full --fastcmp-bytes" "--full --fastcmp-shorts"; do
    got=$("$D2" cmp -k 21 --parse-by-seq $opts ab.fa 2>/dev/null | grep -v '^#' | awk 'NR == 1 {print $3}')
    rm -f ./*.kmer*
    if python3 -c "import sys; sys.exit(abs($got - $exp) > 0.05)"; then echo "PASS by-seq $opts: got $got expected $exp +- 0.05"
    else echo "FAIL by-seq $opts: got $got expected $exp +- 0.05"; fail=1; fi
done
exit $fail
