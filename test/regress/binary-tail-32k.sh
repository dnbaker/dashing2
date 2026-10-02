#!/usr/bin/env bash
# Rows still queued when the comparisons finish are written in 128 KiB chunks.
# The last chunk must hold the remainder, which is a whole chunk when a queued
# block is an exact multiple of 128 KiB. Fixtures: 256 inputs of 300 bp, compared
# with --square -p 128 --batch-size 128, so each queued block is 128 rows.
# Binary: a row is 256 floats, a block 128 * 256 * 4 = 131072 bytes, and the whole
# matrix 256 * 256 * 4 = 262144 bytes. Text: the inputs are identical and their
# paths 511 characters long, so a row is 511 + 256 * len("\t1") + 1 = 1024 bytes,
# a block 131072 bytes, and all 256 rows must be present.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import os, random
r = random.Random(2)
with open('rand.txt', 'w') as f:
    for i in range(256):
        open('g%03d.fa' % i, 'w').write('>g\n' + ''.join(r.choice('ACGT') for _ in range(300)) + '\n')
        f.write('g%03d.fa\n' % i)
s = ''.join(r.choice('ACGT') for _ in range(300))
d = 'd' * 200 + '/' + 'e' * 200
os.makedirs(d)
with open('same.txt', 'w') as f:
    for i in range(256):
        name = d + '/' + 'x' * (511 - len(d) - 8) + '%04d.fa' % i
        open(name, 'w').write('>s\n' + s + '\n')
        f.write(name + '\n')
PY
fail=0
check() { # name, expected, got
    if [ "$2" = "$3" ]; then echo "PASS $1: got $3 expected $2"; else echo "FAIL $1: got $3 expected $2"; fail=1; fi
}
for p in 64 128; do
    "$D" cmp -k 15 -S 64 --square --binary-output -p $p --batch-size 128 --cmpout out.bin -F rand.txt 2>/dev/null
    check "binary -p $p bytes" 262144 "$(wc -c < out.bin | tr -d ' ')"
    cp out.bin out$p.bin
    "$D" cmp -k 15 -S 64 --square -p $p --batch-size 128 --cmpout out.txt -F same.txt 2>/dev/null
    check "text -p $p rows" 256 "$(grep -vc '^#' out.txt)"
    check "text -p $p row length" 1024 "$(grep -v '^#' out.txt | awk '{print length($0) + 1}' | sort -u | tr '\n' ' ' | sed 's/ $//')"
done
check "binary -p 128 equals -p 64" same "$(cmp -s out64.bin out128.bin && echo same || echo differ)"
rm -f ./*.kmerset* ./*.kmercounts*
exit $fail
