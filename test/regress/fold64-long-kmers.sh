#!/bin/bash
# With -2 (128-bit k-mers), sketches must use all bits of each k-mer. Two 40-mers that share their
# last 32 bases but differ in the first 8 are different k-mers, so their similarity must be 0.
# Usage: DASHING2=/path/to/dashing2 test/regress/fold64-long-kmers.sh
D=${DASHING2:-$(dirname "$0")/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
python3 - "$T" <<'PY'
import random, sys
r = random.Random(4)
tail = "".join(r.choice("ACGT") for _ in range(32))
open(f"{sys.argv[1]}/h1.fa", "w").write(">h1\nACGTACGT" + tail + "\n")
open(f"{sys.argv[1]}/h2.fa", "w").write(">h2\nTTGCATGA" + tail + "\n")
PY
fail=0
for mode in --full ""; do
  got=$(cd "$T" && "$D" cmp -2 --no-canon -k 40 $mode h1.fa h2.fa 2>/dev/null | awk -F'\t' '!/^#/{print $3; exit}')
  if [ "$got" = 0 ]; then echo "PASS ${mode:-oph}: got $got expected 0"; else echo "FAIL ${mode:-oph}: got $got expected 0"; fail=1; fi
done
exit $fail
