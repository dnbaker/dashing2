#!/bin/bash
# --union-size of two inputs with no k-mers in common must be |A| + |B|, not 0.
# Usage: DASHING2=/path/to/dashing2 test/regress/union-size-disjoint.sh
D=${DASHING2:-$(dirname "$0")/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
python3 - "$T" <<'PY'
import random, sys
r = random.Random(1)
for name in "ab":
    open(f"{sys.argv[1]}/{name}.fa", "w").write(f">{name}\n" + "".join(r.choice("ACGT") for _ in range(2000)) + "\n")
PY
fail=0
for mode in --full; do
  got=$("$D" cmp $mode -k 21 -S 1024 --union-size "$T/a.fa" "$T/b.fa" 2>/dev/null | awk -F'\t' '!/^#/{print $3; exit}')
  # Each 2 kbp random sequence has 1980 distinct 21-mers and they share none, so the union is about 3960.
  if python3 -c "import sys; sys.exit(0 if abs(float('$got') - 3960) < 400 else 1)"; then
    echo "PASS union-size $mode: got $got expected about 3960"
  else
    echo "FAIL union-size $mode: got $got expected about 3960"; fail=1
  fi
done
exit $fail
