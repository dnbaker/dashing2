#!/bin/bash
# --set --union-size must print |A u B| of the exact canonical k-mer sets, not |A n B|.
# Usage: DASHING2=/path/to/dashing2 test/regress/union-size-set.sh
D=${DASHING2:-$(dirname "$0")/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
exp=$(python3 - "$T" <<'PY'
import random, sys
r = random.Random(2)
a = "".join(r.choice("ACGT") for _ in range(2000))
b = list(a)
for i in range(0, 2000, 50):  # one substitution every 50 bp
    b[i] = "ACGT"[("ACGT".index(b[i]) + 1) % 4]
b = "".join(b)
open(f"{sys.argv[1]}/a.fa", "w").write(">a\n" + a + "\n")
open(f"{sys.argv[1]}/b.fa", "w").write(">b\n" + b + "\n")
rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
kset = lambda s: {min(s[i:i+21], rc(s[i:i+21])) for i in range(len(s) - 20)}
print(len(kset(a) | kset(b)))
PY
)
got=$(cd "$T" && "$D" cmp -k 21 --set --union-size a.fa b.fa 2>/dev/null | awk -F'\t' '!/^#/{print $3; exit}')
rm -f "$T"/*.kmerset*
if python3 -c "import sys; sys.exit(0 if abs(float('$got') - $exp) < 0.5 else 1)"; then
  echo "PASS set union-size: got $got expected $exp"; exit 0
else
  echo "FAIL set union-size: got $got expected $exp"; exit 1
fi
