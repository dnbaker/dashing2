#!/bin/bash
# With k > 32 (rolling hash path) and several threads, a file and an identical copy must have similarity 1.
# Usage: DASHING2=/path/to/dashing2 test/regress/rolling-hasher-threads.sh
D=${DASHING2:-$(dirname "$0")/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
python3 - "$T" <<'PY'
import random, sys
r = random.Random(3)
x, y = ("".join(r.choice("ACGT") for _ in range(2000)) for _ in range(2))
for name, s in [("x", x), ("x_copy", x), ("y", y), ("z", y[:1000] + x[1000:])]:
    open(f"{sys.argv[1]}/{name}.fa", "w").write(f">{name}\n{s}\n")
PY
fail=0
for trial in 1 2 3 4 5; do
  got=$(cd "$T" && rm -f *.kmerset* && "$D" cmp -p 4 -k 40 --set x.fa x_copy.fa y.fa z.fa 2>/dev/null | awk -F'\t' '!/^#/{print $3; exit}')
  if [ "$got" = 1 ]; then echo "PASS trial $trial: J(x, copy of x) = $got"; else echo "FAIL trial $trial: J(x, copy of x) = $got, expected 1"; fail=1; fi
done
exit $fail
