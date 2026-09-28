#!/bin/bash
# --protein must parse amino acids: exact k-mer Jaccard must match a direct count.
# Usage: DASHING2=/path/to/dashing2 test/regress/protein-alphabet.sh
D=${DASHING2:-$(dirname "$0")/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
read -r jm_exp jq_exp < <(python3 - "$T" <<'PY'
import random, sys
r = random.Random(5)
aa = "ACDEFGHIKLMNPQRSTVWY"
p, q = ("".join(r.choice(aa) for _ in range(2000)) for _ in range(2))
m = "".join(r.choice(aa.replace(c, "")) if r.random() < 0.03 else c for c in p)  # 3% substitutions
for name, s in [("p", p), ("q", q), ("m", m)]:
    open(f"{sys.argv[1]}/{name}.fa", "w").write(f">{name}\n{s}\n")
ks = lambda s: {s[i:i+5] for i in range(len(s) - 4)}
J = lambda a, b: len(ks(a) & ks(b)) / len(ks(a) | ks(b))
print(J(p, m), J(p, q))
PY
)
v() { (cd "$T" && rm -f *.kmerset* && "$D" cmp --protein -k 5 --set "$@" 2>/dev/null | awk -F'\t' '!/^#/{print $3; exit}'); }
fail=0
check() {
  if python3 -c "import sys; sys.exit(0 if abs(float('$2') - $3) < 1e-4 else 1)" 2>/dev/null; then echo "PASS $1: got $2 expected $3"; else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
check "J(p, p with 3% substitutions)" "$(v p.fa m.fa)" "$jm_exp"
check "J(p, unrelated q)" "$(v p.fa q.fa)" "$jq_exp"
exit $fail
