#!/usr/bin/env bash
# A -2 (128-bit k-mer) run with --cache must not load sketches cached by a
# 64-bit run on the same inputs. Invariant: the cached -2 result equals a fresh
# -2 result computed without any cache.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
a = ''.join(r.choice('ACGT') for _ in range(5000))
b = list(a)
for i in range(0, len(b), 97):
    b[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[b[i]]
open('a.fa', 'w').write('>a\n' + a + '\n')
open('b.fa', 'w').write('>b\n' + ''.join(b) + '\n')
PY
val() { "$D" cmp "$@" a.fa b.fa 2>/dev/null | awk 'NR==4{print $3}'; }
fail=0
for mode in --full --bagminhash; do
    fresh=$(val -2 $mode -k 40)
    val $mode -k 40 --cache >/dev/null      # populate the cache with 64-bit sketches
    cached=$(val -2 $mode -k 40 --cache)
    rm -f a.fa.* b.fa.*
    if [ -n "$fresh" ] && [ "$cached" = "$fresh" ]; then
        echo "PASS -2 $mode after 64-bit cache: got $cached expected $fresh"
    else
        echo "FAIL -2 $mode after 64-bit cache: got $cached expected $fresh"; fail=1
    fi
done
exit $fail
