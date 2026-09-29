#!/usr/bin/env bash
# --cache must not reuse a sketch when its input file changed after the sketch
# was cached. Here b.fa is overwritten with a copy of a.fa after caching, so
# the comparison must report 1 (a file against an identical copy).
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
for name in 'ab':
    open(name + '.fa', 'w').write('>' + name + '\n' + ''.join(r.choice('ACGT') for _ in range(5000)) + '\n')
PY
cp b.fa b.orig
val() { "$D" cmp "$@" a.fa b.fa 2>/dev/null | awk 'NR==4{print $3}'; }
fail=0
for mode in --full --set; do
    before=$(val $mode -k 21 --cache)
    # Replace b.fa with a copy of a.fa and date it later than the cache files.
    python3 -c 'import os, shutil, time; shutil.copy("a.fa", "b.fa"); t = time.time() + 5; os.utime("b.fa", (t, t))'
    after=$(val $mode -k 21 --cache)
    rm -f a.fa.* b.fa.*; cp b.orig b.fa
    if [ "$after" = 1 ]; then
        echo "PASS $mode after input changed: got $after expected 1 (before change: $before)"
    else
        echo "FAIL $mode after input changed: got $after expected 1 (before change: $before)"; fail=1
    fi
done
exit $fail
