#!/usr/bin/env bash
# `sketch --cache -o X` must write the same stacked file X whether the sketches
# are computed or loaded from the cache. With --set and --countdict the cached
# path left each input's cardinality at -1 and its bottom-k registers at zero,
# so LSH-assisted --topk and --similarity-threshold found different neighbors
# for cached inputs. Fixtures: two random 6 kbp genomes (the second has 3%
# substitutions), each with a 300 bp block repeated four times so that k-mer
# counts differ. For each mode, the stacked file from a run without --cache must
# be byte-identical to the one from a second --cache run, which loads every
# sketch from the cache written by the first. -B and -P are left out because
# v2.1.20 writes them an all-zero stacked file without --cache, a separate bug.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(11)
rep = ''.join(r.choice('ACGT') for _ in range(300))
a = ''.join(r.choice('ACGT') for _ in range(6000)) + rep * 4
b = list(a)
for p in r.sample(range(len(b)), len(b) * 3 // 100):
    b[p] = r.choice('ACGT'.replace(b[p], ''))
open('a.fa', 'w').write('>a\n' + a + '\n')
open('b.fa', 'w').write('>b\n' + ''.join(b) + '\n')
PY
cards() {  # cardinalities in a stacked file: uint64 count, uint64 sketch size, then one double per input
    python3 -c 'import struct,sys; d=open(sys.argv[1],"rb").read(); n=struct.unpack_from("<Q",d)[0]; print(",".join("%g" % c for c in struct.unpack_from("<%dd" % n, d, 16)))' "$1"
}
fail=0
for mode in "" --full --set --countdict "--set -2" "--countdict -m 2"; do
    name=${mode:-oneperm}; name=${name//[- ]/}
    mkdir -p "$name/fresh" "$name/cache"
    cp a.fa b.fa "$name/fresh"; cp a.fa b.fa "$name/cache"
    (cd "$name/fresh" && "$D" sketch -k 21 -S 256 $mode -o X a.fa b.fa) >/dev/null 2>&1
    (cd "$name/cache" && "$D" sketch -k 21 -S 256 $mode --cache -o X1 a.fa b.fa &&
                         "$D" sketch -k 21 -S 256 $mode --cache -o X a.fa b.fa) >/dev/null 2>&1
    got=$(cards "$name/cache/X"); want=$(cards "$name/fresh/X")
    if cmp -s "$name/fresh/X" "$name/cache/X"; then
        echo "PASS $name: cached stacked file identical to fresh (cardinalities got $got expected $want)"
    else
        echo "FAIL $name: cached stacked file differs from fresh (cardinalities got $got expected $want)"; fail=1
    fi
done
exit $fail
