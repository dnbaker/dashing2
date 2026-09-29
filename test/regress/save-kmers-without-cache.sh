#!/usr/bin/env bash
# -s/--save-kmers and -N/--save-kmercounts with -B (BagMinHash) or --prob
# (ProbMinHash) must work without --cache and write the same k-mer and count
# files as the --cache run, which opens the output stream beforehand.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random
r = random.Random(2)
s = ''.join(r.choice('ACGT') for _ in range(5000))
open('a.fa', 'w').write(f'>a\n{s}\n')
open('b.fa', 'w').write(f'>b0\n{s}\n>b1\n{s[:2000]}\n')
PY
clean() { rm -f ./*.kmer* ./*.bmh ./*.pmh o.*; }
for mode in -B --prob; do
    for save in -s -N; do
        name="$mode $save"
        clean
        "$D2" sketch -k 21 $mode $save --cache -o o.ref a.fa b.fa >/dev/null 2>&1
        ref=$(cat o.ref.kmer64 o.ref.kmercounts.f64 2>/dev/null | cksum)
        clean
        "$D2" sketch -k 21 $mode $save -o o.new a.fa b.fa >/dev/null 2>&1; rc=$?
        new=$(cat o.new.kmer64 o.new.kmercounts.f64 2>/dev/null | cksum)
        if [ $rc -eq 0 ]; then echo "PASS $name exit status: got $rc expected 0"
        else echo "FAIL $name exit status: got $rc expected 0"; fail=1; fi
        if [ "$new" = "$ref" ]; then echo "PASS $name saved files match --cache run"
        else echo "FAIL $name saved files: got checksum $new expected $ref (--cache run)"; fail=1; fi
    done
done
clean
exit $fail
