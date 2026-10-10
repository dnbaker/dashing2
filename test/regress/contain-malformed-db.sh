#!/usr/bin/env bash
# contain must reject malformed k-mer databases with a diagnostic, not crash.
# The 24-byte database header was guarded by an assert only, so release builds
# read (and divided by) whatever a truncated or garbage file held: a 3-byte
# file died with SIGFPE, a zero sketch size in the header died the same way.
# Fixtures: a 3-byte file, a 24-byte header with sketchsize 0, and a valid
# database as a control.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
printf 'ABC' > tiny.kmer64
python3 - <<'PY'
import struct
# inputtype DNA, sketchsize 0, k 0, w 0, seed 0: a header that used to reach
# (size - 24) % sketchsize and size / sketchsize / 8
open('zeross.kmer64', 'wb').write(struct.pack('<IIIIQ', 0, 0, 0, 0, 0))
PY
python3 -c "
import random
r = random.Random(3)
open('x.fa', 'w').write('>x\n' + ''.join(r.choice('ACGT') for _ in range(5000)) + '\n')
"
fail=0
check() { # name, file, expected message fragment
    set +e
    out=$("$D" contain "$2" x.fa 2>&1 >/dev/null)
    rc=$?
    set -e
    if [ $rc -ne 0 ] && [ $rc -ne 139 ] && [ $rc -ne 136 ] && [ $rc -ne 134 ] && echo "$out" | grep -q "$3"; then
        echo "FAIL $1: unexpected exit $rc"
        fail=1
    elif [ $rc -eq 0 ]; then
        echo "FAIL $1: exit 0 on a malformed database"
        fail=1
    elif [ $rc -eq 139 ] || [ $rc -eq 136 ]; then
        echo "FAIL $1: died on signal (exit $rc): $(echo "$out" | tail -1)"
        fail=1
    elif ! echo "$out" | grep -q "$3"; then
        echo "FAIL $1: exit $rc without the expected message"
        fail=1
    else
        echo "PASS $1: rejected with a diagnostic"
    fi
}
set -e
check "3-byte database" tiny.kmer64 "too small to be a k-mer database"
check "zero sketchsize header" zeross.kmer64 "zero sketch size"
"$D" sketch --full -k 15 --save-kmers -o db x.fa 2>/dev/null
set +e
val=$("$D" contain db.kmer64 x.fa 2>/dev/null | grep -c "100%")
set -e
if [ "$val" = "1" ]; then echo "PASS valid database still reports 100% self-coverage"; else echo "FAIL valid database control: got $val"; fail=1; fi
exit $fail
