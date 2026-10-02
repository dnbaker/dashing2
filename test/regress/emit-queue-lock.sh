#!/bin/bash
# Every row of a cmp matrix must be printed. The writer thread that formats rows
# shares a queue with the compute threads; before that queue was read under its
# lock, a row could rarely be dropped with exit status 0 (seen on Apple Silicon).
# The failure is timing-dependent, so this repeats a multi-threaded comparison
# and checks that each run prints exactly one row per input, in order. A clean
# run does not prove the absence of the race; building with -fsanitize=thread
# shows it directly.
# Usage: DASHING2=/path/to/dashing2 [REPS=30] test/regress/emit-queue-lock.sh
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
REPS=${REPS:-30}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(13)
for i in range(40):
    open(f"f{i:02d}.fa", "w").write(f">f{i}\n" + "".join(r.choice("ACGT") for _ in range(2000)) + "\n")
PY
expected=$(printf 'f%02d.fa\n' $(seq 0 39))
fail=0
for mode in "--set" "--full" "--full --square"; do
    bad=0
    for i in $(seq 1 "$REPS"); do
        rm -f ./*.kmerset* 2>/dev/null
        got=$("$D" cmp -p 4 -k 21 $mode f*.fa 2>/dev/null | awk -F'\t' '!/^#/{sub(/ +$/, "", $1); print $1}')
        [ "$got" = "$expected" ] || bad=$((bad + 1))
    done
    if [ $bad = 0 ]; then echo "PASS cmp $mode: all rows present in $REPS runs"; else echo "FAIL cmp $mode: rows missing or out of order in $bad of $REPS runs"; fail=1; fi
done
exit $fail
