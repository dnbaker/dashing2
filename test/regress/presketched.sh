#!/usr/bin/env bash
# cmp --presketched on files written by dashing2 itself. Each loaded result must
# equal the same comparison run directly on the FASTA inputs, and files that
# cannot be compared must be rejected with a message rather than a crash.
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
fail=0
report() {  # report <name> <got> <expected>
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got ${2:-nothing} expected $3"; fail=1; fi
}
cell() { awk 'NR==4{print $3}'; }
# rejects <message pattern> <args>: "rejected" if dashing2 exits with an error status (not a signal) and the message
rejects() {
    local pat=$1; shift
    local out; out=$("$D" "$@" 2>&1); local rc=$?
    if [ $rc -gt 0 ] && [ $rc -lt 128 ] && echo "$out" | grep -q "$pat"; then echo rejected; else echo "exit $rc"; fi
}
# 1. Stacked sketches written to a path without an extension (the example in cmp's help).
direct=$("$D" cmp --full -k 21 a.fa b.fa 2>/dev/null | cell)
"$D" sketch --full -k 21 -o stacked a.fa b.fa >/dev/null 2>&1
report "stacked file without extension" "$("$D" cmp --full --presketched stacked 2>/dev/null | cell)" "$direct"
# 2. Per-input sketch files keep their paths as row names.
"$D" sketch --full -k 21 --cache a.fa b.fa >/dev/null 2>&1
out=$("$D" cmp --presketched a.fa.*.ss b.fa.*.ss 2>/dev/null)
report "per-input .ss files" "$(echo "$out" | cell)" "$direct"
report "per-input .ss row name" "$(echo "$out" | awk 'NR==4{print $1}')" "$(ls a.fa.*.ss)"
# 3. Per-input exact k-mer set files.
direct=$("$D" cmp --set -k 21 a.fa b.fa 2>/dev/null | cell)
"$D" cmp --set -k 21 --cache a.fa b.fa >/dev/null 2>&1
report "per-input .kmerset64 files" "$("$D" cmp --presketched a.fa.*.kmerset64 b.fa.*.kmerset64 2>/dev/null | cell)" "$direct"
# 4. Stacked --set output holds only bottom-k hashes: rejected, not a segfault.
"$D" sketch --set -k 21 -o stacked.kmerset64 a.fa b.fa >/dev/null 2>&1
report "stacked --set file" "$(rejects 'cannot be compared' cmp --presketched stacked.kmerset64)" rejected
# 5. Sketches of different sizes: rejected, not an out_of_range abort.
rm -f a.fa.*.ss
"$D" sketch --full -k 21 -S 512 --cache a.fa >/dev/null 2>&1
report "sketch size mismatch" "$(rejects 'different sizes' cmp --presketched a.fa.*.ss b.fa.*.ss)" rejected
# 6. Files of different sketch types: rejected rather than read as the first type.
report "mixed sketch types" "$(rejects 'different sketch types' cmp --presketched b.fa.*.ss b.fa.*.kmerset64)" rejected
exit $fail
