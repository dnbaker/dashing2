#!/usr/bin/env bash
# Edge cases of the ratio measures. Two inputs without k-mers are identical
# (similarity, containment and symmetric containment 1, Mash distance 0); an
# empty input shares nothing with a non-empty one (0). The Mash distance of
# identical inputs is 0, printed without a minus sign.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
: > e1.fa; : > e2.fa
python3 -c "import random; R = random.Random(3); open('x.fa', 'w').write('>x\n' + ''.join(R.choice('ACGT') for _ in range(5000)) + '\n')"
cp x.fa x2.fa
fail=0
check() { # name got expected (compared as printed)
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"; else echo "FAIL $1: got '$2' expected $3"; fail=1; fi
}
val() { # row file, column file, options...
    local r=$1 c=$2; shift 2
    "$D2" cmp -k 21 --square "$@" "$r" "$c" 2>/dev/null | awk -v f="$r" '$1 == f {print $3}'
}
for mode in --set -J -B; do
    check "empty_empty_similarity_$mode" "$(val e1.fa e2.fa $mode)" 1
    check "empty_empty_mash_$mode" "$(val e1.fa e2.fa $mode --mash-distance)" 0
    check "empty_x_containment_$mode" "$(val e1.fa x.fa $mode --containment)" 0
    check "x_empty_symmetric_containment_$mode" "$(val x.fa e1.fa $mode --symmetric-containment)" 0
    check "identical_mash_$mode" "$(val x.fa x2.fa $mode --mash-distance)" 0
done
check identical_mash_default "$(val x.fa x2.fa --mash-distance)" 0
rm -f ./*.kmerset* ./*.kmer*
exit $fail
