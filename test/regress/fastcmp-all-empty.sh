#!/usr/bin/env bash
# --fastcmp N log-compresses sketch registers with a and b fitted to the
# smallest and largest register over all inputs. When no input has a k-mer
# there is nothing to fit, and the comparison must still print what the
# uncompressed sketches print (two inputs without k-mers are identical:
# similarity 1, and with --full union and intersection 0) instead of NaN or inf.
# Fixtures: e1.fa and e2.fa hold sequences shorter than k.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
printf '>e1\nACGTACGTAC\n' > e1.fa
printf '>e2\nTTGCA\n' > e2.fa
fail=0
check() { # name got expected
    if [ -n "$2" ] && [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got '$2' expected '$3'"; fail=1; fi
}
# Value printed for the pair (e1, e2).
val() { "$D2" cmp -k 21 --square "$@" e1.fa e2.fa 2>/dev/null | awk '$1 == "e1.fa" {print $3}'; }
check "full_union" "$(val --full --fastcmp 1 --union-size)" 0
check "full_similarity" "$(val --full --fastcmp 1)" 1
for mode in "" --full -B --prob; do
    for meas in "" --union-size --intersection --mash-distance; do
        [ "$mode" = --prob ] && [ -n "$meas" ] && [ "$meas" != --mash-distance ] && continue
        expected=$(val $mode $meas)
        for fc in 1 2 4; do
            check "${mode:-oph}${meas:- --similarity} --fastcmp $fc" "$(val $mode $meas --fastcmp $fc)" "$expected"
        done
    done
done
rm -f ./*.kmerset* ./*.ss ./*.opss
exit $fail
