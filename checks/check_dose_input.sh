#!/bin/bash
#
# This script tests whether dose data without a MaCH/minimac-style
# arrow is read correctly by palinear (and by palogist, since reading
# genetic data is done using the same code).

echo "Checking palinear with dose data without '->'"
# Exit with error when one of the steps in the script fails
set -e

if [ -z ${srcdir} ]; then
    srcdir="."
fi
inputdir="${srcdir}/inputfiles/"
results="${srcdir}/verified_results/"
outfile="height_base_add.out.txt"

sed 's/^[[:digit:]]*->//' $inputdir/test.mldose > test.mldose


../src/palinear \
    -p ${inputdir}/height.txt \
    -d test.mldose \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 \
    -o height_base > /dev/null

blanks="                                                          "
echo -n "  Verifying "
if diff $outfile $results/$outfile; then
    echo -e "${outfile}${blanks:${#outfile}} OK"
else
    echo -e "${outfile}${blanks:${#outfile}} FAILED"
    exit 1
fi
