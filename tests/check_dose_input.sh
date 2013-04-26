#!/bin/sh
#
# This script tests whether dose data without a MaCH/minimac-style
# arrow is read correctly by palinear (and by palogist, since reading
# genetic data is done using the same code).

echo "Checking palinear with dose data without '->'"
if [ -z ${srcdir} ]; then
    srcdir="."
fi

exampledir="${srcdir}/../examples/"
results="${srcdir}/verified_results/"
outfile="height_base_add.out.txt"

sed 's/^[[:digit:]]*->//' $exampledir/test.mldose > test.mldose


../src/palinear \
    -p ${exampledir}/height.txt \
    -d test.mldose \
    -i ${exampledir}/test.mlinfo \
    -m ${exampledir}/test.map \
    -c 19 \
    -o height_base > /dev/null

echo -n "  Verifying $outfile: "
if diff $outfile $results/$outfile; then
    echo -e "\t\tOK"
else
    echo -e "\t\tFAILED"
    exit 1
fi
