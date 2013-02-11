#!/bin/sh
# This script runs checks on ProbABEL's palinear module for
# quantitative traits.

run_diff()
{
    # This function is run after each check. It needs two arguments:
    # $1: first file to compare
    # $2: second file to compare
    if diff "$1" "$2"; then
        echo -e "\t\tOK"
    else
        echo -e "\t\tFAILED"
        exit 1
    fi
}


# ---- The checks start here ----
echo "analysing BT"
if [ -z ${srcdir} ]; then
    srcdir="."
fi

# Redirect all output to file descriptor 3 to /dev/null except if
# the first argument is "verbose" then redirect handle 3 to stdout
exec 3>/dev/null
if [ "$1" = "verbose" ]; then
    echo "Verbose mode ON"
    exec 3>&1
fi

../src/palogist \
    -p ${srcdir}/logist_data.txt \
    -d ${srcdir}/test.mldose \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 \
    -o logist \
    >& 3

../src/palogist \
    -p ${srcdir}/logist_data.txt \
    -d ${srcdir}/test.dose.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 \
    -o logist_fv \
    >& 3

echo "BT check: dose vs. dose_fv"
run_diff logist_add.out.txt logist_fv_add.out.txt

../src/palogist \
    -p ${srcdir}/logist_data.txt \
    -d ${srcdir}/test.mlprob \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o logist_prob \
    >& 3

../src/palogist \
    -p ${srcdir}/logist_data.txt \
    -d ${srcdir}/test.prob.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o logist_prob_fv \
    >& 3

for model in add domin over_domin recess 2df; do
    echo -n "BT check ($model model): prob vs. prob_fv"
    run_diff logist_prob_${model}.out.txt logist_prob_fv_${model}.out.txt
done

# Commented out because of slightly different output formats. We need
# something smart here.
# echo "BT check: prob vs. dose"
# diff logist_prob_add.out.txt logist_add.out.txt
# echo "BT check: prob_fv vs. dose_fv"
# diff logist_prob_fv_add.out.txt logist_fv_add.out.txt
