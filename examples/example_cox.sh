#!/bin/sh
# This script runs checks on ProbABEL's pacoxph module

run_diff()
{
    # This function is run after each check. It needs two arguments:
    # $1: first file to compare
    # $2: second file to compare
    if diff "$1" "$2"; then
        echo -e "\t\tOK"
    else
        echo -e "\t\tFAILED"
#        exit 1
    fi
}


# ---- The checks start here ----
echo "Analysing Cox model..."
if [ -z ${srcdir} ]; then
    srcdir="."
fi

../src/pacoxph \
    -p ${srcdir}/coxph_data.txt \
    -d ${srcdir}/test.mldose \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 \
    -o coxph

../src/pacoxph \
    -p ${srcdir}/coxph_data.txt \
    -d ${srcdir}/test.dose.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 \
    -o coxph_fv

echo -n "pacoxph check: dose vs. dose_fv"
run_diff coxph_add.out.txt coxph_fv_add.out.txt


../src/pacoxph \
    -p ${srcdir}/coxph_data.txt \
    -d ${srcdir}/test.mlprob \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o coxph_prob

echo -n "pacoxph check: dose vs. prob"
run_diff coxph_add.out.txt coxph_prob_add.out.txt


../src/pacoxph \
    -p ${srcdir}/coxph_data.txt \
    -d ${srcdir}/test.prob.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o coxph_fv_prob

echo -n "pacoxph check: prob vs. prob_fv"
run_diff coxph_prob_add.out.txt coxph_fv_prob_add.out.txt
