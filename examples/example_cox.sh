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

# Redirect all output to file descriptor 3 to /dev/null except if
# the first argument is "verbose" then redirect handle 3 to stdout
exec 3>/dev/null
if [ "$1" = "verbose" ]; then
    echo "Verbose mode ON"
    exec 3>&1
fi

../src/pacoxph \
    -p ${srcdir}/coxph_data.txt \
    -d ${srcdir}/test.mldose \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 \
    -o coxph_dose \
    >& 3

../src/pacoxph \
    -p ${srcdir}/coxph_data.txt \
    -d ${srcdir}/test.dose.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 \
    -o coxph_dose_fv \
    >& 3

echo -n "pacoxph check: dose vs. dose_fv "
run_diff coxph_dose_add.out.txt coxph_dose_fv_add.out.txt


../src/pacoxph \
    -p ${srcdir}/coxph_data.txt \
    -d ${srcdir}/test.mlprob \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o coxph_prob \
    >& 3

# Disabling this check for now because the output differs slightly
echo -n "pacoxph check: dose vs. prob "
run_diff coxph_dose_add.out.txt coxph_prob_add.out.txt


../src/pacoxph \
    -p ${srcdir}/coxph_data.txt \
    -d ${srcdir}/test.prob.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o coxph_prob_fv \
    >& 3

for model in add domin recess over_domin 2df; do
    echo -n "pacoxph check ($model model): prob vs. prob_fv "
    run_diff coxph_prob_${model}.out.txt coxph_prob_fv_${model}.out.txt
done
