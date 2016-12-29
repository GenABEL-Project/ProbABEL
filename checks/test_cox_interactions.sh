#!/bin/bash
# This script runs checks on ProbABEL's pacoxph module

echo "Analysing Cox model with interactions..."

scriptdir=$(dirname $0)

if [ -z ${PA_BINDIR} ]; then
    PA_BINDIR="${scriptdir}/../src/"
fi
if [ -z ${srcdir} ]; then
    srcdir="."
    PA_BINDIR=${scriptdir}/../src/
fi

. ${scriptdir}/run_diff.sh

inputdir=${scriptdir}/inputfiles
pacoxph=${PA_BINDIR}/pacoxph

# Redirect all output to file descriptor 3 to /dev/null except if
# the first argument is "verbose" then redirect handle 3 to stdout
exec 3>/dev/null
if [ "$1" = "verbose" ]; then
    echo "Verbose mode ON"
    exec 3>&1
fi

for covar in {1..3}; do
    echo "Option --interaction=${covar}"
    $pacoxph \
        -p ${inputdir}/coxph_data.txt \
        -d ${inputdir}/test.mldose \
        -i ${inputdir}/test.mlinfo \
        -m ${inputdir}/test.map \
        -c 19 \
        --interaction=${covar} \
        -o coxph_dose_int${covar} \
        >& 3
    $pacoxph \
        -p ${inputdir}/coxph_data.txt \
        -d ${inputdir}/test.mlprob \
        -i ${inputdir}/test.mlinfo \
        -m ${inputdir}/test.map \
        --ngpreds=2 \
        -c 19 \
        --interaction=${covar} \
        -o coxph_prob_int${covar} \
        >& 3

    run_diff coxph_dose_int${covar}_add.out.txt \
             coxph_prob_int${covar}_add.out.txt \
             "pacoxph check --interaction=${covar}: dose vs. prob"
done
