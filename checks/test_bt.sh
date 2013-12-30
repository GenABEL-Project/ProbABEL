#!/bin/bash
# This script runs checks on ProbABEL's palogist module for
# binary traits.

echo "Analysing BT..."

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
palogist=${PA_BINDIR}/palogist

# Redirect all output to file descriptor 3 to /dev/null except if
# the first argument is "verbose" then redirect handle 3 to stdout
exec 3>/dev/null
if [ "$1" = "verbose" ]; then
    echo "Verbose mode ON"
    exec 3>&1
fi

$palogist \
    -p ${inputdir}/logist_data.txt \
    -d ${inputdir}/test.mldose \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 \
    -o logist \
    >& 3

$palogist \
    -p ${inputdir}/logist_data.txt \
    -d ${inputdir}/test.dose.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 \
    -o logist_fv \
    >& 3


run_diff logist_add.out.txt \
    logist_fv_add.out.txt \
    "BT check: dose vs. dose_fv"

$palogist \
    -p ${inputdir}/logist_data.txt \
    -d ${inputdir}/test.mlprob \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o logist_prob \
    >& 3

$palogist \
    -p ${inputdir}/logist_data.txt \
    -d ${inputdir}/test.prob.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o logist_prob_fv \
    >& 3

for model in add domin over_domin recess 2df; do
    run_diff logist_prob_${model}.out.txt \
        logist_prob_fv_${model}.out.txt \
        "BT check ($model model): prob vs. prob_fv"
done

run_diff logist_prob_add.out.txt \
    logist_add.out.txt \
    "BT check: prob vs. dose" \
    -I beta_SNP

run_diff logist_prob_fv_add.out.txt \
    logist_fv_add.out.txt \
    "BT check: prob_fv vs. dose_fv" \
    -I beta_SNP
