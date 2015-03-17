#!/bin/bash
# This script runs checks on ProbABEL's pacoxph module
# These are almost identical to the checks in test_cox.sh, but this
# time without covariate data (see bug #1266).


echo "Analysing Cox model without covariates..."

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

$pacoxph \
    -p ${inputdir}/coxph_data_nocovar.txt \
    -d ${inputdir}/test.mldose \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 \
    -o coxph_dose_nocovar \
    >& 3

$pacoxph \
    -p ${inputdir}/coxph_data_nocovar.txt \
    -d ${inputdir}/test.dose.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 \
    -o coxph_dose_nocovar_fv \
    >& 3

run_diff coxph_dose_nocovar_add.out.txt coxph_dose_nocovar_fv_add.out.txt \
    "pacoxph check: dose vs. dose_fv"


$pacoxph \
    -p ${inputdir}/coxph_data_nocovar.txt \
    -d ${inputdir}/test.mlprob \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o coxph_prob_nocovar \
    >& 3

run_diff coxph_dose_nocovar_add.out.txt coxph_prob_nocovar_add.out.txt \
    "pacoxph check: dose vs. prob" -I SNP


$pacoxph \
    -p ${inputdir}/coxph_data_nocovar.txt \
    -d ${inputdir}/test.prob.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o coxph_prob_nocovar_fv \
    >& 3

for model in add domin recess over_domin 2df; do
    run_diff coxph_prob_nocovar_${model}.out.txt \
        coxph_prob_nocovar_fv_${model}.out.txt \
        "pacoxph check ($model model): prob vs. prob_fv"
done
