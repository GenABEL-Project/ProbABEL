#!/bin/bash
# This script runs checks on ProbABEL's palinear module for
# quantitative traits combined with the mmscore option.

echo "Analysis using MMScore..."

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
palinear=${PA_BINDIR}/palinear

# Redirect all output to file descriptor 3 to /dev/null except if
# the first argument is "verbose" then redirect handle 3 to stdout
exec 3>/dev/null
if [ "$1" = "verbose" ]; then
    echo "Verbose mode ON"
    exec 3>&1
fi

$palinear \
    -p ${inputdir}/mmscore_pheno.PHE \
    -i ${inputdir}/mmscore_gen.mlinfo \
    -d ${inputdir}/mmscore_gen.mldose \
    --sep="," \
    -o mmscore_dose \
    --mmscore ${inputdir}/mmscore_InvSigma_aj.sex.age.dat \
    >& 3

$palinear \
    -p ${inputdir}/mmscore_pheno.PHE \
    -i ${inputdir}/mmscore_gen.mlinfo \
    -d ${inputdir}/mmscore_gen.dose.fvi \
    --sep="," \
    -o mmscore_dose_fv \
    --mmscore ${inputdir}/mmscore_InvSigma_aj.sex.age.dat \
    >& 3


run_diff mmscore_dose_add.out.txt \
    mmscore_dose_fv_add.out.txt \
    "mmscore check: dose vs. dose_fv"


$palinear \
    -p ${inputdir}/mmscore_pheno.PHE \
    -i ${inputdir}/mmscore_gen.mlinfo \
    -d ${inputdir}/mmscore_gen.mlprob \
    --ngpreds=2 --sep="," \
    -o mmscore_prob \
    --mmscore ${inputdir}/mmscore_InvSigma_aj.sex.age.dat \
    >& 3

$palinear \
    -p ${inputdir}/mmscore_pheno.PHE \
    -i ${inputdir}/mmscore_gen.mlinfo \
    -d ${inputdir}/mmscore_gen.prob.fvi \
    --ngpreds=2 --sep="," \
    -o mmscore_prob_fv \
    --mmscore ${inputdir}/mmscore_InvSigma_aj.sex.age.dat \
    >& 3

for model in add domin over_domin recess 2df; do
    run_diff mmscore_prob_${model}.out.txt \
        mmscore_prob_fv_${model}.out.txt \
        "mmscore check ($model model): prob vs. prob_fv"
done

run_diff mmscore_prob_add.out.txt \
    mmscore_dose_add.out.txt \
    "mmscore check: prob vs. dose" \
    -I SNP

run_diff mmscore_prob_fv_add.out.txt \
    mmscore_dose_fv_add.out.txt \
    "mmscore check: prob_fv vs. dose_fv" \
    -I SNP
