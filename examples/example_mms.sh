#!/bin/sh
# This script runs checks on ProbABEL's palinear module for
# quantitative traits combined with the mmscore option.

. ./run_diff.sh

echo "analysis using MMScore"
if [ -z ${srcdir} ]; then
    srcdir="."
fi

# Redirect all output to file descriptor 3 to /dev/null except if
# the first argument is "verbose" then redirect handle 3 to stdout
exec 3>/dev/null
if [ "$1" = "verbose" ]; then
    exec 3>&1
fi

../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.mldose \
    --sep="," \
    -o mmscore_dose \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat \
    >& 3

../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.dose.fvi \
    --sep="," \
    -o mmscore_dose_fv \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat \
    >& 3


run_diff mmscore_dose_add.out.txt \
    mmscore_dose_fv_add.out.txt \
    "mmscore check: dose vs. dose_fv"


../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.mlprob \
    --ngpreds=2 --sep="," \
    -o mmscore_prob \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat \
    >& 3

../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.prob.fvi \
    --ngpreds=2 --sep="," \
    -o mmscore_prob_fv \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat \
    >& 3

for model in add domin over_domin recess 2df; do
    run_diff mmscore_prob_${model}.out.txt \
        mmscore_prob_fv_${model}.out.txt \
        "mmscore check ($model model): prob vs. prob_fv"
done

# The following checks are disabled because of the missing LogLik
# column in the prob data
# run_diff mmscore_prob_add.out.txt \
#     mmscore_add.out.txt \
#     "mmscore check: prob vs. dose" \
#     -I beta_SNP

# run_diff mmscore_prob_fv_add.out.txt \
#     mmscore_fv_add.out.txt \
#     "mmscore check: prob_fv vs. dose_fv" \
#     -I beta_SNP
