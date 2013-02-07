#!/bin/sh
# This script runs checks on ProbABEL's palinear module with MMScore

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
echo "Analysis using MMSCORE"
if [ -z ${srcdir} ]; then
    srcdir="."
fi

../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.mldose \
    --sep="," \
    -o mmscore_dose \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat

../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.dose.fvi \
    --sep="," \
    -o mmscore_dose_fv \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat

echo "mmscore check: dose vs. dose_fv"
diff mmscore_add.out.txt mmscore_fv_add.out.txt


../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.mlprob \
    --ngpreds=2 --sep="," \
    -o mmscore_prob \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat

../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.prob.fvi \
    --ngpreds=2 --sep="," \
    -o mmscore_prob_fv \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat

echo -n "mmscore check (additive model): prob vs. prob_fv "
run_diff mmscore_prob_add.out.txt mmscore_prob_fv_add.out.txt
echo -n "mmscore check (dominant model): prob vs. prob_fv "
run_diff mmscore_prob_domin.out.txt mmscore_prob_fv_domin.out.txt
echo -n "mmscore check (recessive model): prob vs. prob_fv "
run_diff mmscore_prob_recess.out.txt mmscore_prob_fv_recess.out.txt
echo -n "mmscore check (over-dominant model): prob vs. prob_fv "
run_diff mmscore_prob_over_domin.out.txt mmscore_prob_fv_over_domin.out.txt
echo -n "mmscore check (2df model): prob vs. prob_fv "
run_diff mmscore_prob_2df.out.txt mmscore_prob_fv_2df.out.txt

# Commented out because of slightly different output formats. We need
# something smart here.
# echo "mmscore check: prob vs. dose"
# diff mmscore_prob_add.out.txt mmscore_dose_add.out.txt

# echo "mmscore check: prob_fv vs. dose_fv"
# diff mmscore_prob_fv_add.out.txt mmscore_dose_fv_add.out.txt
