#!/bin/sh
# This script runs checks on ProbABEL's palinear module for
# quantitative traits combined with the mmscore option.

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
    -o mmscore \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat \
    >& 3

../src/palinear \
    -p ${srcdir}/mmscore_pheno.PHE \
    -i ${srcdir}/mmscore_gen.mlinfo \
    -d ${srcdir}/mmscore_gen.dose.fvi \
    --sep="," \
    -o mmscore_fv \
    --mmscore ${srcdir}/mmscore_InvSigma_aj.sex.age.dat \
    >& 3

echo "mmscore check: dose vs. dose_fv"
diff mmscore_add.out.txt mmscore_fv_add.out.txt


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
    echo -n "mmscore check ($model model): prob vs. prob_fv"
    run_diff mmscore_prob_${model}.out.txt mmscore_prob_fv_${model}.out.txt
done

# Commented out because of slightly different output formats. We need
# something smart here.
# echo "mmscore check: prob vs. dose"
# diff mmscore_prob_add.out.txt mmscore_add.out.txt

# echo "mmscore check: prob_fv vs. dose_fv"
# diff mmscore_prob_fv_add.out.txt mmscore_fv_add.out.txt
