#!/bin/bash
# This script runs checks on ProbABEL's palinear module for
# quantitative traits.

echo "analysing QT"
if [ -z ${srcdir} ]; then
    srcdir="."
fi

. ${srcdir}/run_diff.sh

inputdir=${srcdir}/inputfiles

# Redirect all output to file descriptor 3 to /dev/null except if
# the first argument is "verbose" then redirect handle 3 to stdout
exec 3>/dev/null
if [ "$1" = "verbose" ]; then
    exec 3>&1
fi

echo "base analysis"
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.mldose \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 \
    -o height_base \
    >& 3
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.dose.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 \
    -o height_base_fv \
    >& 3

#echo -n "QT check: dose vs. dose_fv"
run_diff height_base_add.out.txt \
    height_base_fv_add.out.txt \
    "QT check: dose vs. dose_fv"


echo "Option --allcov"
../src/palinear -p ${inputdir}/height.txt \
    -d ${inputdir}/test.mldose \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --allcov \
    -o height_allcov \
    >& 3
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.dose.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --allcov \
    -o height_allcov_fv \
    >& 3

run_diff height_allcov_add.out.txt \
    height_allcov_fv_add.out.txt \
    "QT check: allcov: dose vs. dose_fv"


echo "Option --interaction=1"
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.mldose \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --interaction=1 \
    -o height_int1 \
    >& 3
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.dose.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --interaction=1 \
    -o height_int1_fv \
    >& 3

run_diff height_int1_add.out.txt \
    height_int1_fv_add.out.txt \
    "QT check: interactions: dose vs. dose_fv"


echo "Option --robust"
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.mldose \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --robust \
    -o height_robust \
    >& 3
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.dose.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --robust \
    -o height_robust_fv \
    >& 3

run_diff height_robust_add.out.txt \
    height_robust_fv_add.out.txt \
    "QT check: robust: dose vs. dose_fv"


echo "Option --robust --interaction=1"
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.mldose \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --robust --interaction=1 \
    -o height_robust_int1 \
    >& 3
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.dose.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --robust --interaction=1 \
    -o height_robust_int1_fv \
    >& 3

run_diff height_robust_int1_add.out.txt \
    height_robust_int1_fv_add.out.txt \
    "QT check: robust & interaction: dose vs. dose_fv"


echo "Option --ngp=2, mlprob file"
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.mlprob \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --ngp=2 \
    -o height_ngp2 \
    >& 3
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.prob.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --ngp=2 \
    -o height_ngp2_fv \
    >& 3


# Remove header from the outputs, because they differ
run_diff height_base_add.out.txt \
    height_ngp2_add.out.txt \
    "QT check: dose vs. prob (additive model)" -I SNP

for model in add domin over_domin recess 2df; do
    run_diff height_ngp2_${model}.out.txt \
        height_ngp2_fv_${model}.out.txt \
        "QT check ($model model): prob vs. prob_fv"
done


echo "Option --ngp=2 --allcov"
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.mlprob \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --ngp=2 --allcov \
    -o height_ngp2_allcov \
    >& 3
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.prob.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --ngp=2 --allcov \
    -o height_ngp2_allcov_fv \
    >& 3

for model in add domin over_domin recess 2df; do
    run_diff height_ngp2_allcov_${model}.out.txt \
        height_ngp2_allcov_fv_${model}.out.txt \
        "QT check --allcov ($model model): prob vs. prob_fv"
done


echo "Option --ngp=2 --interaction=1"
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.mlprob \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --ngp=2 --interaction=1 \
    -o height_ngp2_int1 \
    >& 3
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.prob.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --ngp=2 --interaction=1 \
    -o height_ngp2_int1_fv \
    >& 3

for model in add domin over_domin recess 2df; do
    run_diff height_ngp2_int1_${model}.out.txt \
        height_ngp2_int1_fv_${model}.out.txt \
        "QT check --interactions ($model model): prob vs. prob_fv"
done


echo "Option --ngp=2 --robust"
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.mlprob \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --ngp=2 --robust \
    -o height_ngp2_robust \
    >& 3
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.prob.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --ngp=2 --robust \
    -o height_ngp2_robust_fv \
    >& 3

for model in add domin over_domin recess 2df; do
    run_diff height_ngp2_robust_${model}.out.txt \
        height_ngp2_robust_fv_${model}.out.txt \
        "QT check --robust ($model model): prob vs. prob_fv"
done


echo "Option --ngp=2 --robust --interaction=1"
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.mlprob \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --ngp=2 --robust --interaction=1 \
    -o height_ngp2_robust_int1 \
    >& 3
../src/palinear \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.prob.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --ngp=2 --robust --interaction=1 \
    -o height_ngp2_robust_int1_fv \
    >& 3

for model in add domin over_domin recess 2df; do
    run_diff height_ngp2_robust_int1_${model}.out.txt \
        height_ngp2_robust_int1_fv_${model}.out.txt \
        "QT check --robust --interactions ($model model): prob vs. prob_fv"
done
