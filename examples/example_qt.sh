#!/bin/sh
# This script runs checks on ProbABEL's palinear module for
# quantitative traits.

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
echo "analysing QT"
if [ -z ${srcdir} ]; then
    srcdir="."
fi

echo "base analysis"
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.mldose \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 \
    -o height_base

../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.dose.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 \
    -o height_base_fv

echo -n "QT check: dose vs. dose_fv"
run_diff height_base_add.out.txt height_base_fv_add.out.txt


echo "Option --allcov"
../src/palinear -p ${srcdir}/height.txt \
    -d ${srcdir}/test.mldose \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --allcov \
    -o height_allcov
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.dose.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --allcov \
    -o height_allcov_fv

echo -n "QT check: allcov: dose vs. dose_fv"
run_diff height_allcov_add.out.txt height_allcov_fv_add.out.txt


echo "Option --interaction=1"
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.mldose \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --interaction=1 \
    -o height_int1
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.dose.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --interaction=1 \
    -o height_int1_fv

echo -n "QT check: interactions: dose vs. dose_fv"
run_diff height_int1_add.out.txt height_int1_fv_add.out.txt


echo "Option --robust"
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.mldose \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --robust \
    -o height_robust
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.dose.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --robust \
    -o height_robust_fv

echo -n "QT check: robust: dose vs. dose_fv"
run_diff height_robust_add.out.txt height_robust_fv_add.out.txt


echo "Option --robust --interaction=1"
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.mldose \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --robust --interaction=1 \
    -o height_robust_int1
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.dose.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --robust --interaction=1 \
    -o height_robust_int1_fv

echo -n "QT check: robust & interaction: dose vs. dose_fv"
run_diff height_robust_int1_add.out.txt height_robust_int1_fv_add.out.txt


echo "Option --ngp=2, mlprob file"
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.mlprob \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --ngp=2 \
    -o height_ngp2
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.prob.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --ngp=2 \
    -o height_ngp2_fv

# echo -n "QT check: dose vs prob (additive model)"
# run_diff height_base_add.out.txt height_ngp2_add.out.txt

for model in add domin over_domin recess 2df; do
    echo -n "QT check ($model model): prob vs. prob_fv"
    run_diff height_ngp2_${model}.out.txt height_ngp2_fv_${model}.out.txt
done


echo "Option --ngp=2 --allcov"
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.mlprob \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --ngp=2 --allcov \
    -o height_ngp2_allcov
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.prob.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --ngp=2 --allcov \
    -o height_ngp2_allcov_fv

for model in add domin over_domin recess 2df; do
    echo -n "QT check --allcov ($model model): prob vs. prob_fv"
    run_diff height_ngp2_allcov_${model}.out.txt height_ngp2_allcov_fv_${model}.out.txt
done


echo "Option --ngp=2 --interaction=1"
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.mlprob \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --ngp=2 --interaction=1 \
    -o height_ngp2_int1
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.prob.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --ngp=2 --interaction=1 \
    -o height_ngp2_int1_fv

for model in add domin over_domin recess 2df; do
    echo -n "QT check --interactions ($model model): prob vs. prob_fv"
    run_diff height_ngp2_int1_${model}.out.txt height_ngp2_int1_fv_${model}.out.txt
done


echo "Option --ngp=2 --robust"
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.mlprob \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --ngp=2 --robust \
    -o height_ngp2_robust
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.prob.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --ngp=2 --robust \
    -o height_ngp2_robust_fv

for model in add domin over_domin recess 2df; do
    echo -n "QT check --robust ($model model): prob vs. prob_fv"
    run_diff height_ngp2_robust_${model}.out.txt height_ngp2_robust_fv_${model}.out.txt
done


echo "Option --ngp=2 --robust --interaction=1"
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.mlprob \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --ngp=2 --robust --interaction=1 \
    -o height_ngp2_robust_int1
../src/palinear \
    -p ${srcdir}/height.txt \
    -d ${srcdir}/test.prob.fvi \
    -i ${srcdir}/test.mlinfo \
    -m ${srcdir}/test.map \
    -c 19 --ngp=2 --robust --interaction=1 \
    -o height_ngp2_robust_int1_fv

for model in add domin over_domin recess 2df; do
    echo -n "QT check --robust --interactions ($model model): prob vs. prob_fv"
    run_diff height_ngp2_robust_int1_${model}.out.txt height_ngp2_robust_int1_fv_${model}.out.txt
done
