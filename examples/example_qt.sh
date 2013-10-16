#!/bin/bash
# This script shows some examples on how to run ProbABEL's palinear
# module for quantitative traits.


# Set the name of the directory with the genotype data files
gtdatadir="gtdata"

# Set the path to the executables (palinear/palogist/pacoxph). If you
# installed ProbABEL according to the instructions in the doc/INSTALL
# file, or installed it via your Linux distribution's package manager,
# you can leave this variable empty.
padir="${padir}"

# Using text-based dosage genotype files as input
echo "basic analysis"
${padir}palinear \
    -p height.txt \
    -d ${gtdatadir}/test.mldose \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 \
    -o height_base


# Using filevector (DatABEL) files as dosage genotype input
${padir}palinear \
    -p height.txt \
    -d ${gtdatadir}/test.dose.fvi \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 \
    -o height_base_fv


echo "Option --allcov"
${padir}palinear \
    -p height.txt \
    -d ${gtdatadir}/test.dose.fvi \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 --allcov \
    -o height_allcov_fv


echo "Option --interaction=1"
${padir}palinear \
    -p height.txt \
    -d ${gtdatadir}/test.dose.fvi \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 --interaction=1 \
    -o height_int1_fv


echo "Option --robust"
${padir}palinear \
    -p height.txt \
    -d ${gtdatadir}/test.dose.fvi \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 --robust \
    -o height_robust_fv


echo "Option --robust --interaction=1"
${padir}palinear \
    -p height.txt \
    -d ${gtdatadir}/test.dose.fvi \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 --robust --interaction=1 \
    -o height_robust_int1_fv


# Using text-based probability files as genotype input
echo "Option --ngp=2, mlprob file"
${padir}palinear \
    -p height.txt \
    -d ${gtdatadir}/test.mlprob \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 --ngp=2 \
    -o height_ngp2


# Using filevector (DatABEL) probability files as genotype input
${padir}palinear \
    -p height.txt \
    -d ${gtdatadir}/test.prob.fvi \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 --ngp=2 \
    -o height_ngp2_fv


echo "Option --ngp=2 --allcov"
${padir}palinear \
    -p height.txt \
    -d ${gtdatadir}/test.prob.fvi \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 --ngp=2 --allcov \
    -o height_ngp2_allcov_fv


echo "Option --ngp=2 --interaction=1"
${padir}palinear \
    -p height.txt \
    -d ${gtdatadir}/test.prob.fvi \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 --ngp=2 --interaction=1 \
    -o height_ngp2_int1_fv


echo "Option --ngp=2 --robust"
${padir}palinear \
    -p height.txt \
    -d ${gtdatadir}/test.prob.fvi \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 --ngp=2 --robust \
    -o height_ngp2_robust_fv


echo "Option --ngp=2 --robust --interaction=1"
${padir}palinear \
    -p height.txt \
    -d ${gtdatadir}/test.prob.fvi \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 --ngp=2 --robust --interaction=1 \
    -o height_ngp2_robust_int1_fv
