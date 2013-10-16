#!/bin/bash
# This script shows some examples on how to run ProbABEL's palogist
# module for binary traits.

# Set the name of the directory with the genotype data files
gtdatadir="gtdata"

# Set the path to the executables (palinear/palogist/pacoxph). If you
# installed ProbABEL according to the instructions in the doc/INSTALL
# file, or installed it via your Linux distribution's package manager,
# you can leave this variable empty.
padir="../src/"


# Using text-based dosage genotype files as input
${padir}palogist \
    -p logist_data.txt \
    -d ${gtdatadir}/test.mldose \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 \
    -o logist


# Using filevector (DatABEL) files as dosage genotype input
${padir}palogist \
    -p logist_data.txt \
    -d ${gtdatadir}/test.dose.fvi \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    -c 19 \
    -o logist_fv


# Using text-based probability files as genotype input
${padir}palogist \
    -p logist_data.txt \
    -d ${gtdatadir}/test.mlprob \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o logist_prob


# Using filevector (DatABEL) probability files as genotype input
${padir}palogist \
    -p logist_data.txt \
    -d ${gtdatadir}/test.prob.fvi \
    -i ${gtdatadir}/test.mlinfo \
    -m ${gtdatadir}/test.map \
    --ngpreds=2 \
    -c 19 \
    -o logist_prob_fv
