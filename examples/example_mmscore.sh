#!/bin/bash
# This script shows some examples on how to run ProbABEL's palinear
# module for quantitative traits combined with the mmscore option. See
# also the mmscore.R file on how to prepare the phenotype and
# inverse-sigma files.

# Set the name of the directory with the genotype data files
gtdatadir="gtdata"

# Set the path to the executables (palinear/palogist/pacoxph). If you
# installed ProbABEL according to the instructions in the doc/INSTALL
# file, or installed it via your Linux distribution's package manager,
# you can leave this variable empty.
padir="../src/"


# Using text-based dosage genotype files as input
${padir}palinear \
    -p mmscore_pheno.PHE \
    -i ${gtdatadir}/mmscore_gen.mlinfo \
    -d ${gtdatadir}/mmscore_gen.mldose \
    --sep="," \
    -o mmscore_dose \
    --mmscore mmscore_InvSigma_aj.sex.age.dat


# Using filevector (DatABEL) files as dosage genotype input
${padir}palinear \
    -p mmscore_pheno.PHE \
    -i ${gtdatadir}/mmscore_gen.mlinfo \
    -d ${gtdatadir}/mmscore_gen.dose.fvi \
    --sep="," \
    -o mmscore_dose_fv \
    --mmscore mmscore_InvSigma_aj.sex.age.dat


# Using text-based probability files as genotype input
${padir}palinear \
    -p mmscore_pheno.PHE \
    -i ${gtdatadir}/mmscore_gen.mlinfo \
    -d ${gtdatadir}/mmscore_gen.mlprob \
    --ngpreds=2 --sep="," \
    -o mmscore_prob \
    --mmscore mmscore_InvSigma_aj.sex.age.dat


# Using filevector (DatABEL) probability files as genotype input
${padir}palinear \
    -p mmscore_pheno.PHE \
    -i ${gtdatadir}/mmscore_gen.mlinfo \
    -d ${gtdatadir}/mmscore_gen.prob.fvi \
    --ngpreds=2 --sep="," \
    -o mmscore_prob_fv \
    --mmscore mmscore_InvSigma_aj.sex.age.dat
