#!/bin/bash
# L.C. Karssen
# This script is used to test whether the Cox PH regression works
# correctly when no covariate is present in the phenotype input file.
# Currently this test fails, see bug #1266.

echo "Testing pacoxph without covariates..."

if [ -z ${AWK} ]; then
    AWK=awk
fi

# Exit with error when one of the steps in the script fails
set -e

# -------- Set some default paths and file names -------
if [ -z ${srcdir} ]; then
    srcdir="."
fi
inputdir="${srcdir}/inputfiles/"
padir="${srcdir}/../src/"

dosefile="$inputdir/test.mldose"
infofile="$inputdir/test.mlinfo"
mapfile="$inputdir/test.map"
orig_phenofile="$inputdir/coxph_data.txt"
phenofile="coxph_data.txt"
outfile="pacoxph_nocovar"
pacoxph="${padir}/pacoxph"

# ------ Prepare the phenotype file by removing the covariate column
# from the existing phenotype file ------
${AWK} '{print $1, $2, $3}' $orig_phenofile > $phenofile


# ---------- function definitions ----------
run_test ()
{
    ## When bug #1266 is fixed, this function should be expanded to
    ## include a verification against known-good results.
    echo "Checking whether Cox PH regression works without covariates..."
    $pacoxph -p $phenofile -d $dosefile -i $infofile -m $mapfile \
        -o $outfile

}


# ---------- Continuation of the main function ----------
run_test
echo "---- Finished check for Cox regression without covariates ----"
