#!/bin/bash
#
# This script checks the flipmaf functionality of ProbABEL.

echo "Checking --flipmaf option with QT..."

scriptdir=$(dirname $0)

if [ -z ${PA_BINDIR} ]; then
    PA_BINDIR="${scriptdir}/../src/"
fi
if [ -z ${srcdir} ]; then
    srcdir="."
    PA_BINDIR=${scriptdir}/../src/
fi
if [ -z ${AWK} ]; then
    AWK=awk
fi

. ${scriptdir}/run_diff.sh

inputdir=${scriptdir}/inputfiles
pabin=${PA_BINDIR}/palinear

# Redirect all output to file descriptor 3 to /dev/null except if
# the first argument is "verbose" then redirect handle 3 to stdout
exec 3>/dev/null
if [ "$1" = "verbose" ]; then
    echo "Verbose mode ON"
    exec 3>&1
fi


verify_results() {
    # This function compares the outcomes of two files, one with a
    # allele flip, one without.
    # It requires two arguments:
    # $2: first model to compare
    # $3: second model to compare
    model1=$1
    model2=$2

    blanks="                                                                      "

    text="   QT check --flipmaf (${model1} vs. ${model2}): "
    if ${AWK} -f verify_flipmaf.awk \
              linear_ngp2_${model1}.out.txt \
              linear_ngp2_flipmaf_${model2}.out.txt 2>&1 |
            grep -w "^Discordant" ; then
        echo -e "${text}${blanks:${#text}} FAILED"
        exit 1
    else
        echo -e "${text}${blanks:${#text}} OK"
    fi
}

echo "Option --flipmaf"
$pabin -p ${inputdir}/height.txt \
    -d ${inputdir}/test.mldose \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --flipmaf \
    -o linear_flipmaf \
    >& 3
$pabin \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.dose.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --flipmaf \
    -o linear_flipmaf_fv \
    >& 3

run_diff linear_flipmaf_add.out.txt \
    linear_flipmaf_fv_add.out.txt \
    "   QT check: flipmaf: dose vs. dose_fv"

echo "Check signs of beta_SNP, size of SE_beta"
${AWK} -f verify_flipmaf.awk \
       linear_base_add.out.txt linear_flipmaf_add.out.txt


echo "Option --ngp=2 --flipmaf"
$pabin \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.mlprob \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --ngp=2 --flipmaf \
    -o linear_ngp2_flipmaf \
    >& 3
$pabin \
    -p ${inputdir}/height.txt \
    -d ${inputdir}/test.prob.fvi \
    -i ${inputdir}/test.mlinfo \
    -m ${inputdir}/test.map \
    -c 19 --ngp=2 --flipmaf \
    -o linear_ngp2_flipmaf_fv \
    >& 3

for model in add domin over_domin recess 2df; do
    run_diff linear_ngp2_flipmaf_${model}.out.txt \
        linear_ngp2_flipmaf_fv_${model}.out.txt \
        "   QT check --flipmaf ($model model): prob vs. prob_fv"
done

echo "Check signs of beta_SNP, size of SE_beta"
# After flipping alleles: beta_add -> -beta_add
verify_results add add
# After flipping alleles: beta_A1A1 -> -beta_A1A1
verify_results 2df 2df
# After flipping alleles: beta_domin -> -beta_recess
verify_results domin recess
# After flipping alleles: beta_recess -> -beta_domin
verify_results recess domin
# For the overdominant model we don't expect changes
run_diff linear_ngp2_over_domin.out.txt \
         linear_ngp2_flipmaf_over_domin.out.txt \
         "   QT check --flipmaf (overdomin): "
