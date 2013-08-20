#!/bin/bash
#
# This script runs the R-based tests for logistic regression

Rcommand="R --vanilla --slave"

if [ -z ${srcdir} ]; then
    srcdir="."
fi

$Rcommand -f ${srcdir}/run_models_in_R_palogist.R --args ${srcdir}/
