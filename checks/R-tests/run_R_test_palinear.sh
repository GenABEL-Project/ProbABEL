#!/bin/bash
#
# This script runs the R-based tests for linear regression

Rcommand="R --vanilla --slave"

if [ -z ${srcdir} ]; then
    srcdir="."
fi
echo "SRCDIR=$srcdir"
$Rcommand -f ${srcdir}/run_models_in_R_palinear.R --args ${srcdir}/
