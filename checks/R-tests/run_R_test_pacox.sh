#!/bin/bash
#
# This script runs the R-based tests for Cox PH regression

Rcommand="R --vanilla --slave --quiet"

if [ -z ${srcdir} ]; then
    srcdir="."
fi

$Rcommand -f ${srcdir}/run_models_in_R_pacox.R --args ${srcdir}/
