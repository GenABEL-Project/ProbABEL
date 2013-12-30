## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = foreign color-tests

R_test_files = run_models_in_R_palinear.R run_model_linear.R	\
 run_models_in_R_palogist.R run_model_logist.R			\
 run_models_in_R_pacox.R run_model_coxph.R			\
 initial_checks.R

check_SCRIPTS =

if HAVE_R
if BUILD_palinear
check_SCRIPTS += run_R_test_palinear.sh
endif
if BUILD_palogist
check_SCRIPTS += run_R_test_palogist.sh
endif
if BUILD_pacoxph
check_SCRIPTS += run_R_test_pacox.sh
endif
endif


testsdir = $(pkgdatadir)/checks/R-tests

TESTS = $(check_SCRIPTS)

## The palogist R test still doesn't run correctly.
XFAIL_TESTS = run_R_test_palogist.sh

EXTRA_DIST = $(check_SCRIPTS) $(R_test_files)
