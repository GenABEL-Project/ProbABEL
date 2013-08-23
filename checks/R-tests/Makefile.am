## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = foreign color-tests

R_test_files =  run_models_in_R_palinear.R	\
 run_models_in_R_palogist.R			\
 run_models_in_R_pacox.R			\
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
dist_tests_DATA = $(R_test_files) $(check_SCRIPTS)
#dist_tests_SCRIPTS = $(check_SCRIPTS)

TESTS_ENVIRONMENT = bash



TESTS = $(check_SCRIPTS)
