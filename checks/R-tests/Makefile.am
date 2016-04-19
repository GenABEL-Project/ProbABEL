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
AM_TESTS_ENVIRONMENT= \
	PA_BINDIR=$(top_builddir)/src
TESTS = $(check_SCRIPTS)

## The palogist R test still doesn't run correctly.
XFAIL_TESTS = run_R_test_palogist.sh

EXTRA_DIST = $(check_SCRIPTS) $(R_test_files)


##################################################
# Files that need to be cleaned up after testing
##################################################

cleanfiles_bt = logist_add.out.txt logist_fv_add.out.txt	\
logist_prob_2df.out.txt logist_prob_fv_2df.out.txt		\
logist_prob_fv_domin.out.txt logist_prob_add.out.txt		\
logist_prob_fv_add.out.txt logist_prob_over_domin.out.txt	\
logist_prob_domin.out.txt logist_prob_fv_recess.out.txt		\
logist_prob_fv_over_domin.out.txt logist_prob_recess.out.txt

cleanfiles_qt = linear_base_add.out.txt linear_base_fv_add.out.txt	\
linear_allcov_add.out.txt linear_allcov_fv_add.out.txt			\
linear_int1_add.out.txt linear_int1_fv_add.out.txt			\
linear_robust_add.out.txt linear_robust_fv_add.out.txt			\
linear_robust_int1_add.out.txt linear_robust_int1_fv_add.out.txt	\
linear_ngp2_add.out.txt linear_ngp2_fv_add.out.txt			\
linear_ngp2_domin.out.txt linear_ngp2_fv_domin.out.txt			\
linear_ngp2_over_domin.out.txt linear_ngp2_fv_over_domin.out.txt	\
linear_ngp2_recess.out.txt linear_ngp2_fv_recess.out.txt		\
linear_ngp2_2df.out.txt linear_ngp2_fv_2df.out.txt			\
linear_ngp2_allcov_add.out.txt linear_ngp2_allcov_fv_add.out.txt	\
linear_ngp2_allcov_domin.out.txt linear_ngp2_allcov_fv_domin.out.txt	\
linear_ngp2_allcov_over_domin.out.txt					\
linear_ngp2_allcov_fv_over_domin.out.txt				\
linear_ngp2_allcov_recess.out.txt linear_ngp2_allcov_fv_recess.out.txt	\
linear_ngp2_allcov_recess.out.txt linear_ngp2_allcov_fv_recess.out.txt	\
linear_ngp2_int1_add.out.txt linear_ngp2_int1_fv_add.out.txt		\
linear_ngp2_int1_domin.out.txt linear_ngp2_int1_fv_domin.out.txt	\
linear_ngp2_int1_over_domin.out.txt					\
linear_ngp2_int1_fv_over_domin.out.txt linear_ngp2_int1_recess.out.txt	\
linear_ngp2_int1_fv_recess.out.txt linear_ngp2_int1_2df.out.txt		\
linear_ngp2_int1_fv_2df.out.txt linear_ngp2_robust_add.out.txt		\
linear_ngp2_robust_fv_add.out.txt linear_ngp2_robust_domin.out.txt	\
linear_ngp2_robust_fv_domin.out.txt					\
linear_ngp2_robust_over_domin.out.txt					\
linear_ngp2_robust_fv_over_domin.out.txt				\
linear_ngp2_robust_recess.out.txt linear_ngp2_robust_fv_recess.out.txt	\
linear_ngp2_robust_2df.out.txt linear_ngp2_robust_fv_2df.out.txt	\
linear_ngp2_robust_int1_add.out.txt					\
linear_ngp2_robust_int1_fv_add.out.txt					\
linear_ngp2_robust_int1_domin.out.txt					\
linear_ngp2_robust_int1_fv_domin.out.txt				\
linear_ngp2_robust_int1_over_domin.out.txt				\
linear_ngp2_robust_int1_fv_over_domin.out.txt				\
linear_ngp2_robust_int1_recess.out.txt					\
linear_ngp2_robust_int1_fv_recess.out.txt				\
linear_ngp2_robust_int1_2df.out.txt					\
linear_ngp2_robust_int1_fv_2df.out.txt					\
linear_ngp2_allcov_fv_2df.out.txt linear_ngp2_robust_fv_add.out.txt	\
linear_ngp2_allcov_2df.out.txt

cleanfiles_mms = mmscore_dose_add.out.txt mmscore_dose_fv_add.out.txt	\
mmscore_prob_fv_over_domin.out.txt mmscore_prob_fv_domin.out.txt	\
mmscore_prob_over_domin.out.txt mmscore_prob_fv_add.out.txt		\
mmscore_prob_fv_recess.out.txt mmscore_prob_domin.out.txt		\
mmscore_prob_recess.out.txt mmscore_prob_add.out.txt			\
mmscore_prob_2df.out.txt mmscore_prob_fv_2df.out.txt

cleanfiles_cox = coxph_dose_add.out.txt coxph_dose_fv_add.out.txt	\
coxph_prob_fv_over_domin.out.txt coxph_prob_fv_domin.out.txt		\
coxph_prob_over_domin.out.txt coxph_prob_fv_add.out.txt			\
coxph_prob_fv_recess.out.txt coxph_prob_domin.out.txt			\
coxph_prob_recess.out.txt coxph_prob_add.out.txt			\
coxph_prob_2df.out.txt coxph_prob_fv_2df.out.txt			\
coxph_prob_allcov_recess.out.txt					\
coxph_prob_allcov_fv_2df.out.txt					\
coxph_prob_allcov_2df.out.txt						\
coxph_prob_allcov_domin.out.txt						\
coxph_prob_allcov_fv_over_domin.out.txt					\
coxph_prob_allcov_over_domin.out.txt					\
coxph_prob_allcov_fv_add.out.txt					\
coxph_prob_allcov_fv_recess.out.txt					\
coxph_prob_allcov_add.out.txt						\
coxph_prob_allcov_fv_domin.out.txt

cleanfiles_cox_nocovar = coxph_dose_nocovar_add.out.txt	\
coxph_dose_nocovar_fv_add.out.txt			\
coxph_prob_nocovar_fv_over_domin.out.txt		\
coxph_prob_nocovar_fv_domin.out.txt			\
coxph_prob_nocovar_over_domin.out.txt			\
coxph_prob_nocovar_fv_add.out.txt			\
coxph_prob_nocovar_fv_recess.out.txt			\
coxph_prob_nocovar_domin.out.txt			\
coxph_prob_nocovar_recess.out.txt			\
coxph_prob_nocovar_add.out.txt				\
coxph_prob_nocovar_2df.out.txt				\
coxph_prob_nocovar_fv_2df.out.txt

CLEANFILES = $(cleanfiles_bt) $(cleanfiles_qt) $(cleanfiles_mms)	\
 $(cleanfiles_cox) $(cleanfiles_cox_nocovar)
