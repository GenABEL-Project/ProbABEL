## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = foreign color-tests

verified_results_dir = verified_results

verified_results = $(verified_results_dir)/height_add.out.txt	\
 $(verified_results_dir)/height_ngp2_add.out.txt		\
 $(verified_results_dir)/height_ngp2_2df.out.txt		\
 $(verified_results_dir)/height_ngp2_domin.out.txt		\
 $(verified_results_dir)/height_ngp2_over_domin.out.txt		\
 $(verified_results_dir)/height_ngp2_recess.out.txt		\
 $(verified_results_dir)/height_base_add.out.txt

input_files_dir = inputfiles

phenofiles = $(input_files_dir)/height.txt	\
 $(input_files_dir)/logist_data.txt		\
 $(input_files_dir)/coxph_data.txt		\
 $(input_files_dir)/coxph_data_nocovar.txt	\
 $(input_files_dir)/mmscore_pheno.PHE

genofiles = $(input_files_dir)/test.mldose		\
 $(input_files_dir)/test.mlprob				\
 $(input_files_dir)/test.mlinfo				\
 $(input_files_dir)/test.map				\
 $(input_files_dir)/test.dose.fvi			\
 $(input_files_dir)/test.dose.fvd			\
 $(input_files_dir)/test.prob.fvi			\
 $(input_files_dir)/test.prob.fvd			\
 $(input_files_dir)/mmscore_gen.mlinfo			\
 $(input_files_dir)/mmscore_gen.mldose			\
 $(input_files_dir)/mmscore_gen.mlprob			\
 $(input_files_dir)/mmscore_gen.dose.fvi		\
 $(input_files_dir)/mmscore_gen.dose.fvd		\
 $(input_files_dir)/mmscore_gen.prob.fvi		\
 $(input_files_dir)/mmscore_gen.prob.fvd		\
 $(input_files_dir)/mmscore_InvSigma_aj.sex.age.dat

input_files = $(phenofiles) $(genofiles) run_diff.sh test_all.sh

check_SCRIPTS =


if BUILD_palinear
check_SCRIPTS += check_probabel_chunk.sh check_dose_input.sh
check_SCRIPTS += test_qt.sh test_mms.sh test_flipmaf_palinear.sh
endif
if BUILD_palogist
check_SCRIPTS += test_bt.sh test_flipmaf_palogist.sh
endif
if BUILD_pacoxph
check_SCRIPTS += test_cox.sh test_cox_nocovar.sh test_flipmaf_pacoxph.sh
endif

AM_TESTS_ENVIRONMENT= \
	AWK=$(AWK) \
	SED=$(SED) \
	PA_BINDIR=$(top_builddir)/src
TESTS = $(check_SCRIPTS)

## Dependencies for the test scripts
TEST_EXTENSIONS=.sh
test_flipmaf_palinear.log: test_qt.log
test_flipmaf_palogist.log: test_bt.log
test_flipmaf_pacoxph.log: test_cox.log

## The pacoxph nocovar test doesn't run correctly until bug #1266 is
## fixed.
XFAIL_TESTS = check_pacoxph_nocovar.sh

EXTRA_DIST = $(verified_results) $(input_files) $(check_SCRIPTS)


##################################################
# Files that need to be cleaned up after testing
##################################################

cleanfiles_bt = logist_add.out.txt logist_fv_add.out.txt	\
logist_prob_2df.out.txt logist_prob_fv_2df.out.txt		\
logist_prob_fv_domin.out.txt logist_prob_add.out.txt		\
logist_prob_fv_add.out.txt logist_prob_over_domin.out.txt	\
logist_prob_domin.out.txt logist_prob_fv_recess.out.txt		\
logist_prob_fv_over_domin.out.txt logist_prob_recess.out.txt

cleanfiles_bt_flipmaf = logist_flipmaf_add.out.txt	\
logist_flipmaf_fv_add.out.txt				\
logist_prob_flipmaf_2df.out.txt				\
logist_prob_flipmaf_fv_2df.out.txt			\
logist_prob_flipmaf_add.out.txt				\
logist_prob_flipmaf_fv_add.out.txt			\
logist_prob_flipmaf_domin.out.txt			\
logist_prob_flipmaf_fv_domin.out.txt			\
logist_prob_flipmaf_over_domin.out.txt			\
logist_prob_flipmaf_fv_over_domin.out.txt		\
logist_prob_flipmaf_recess.out.txt			\
logist_prob_flipmaf_fv_recess.out.txt

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

cleanfiles_qt_flipmaf = linear_flipmaf_add.out.txt	\
linear_flipmaf_fv_add.out.txt				\
linear_ngp2_flipmaf_2df.out.txt				\
linear_ngp2_flipmaf_add.out.txt				\
linear_ngp2_flipmaf_domin.out.txt			\
linear_ngp2_flipmaf_over_domin.out.txt			\
linear_ngp2_flipmaf_recess.out.txt			\
linear_ngp2_flipmaf_fv_2df.out.txt			\
linear_ngp2_flipmaf_fv_add.out.txt			\
linear_ngp2_flipmaf_fv_domin.out.txt			\
linear_ngp2_flipmaf_fv_over_domin.out.txt		\
linear_ngp2_flipmaf_fv_recess.out.txt

cleanfiles_mms = mmscore_dose_add.out.txt mmscore_dose_fv_add.out.txt	\
mmscore_prob_fv_over_domin.out.txt mmscore_prob_fv_domin.out.txt	\
mmscore_prob_over_domin.out.txt mmscore_prob_fv_add.out.txt		\
mmscore_prob_fv_recess.out.txt mmscore_prob_domin.out.txt		\
mmscore_prob_recess.out.txt mmscore_prob_add.out.txt			\
mmscore_prob_2df.out.txt mmscore_prob_fv_2df.out.txt

cleanfiles_cox = pacoxph_nocovar_add.out.txt			\
coxph_data.txt							\
coxph_dose_add.out.txt coxph_dose_fv_add.out.txt		\
coxph_prob_fv_over_domin.out.txt coxph_prob_fv_domin.out.txt	\
coxph_prob_over_domin.out.txt coxph_prob_fv_add.out.txt		\
coxph_prob_fv_recess.out.txt coxph_prob_domin.out.txt		\
coxph_prob_recess.out.txt coxph_prob_add.out.txt		\
coxph_prob_2df.out.txt coxph_prob_fv_2df.out.txt		\
coxph_prob_allcov_recess.out.txt				\
coxph_prob_allcov_fv_2df.out.txt				\
coxph_prob_allcov_add.out.txt					\
coxph_prob_flipmaf_allcov_over_domin.out.txt			\
coxph_prob_allcov_fv_over_domin.out.txt				\
coxph_prob_flipmaf_allcov_fv_domin.out.txt			\
coxph_prob_flipmaf_allcov_fv_add.out.txt			\
coxph_prob_flipmaf_allcov_add.out.txt				\
coxph_prob_allcov_fv_domin.out.txt				\
coxph_prob_allcov_fv_recess.out.txt

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

cleanfiles_cox_flipmaf = coxph_flipmaf_add.out.txt	\
coxph_flipmaf_fv_add.out.txt				\
coxph_prob_flipmaf_2df.out.txt				\
coxph_prob_flipmaf_add.out.txt				\
coxph_prob_flipmaf_domin.out.txt			\
coxph_prob_flipmaf_over_domin.out.txt			\
coxph_prob_flipmaf_recess.out.txt			\
coxph_prob_flipmaf_fv_2df.out.txt			\
coxph_prob_flipmaf_fv_add.out.txt			\
coxph_prob_flipmaf_fv_domin.out.txt			\
coxph_prob_flipmaf_fv_over_domin.out.txt		\
coxph_prob_flipmaf_fv_recess.out.txt			\
coxph_prob_allcov_over_domin.out.txt			\
coxph_prob_allcov_fv_add.out.txt			\
coxph_prob_allcov_2df.out.txt				\
coxph_prob_flipmaf_allcov_fv_over_domin.out.txt		\
coxph_prob_flipmaf_allcov_2df.out.txt			\
coxph_prob_flipmaf_allcov_fv_2df.out.txt		\
coxph_prob_flipmaf_allcov_fv_recess.out.txt		\
coxph_prob_allcov_domin.out.txt				\
coxph_prob_flipmaf_allcov_recess.out.txt		\
coxph_prob_flipmaf_allcov_domin.out.txt


dose_files = chr1.chunk1.dose chr1.chunk2.dose chr2.chunk1.dose	\
 chr2.chunk2.dose chr1.dose chr2.dose

prob_files = $(dose_files:.dose=.prob)
map_files  = $(dose_files:.dose=.map)
info_files = $(dose_files:.dose=.info)

output_files = height_add.out.txt height_2df.out.txt			\
 height_domin.out.txt height_over_domin.out.txt height_recess.out.txt	\
 height_ngp2_recess.out.txt height_ngp2_domin.out.txt			\
 height_ngp2_over_domin.out.txt height_ngp2_2df.out.txt			\
 height_ngp2_add.out.txt height_base_add.out.txt

other_files = probabel probabel_config.cfg

cleanfiles_probabel_check = height.PHE $(dose_files) $(prob_files)	\
 $(map_files) $(info_files) $(output_files) $(other_files)

cleanfiles_dose_check = test.mldose

CLEANFILES = $(cleanfiles_probabel_check) $(cleanfiles_dose_check)	\
 $(cleanfiles_bt) $(cleanfiles_bt_flipmaf) $(cleanfiles_qt)		\
 $(cleanfiles_qt_flipmaf) $(cleanfiles_mms) $(cleanfiles_cox)		\
 $(cleanfiles_cox_nocovar) $(cleanfiles_cox_flipmaf)
