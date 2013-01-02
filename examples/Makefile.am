## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = foreign color-tests

examplefiles = height_ngp2_robust_int1_over_domin.out.txt.save		\
 test.prob.fvi mmscore_gen.prob.fvd test.map mldose.IDS			\
 impute_with_missing.mlinfo						\
 height_ngp2_robust_int1_recess.out.txt.save				\
 height_ngp2_2df.out.txt.save mmscore_gen.mlprob mmscore_gen.mlinfo	\
 height_ngp2_robust_recess.out.txt.save example_cox.sh test.mlprob	\
 test.prob.fvd impute_with_missing_robust_int1_add.out.txt.save		\
 test.mlinfo regression_add.out.txt.save				\
 height_ngp2_robust_add.out.txt.save mmscore_gen.dose.fvi		\
 impute_with_missing_allcov_add.out.txt.save example_bt.sh		\
 height_ngp2_int1_recess.out.txt.save example_qt.sh coxph_data.txt	\
 height_ngp2_robust_int1_domin.out.txt.save				\
 height_ngp2_allcov_2df.out.txt.save					\
 impute_with_missing_int1_add.out.txt.save mmscore_gen.mldose		\
 example_mms.sh height_ngp2_recess.out.txt.save				\
 height_ngp2_over_domin.out.txt.save					\
 height_ngp2_int1_domin.out.txt.save					\
 height_ngp2_allcov_domin.out.txt.save					\
 height_ngp2_robust_over_domin.out.txt.save				\
 height_ngp2_robust_int1_add.out.txt.save height_int1_add.out.txt.save	\
 height_allcov_add.out.txt.save						\
 impute_with_missing_base_add.out.txt.save example_all.sh		\
 height_base_add.out.txt.save height_ngp2_int1_over_domin.out.txt.save	\
 height_robust_add.out.txt.save test.mldose				\
 height_ngp2_allcov_add.out.txt.save mmscore_InvSigma_aj.sex.age.dat	\
 height_ngp2_robust_int1_2df.out.txt.save mmscore_gen.dose.fvd		\
 mmscore.R test.dose.fvd height_ngp2_domin.out.txt.save test.dose.fvi	\
 impute_with_missing_robust_add.out.txt.save allheight.txt		\
 mmscore_add.out.txt.save coxph.out.txt_add.out.txt.save		\
 mmscore_pheno.PHE height_ngp2_robust_2df.out.txt.save			\
 height_ngp2_int1_add.out.txt.save mmscore_gen.prob.fvi			\
 height_robust_int1_add.out.txt.save logist.out.txt_add.out.txt.save	\
 height_ngp2_allcov_recess.out.txt.save impute_with_missing.PHE		\
 height_ngp2_allcov_over_domin.out.txt.save				\
 height_ngp2_add.out.txt.save logist_data.txt				\
 height_ngp2_int1_2df.out.txt.save					\
 impute_with_missing_add.out.txt.save					\
 height_ngp2_robust_domin.out.txt.save height.txt

exampledir = $(pkgdatadir)/examples
dist_example_DATA = $(examplefiles)

TESTS_ENVIRONMENT = sh
check_SCRIPTS = example_bt.sh example_qt.sh example_mms.sh
# example_cox.sh: don't run this test. It takes several hours...

TESTS = $(check_SCRIPTS)

cleanfiles_bt = logist_add.out.txt logist_fv_add.out.txt	\
logist_prob_2df.out.txt logist_prob_fv_2df.out.txt		\
logist_prob_fv_domin.out.txt logist_prob_add.out.txt		\
logist_prob_fv_add.out.txt logist_prob_over_domin.out.txt	\
logist_prob_domin.out.txt logist_prob_fv_recess.out.txt		\
logist_prob_fv_over_domin.out.txt logist_prob_recess.out.txt

cleanfiles_qt = height_base_add.out.txt height_base_fv_add.out.txt	\
height_allcov_add.out.txt height_allcov_fv_add.out.txt			\
height_int1_add.out.txt height_int1_fv_add.out.txt			\
height_robust_add.out.txt height_robust_fv_add.out.txt			\
height_robust_int1_add.out.txt height_robust_int1_fv_add.out.txt	\
height_ngp2_add.out.txt height_ngp2_fv_add.out.txt			\
height_ngp2_domin.out.txt height_ngp2_fv_domin.out.txt			\
height_ngp2_over_domin.out.txt height_ngp2_fv_over_domin.out.txt	\
height_ngp2_recess.out.txt height_ngp2_fv_recess.out.txt		\
height_ngp2_2df.out.txt height_ngp2_fv_2df.out.txt			\
height_ngp2_allcov_add.out.txt height_ngp2_allcov_fv_add.out.txt	\
height_ngp2_allcov_domin.out.txt height_ngp2_allcov_fv_domin.out.txt	\
height_ngp2_allcov_over_domin.out.txt					\
height_ngp2_allcov_fv_over_domin.out.txt				\
height_ngp2_allcov_recess.out.txt height_ngp2_allcov_fv_recess.out.txt	\
height_ngp2_allcov_recess.out.txt height_ngp2_allcov_fv_recess.out.txt	\
height_ngp2_int1_add.out.txt height_ngp2_int1_fv_add.out.txt		\
height_ngp2_int1_domin.out.txt height_ngp2_int1_fv_domin.out.txt	\
height_ngp2_int1_over_domin.out.txt					\
height_ngp2_int1_fv_over_domin.out.txt height_ngp2_int1_recess.out.txt	\
height_ngp2_int1_fv_recess.out.txt height_ngp2_int1_2df.out.txt		\
height_ngp2_int1_fv_2df.out.txt height_ngp2_robust_add.out.txt		\
height_ngp2_robust_fv_add.out.txt height_ngp2_robust_domin.out.txt	\
height_ngp2_robust_fv_domin.out.txt					\
height_ngp2_robust_over_domin.out.txt					\
height_ngp2_robust_fv_over_domin.out.txt				\
height_ngp2_robust_recess.out.txt height_ngp2_robust_fv_recess.out.txt	\
height_ngp2_robust_2df.out.txt height_ngp2_robust_fv_2df.out.txt	\
height_ngp2_robust_int1_add.out.txt					\
height_ngp2_robust_int1_fv_add.out.txt					\
height_ngp2_robust_int1_domin.out.txt					\
height_ngp2_robust_int1_fv_domin.out.txt				\
height_ngp2_robust_int1_over_domin.out.txt				\
height_ngp2_robust_int1_fv_over_domin.out.txt				\
height_ngp2_robust_int1_recess.out.txt					\
height_ngp2_robust_int1_fv_recess.out.txt				\
height_ngp2_robust_int1_2df.out.txt					\
height_ngp2_robust_int1_fv_2df.out.txt					\
height_ngp2_allcov_fv_2df.out.txt height_ngp2_robust_fv_add.out.txt	\
height_ngp2_allcov_2df.out.txt

cleanfiles_mms = mmscore_add.out.txt mmscore_fv_add.out.txt		\
mmscore_prob_fv_over_domin.out.txt mmscore_prob_fv_domin.out.txt	\
mmscore_prob_over_domin.out.txt mmscore_prob_fv_add.out.txt		\
mmscore_prob_fv_recess.out.txt mmscore_prob_domin.out.txt		\
mmscore_prob_recess.out.txt mmscore_prob_add.out.txt			\
mmscore_prob_2df.out.txt mmscore_prob_fv_2df.out.txt

cleanfiles_cox = coxph_add.out.txt coxph_fv_add.out.txt

CLEANFILES = $(cleanfiles_bt) $(cleanfiles_qt) $(cleanfiles_mms) $(cleanfiles_cox)
