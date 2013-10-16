## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = foreign

phenofiles = height.txt				\
 logist_data.txt				\
 mmscore_pheno.PHE				\
 coxph_data.txt

gtdatadir = gtdata
genofiles = $(gtdatadir)/test.mldose		\
 $(gtdatadir)/test.mlprob			\
 $(gtdatadir)/test.mlinfo			\
 $(gtdatadir)/test.map				\
 $(gtdatadir)/test.dose.fvi			\
 $(gtdatadir)/test.dose.fvd			\
 $(gtdatadir)/test.prob.fvi			\
 $(gtdatadir)/test.prob.fvd			\
 $(gtdatadir)/mmscore_gen.mlinfo		\
 $(gtdatadir)/mmscore_gen.mldose		\
 $(gtdatadir)/mmscore_gen.mlprob		\
 $(gtdatadir)/mmscore_gen.dose.fvi		\
 $(gtdatadir)/mmscore_gen.dose.fvd		\
 $(gtdatadir)/mmscore_gen.prob.fvi		\
 $(gtdatadir)/mmscore_gen.prob.fvd


additional_files = mmscore.R			\
 mmscore_InvSigma_aj.sex.age.dat		\
 README

example_scripts = example_qt.sh			\
 example_bt.sh					\
 example_cox.sh					\
 example_mmscore.sh

exampledir = $(pkgdatadir)/examples
dist_example_DATA = $(phenofiles) $(additional_files)
dist_example_SCRIPTS = $(example_scripts)

gdatadir = $(pkgdatadir)/examples/gtdata
dist_gdata_DATA = $(genofiles)


cleanfiles_example_qt = height_allcov_fv_add.out.txt	\
 height_base_add.out.txt				\
 height_base_fv_add.out.txt				\
 height_int1_fv_add.out.txt				\
 height_ngp2_2df.out.txt				\
 height_ngp2_add.out.txt				\
 height_ngp2_allcov_fv_2df.out.txt			\
 height_ngp2_allcov_fv_add.out.txt			\
 height_ngp2_allcov_fv_domin.out.txt			\
 height_ngp2_allcov_fv_over_domin.out.txt		\
 height_ngp2_allcov_fv_recess.out.txt			\
 height_ngp2_domin.out.txt				\
 height_ngp2_fv_2df.out.txt				\
 height_ngp2_fv_add.out.txt				\
 height_ngp2_fv_domin.out.txt				\
 height_ngp2_fv_over_domin.out.txt			\
 height_ngp2_fv_recess.out.txt				\
 height_ngp2_int1_fv_2df.out.txt			\
 height_ngp2_int1_fv_add.out.txt			\
 height_ngp2_int1_fv_domin.out.txt			\
 height_ngp2_int1_fv_over_domin.out.txt			\
 height_ngp2_int1_fv_recess.out.txt			\
 height_ngp2_over_domin.out.txt				\
 height_ngp2_recess.out.txt				\
 height_ngp2_robust_fv_2df.out.txt			\
 height_ngp2_robust_fv_add.out.txt			\
 height_ngp2_robust_fv_domin.out.txt			\
 height_ngp2_robust_fv_over_domin.out.txt		\
 height_ngp2_robust_fv_recess.out.txt			\
 height_ngp2_robust_int1_fv_2df.out.txt			\
 height_ngp2_robust_int1_fv_add.out.txt			\
 height_ngp2_robust_int1_fv_domin.out.txt		\
 height_ngp2_robust_int1_fv_over_domin.out.txt		\
 height_ngp2_robust_int1_fv_recess.out.txt		\
 height_robust_fv_add.out.txt				\
 height_robust_int1_fv_add.out.txt

cleanfiles_example_bt = logist_add.out.txt	\
 logist_fv_add.out.txt				\
 logist_prob_2df.out.txt			\
 logist_prob_add.out.txt			\
 logist_prob_domin.out.txt			\
 logist_prob_fv_2df.out.txt			\
 logist_prob_fv_add.out.txt			\
 logist_prob_fv_domin.out.txt			\
 logist_prob_fv_over_domin.out.txt		\
 logist_prob_fv_recess.out.txt			\
 logist_prob_over_domin.out.txt			\
 logist_prob_recess.out.txt

cleanfiles_example_cox = coxph_dose_add.out.txt	\
 coxph_dose_fv_add.out.txt			\
 coxph_prob_2df.out.txt				\
 coxph_prob_add.out.txt				\
 coxph_prob_domin.out.txt			\
 coxph_prob_fv_2df.out.txt			\
 coxph_prob_fv_add.out.txt			\
 coxph_prob_fv_domin.out.txt			\
 coxph_prob_fv_over_domin.out.txt		\
 coxph_prob_fv_recess.out.txt			\
 coxph_prob_over_domin.out.txt			\
 coxph_prob_recess.out.txt

cleanfiles_example_mmscore = mmscore_dose_add.out.txt	\
 mmscore_dose_fv_add.out.txt				\
 mmscore_prob_2df.out.txt				\
 mmscore_prob_add.out.txt				\
 mmscore_prob_domin.out.txt				\
 mmscore_prob_fv_2df.out.txt				\
 mmscore_prob_fv_add.out.txt				\
 mmscore_prob_fv_domin.out.txt				\
 mmscore_prob_fv_over_domin.out.txt			\
 mmscore_prob_fv_recess.out.txt				\
 mmscore_prob_over_domin.out.txt			\
 mmscore_prob_recess.out.txt


CLEANFILES =  $(cleanfiles_example_qt) $(cleanfiles_example_bt)	\
 $(cleanfiles_example_cox) $(cleanfiles_example_mmscore)
