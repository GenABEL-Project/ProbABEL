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
