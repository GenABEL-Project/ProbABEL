## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = foreign color-tests

testsfiles = check_probabel.pl_chunk.sh
known_results=known_good_results/height_add.out.txt
testsdir = $(pkgdatadir)/tests
dist_tests_DATA = $(testsfiles) $(known_results)

TESTS_ENVIRONMENT = sh
check_SCRIPTS = check_probabel.pl_chunk.sh

TESTS = $(check_SCRIPTS)

chunk_files=chunk1.chr1.dose chunk2.chr1.dose chunk1.chr2.dose		\
 chunk2.chr2.dose chunk1.chr1.map chunk2.chr1.map chunk1.chr2.map	\
 chunk2.chr2.map chunk1.chr1.info chunk2.chr1.info chunk1.chr2.info	\
 chunk2.chr2.info chr1.dose chr2.dose chr1.info chr2.info chr1.map	\
 chr2.map probabel.pl probabel_config.cfg height.PHE			\
 height_add.out.txt


CLEANFILES = $(chunk_files)
