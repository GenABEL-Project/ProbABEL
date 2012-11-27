## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = foreign color-tests

tests_files = check_probabel.pl_chunk.sh check_dose_input.sh
verified_results = verified_results/height_add.out.txt	\
 verified_results/height_ngp2_add.out.txt			\
 verified_results/height_ngp2_2df.out.txt			\
 verified_results/height_ngp2_domin.out.txt			\
 verified_results/height_ngp2_over_domin.out.txt		\
 verified_results/height_ngp2_recess.out.txt                    \
 verified_results/height_base_add.out.txt

dose_files = chr1.chunk1.dose chr1.chunk2.dose chr2.chunk1.dose	\
 chr2.chunk2.dose chr1.dose chr2.dose

prob_files = chr1.chunk1.prob chr1.chunk2.prob chr2.chunk1.prob	\
 chr2.chunk2.prob chr1.prob chr2.prob

map_files = chr1.chunk1.map chr1.chunk2.map chr2.chunk1.map	\
 chr2.chunk2.map chr1.map chr2.map

info_files = chr1.chunk1.info chr1.chunk2.info chr2.chunk1.info	\
 chr2.chunk2.info chr1.info chr2.info

output_files = height_add.out.txt height_2df.out.txt			\
 height_domin.out.txt height_over_domin.out.txt height_recess.out.txt	\
 height_ngp2_recess.out.txt height_ngp2_domin.out.txt			\
 height_ngp2_over_domin.out.txt height_ngp2_2df.out.txt			\
 height_ngp2_add.out.txt height_base_add.out.txt

other_files = probabel.pl probabel_config.cfg height.PHE test.mldose



testsdir = $(pkgdatadir)/tests
dist_tests_DATA = $(tests_files) $(verified_results)

TESTS_ENVIRONMENT = sh
check_SCRIPTS = $(tests_files)

TESTS = $(check_SCRIPTS)

CLEANFILES = $(dose_files) $(prob_files) $(map_files) ${info_files} \
 ${output_files} ${other_files}
