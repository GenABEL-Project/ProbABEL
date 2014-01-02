## Process this file with automake to produce Makefile.in

dist_doc_DATA = ProbABEL_manual.tex QuickStart.txt			\
 short_coxph_data.txt short_height_base_add.out.txt.save		\
 short_height.txt short_logist_data.txt short_test.mldose test.map	\
 test.mlinfo test_regression.R COPYING LICENSE INSTALL ChangeLog	\
 TODO

man_MANS = probabel.1 probabel.pl.1
if BUILD_palinear
man_MANS += palinear.1
endif
if BUILD_palogist
man_MANS += palogist.1
endif
if BUILD_pacoxph
man_MANS += pacoxph.1
endif

EXTRA_DIST = $(man_MANS)

if HAVE_PDFLATEX
MANNAME = ProbABEL_manual
MANPDF = $(MANNAME).pdf
MANTEXSRC = $(MANNAME).tex
MANAUX = $(MANNAME).aux
doc_DATA = $(MANPDF)

CLEANFILES = $(MANPDF) $(MANNAME).log $(MANNAME).idx $(MANNAME).out \
 $(MANNAME).toc $(MANNAME).idx $(MANNAME).ilg $(MANNAME).ind $(MANAUX) .btmp

endif


# Several make rules to generate the PDF from the LaTeX source
.aux.pdf:
	@echo === Making PDF ===
	$(PDFLATEX) $(srcdir)/$(MANTEXSRC)
	@while ( grep "Rerun to "			\
		$(MANNAME).log  ); do			\
		echo '** Re-running LaTeX **';		\
		$(PDFLATEX) $(srcdir)/$(MANTEXSRC);	\
	done
	@echo === Making index ===
	makeindex $(MANNAME)
	@echo === Making final PDF ===
	$(PDFLATEX) $(srcdir)/$(MANTEXSRC)

.tex.aux:
	@echo === Making $@ file ===
	$(PDFLATEX) $(srcdir)/$(MANTEXSRC)
## Look for citations. Make sure grep never returns an error code.
	@grep "^\\\\citation" $(MANAUX) > .btmp.new || true
## If the citations are not changed, don't do anything. Otherwise replace
## the .btmp file to make sure BibTeX will be run.
	@if ( diff .btmp.new .btmp  > /dev/null ); then \
		rm .btmp.new; \
	else \
		mv .btmp.new .btmp; \
	fi

# A target needed to keep track of the nr. of LaTeX runs
.btmp:
