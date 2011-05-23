MANNAME = ProbABEL_manual
MANTEXSRC = $(MANNAME).tex
MANAUX = $(patsubst %.tex, %.aux, $(MANTEXSRC))
MANPDF =  $(patsubst %.tex, %.pdf, $(MANTEXSRC))

.PHONY:	clean clean_doc

$(MANPDF): $(MANTEXSRC) $(MANAUX)
	@echo === Making final PDF ===
	pdflatex $(MANTEXSRC)
	@while ( grep "Rerun to " 			\
		$(MANNAME).log > /dev/null ); do		\
		echo '** Re-running LaTeX **';		\
		pdflatex $(MANTEXSRC);		\
	done
	pdflatex $(MANTEXSRC)


# A target needed to keep track of the nr. of LaTeX runs
.btmp:


clean_doc:
	rm -f *~ *.log *.idx *.out *.toc *.aux
	rm -rf auto

clean: clean_doc
	rm -f $(MANPDF)
	rm -f .btmp

%.aux:	%.tex $(MANSRC)
	@echo === Making $< file ===
	pdflatex $(MANTEXSRC)
# Look for citations. Make sure grep never returns an error code.
	@grep "^\\\\citation" $(MANAUX) > .btmp.new || true
# If the citations are not changed, don't do anything. Otherwise replace
# the .btmp file to make sure Bibtex will be run.
	@if ( diff .btmp.new .btmp  > /dev/null ); then \
		rm .btmp.new; \
	else \
		mv .btmp.new .btmp; \
	fi

# $(DOCS): $()
# 	cd $(DOCDIR); pdflatex $(DOCSRC)
