
VERSION = 0.1-3-pacoxph-fix-1
BUILDDIRLNK = ProbABEL-$(VERSION)
SRCDIR = src
BINDIR = bin
DOCDIR = doc
LINREG = $(BINDIR)/palinear
LOGREG = $(BINDIR)/palogist
COXREG = $(BINDIR)/pacoxph
EXECUTABLES = $(COXREG)
## in this tree only build pacoxph. Leave these out:  $(LOGREG) $(LINREG)

REGFILES =  $(SRCDIR)/data.h $(SRCDIR)/mematrix.h $(SRCDIR)/reg1.h $(SRCDIR)/usage.h $(SRCDIR)/main.cpp
COXBASE = $(SRCDIR)/chinv2 $(SRCDIR)/cholesky2 $(SRCDIR)/chsolve2 $(SRCDIR)/dmatrix $(SRCDIR)/coxfit2 $(SRCDIR)/coxscore
COXSRC = $(COXBASE:=.c)
COXOBJ = $(COXBASE:=.o)

CPP = g++
CFLAGS = -I $(SRCDIR)/include -O3

all: $(EXECUTABLES)
	cp $(SRCDIR)/extIDS.pl $(BINDIR)/.
	cp $(SRCDIR)/prepare_data.R $(BINDIR)/.
	cp $(SRCDIR)/probabel.pl $(BINDIR)/probabel.pl_example
	cp $(SRCDIR)/probabel_config.cfg $(BINDIR)/probabel_config.cfg_example

$(LINREG): $(REGFILES)
	$(CPP) $(CFLAGS) -DLINEAR $(SRCDIR)/main.cpp -o $(LINREG)

$(LOGREG): $(REGFILES)
	$(CPP) $(CFLAGS) -DLOGISTIC $(SRCDIR)/main.cpp -o $(LOGREG)

$(COXREG): $(COXSRC) $(REGFILES)
	$(CPP) $(CFLAGS) -DCOXPH $(COXSRC) $(SRCDIR)/main.cpp -o $(COXREG)

clean:
	rm -f $(BINDIR)/* $(SRCDIR)/*~ $(SRCDIR)/*.o $(DOCDIR)/*~ *.zip *.tar.gz
	if [ -h ../$(BUILDDIRLNK) ]; then rm ../$(BUILDDIRLNK); fi

linux_distrib: clean
	if [ ! -d ../$(BUILDDIRLNK) ]; then ln -s $(PWD) ../$(BUILDDIRLNK); fi
	cd .. ; tar --exclude-vcs -czhvf ProbABEL_$(VERSION).tar.gz $(BUILDDIRLNK)
	if [ -h ../$(BUILDDIRLNK) ]; then rm ../$(BUILDDIRLNK); fi

win_distrib: all
	if [ ! -d ../$(BUILDDIRLNK) ]; then ln -s $(PWD) ../$(BUILDDIRLNK); fi
	cd .. ; zip -r9 ProbABEL_$(VERSION)_win.zip $(BUILDDIRLNK)
	if [ -h ../$(BUILDDIRLNK) ]; then rm ../$(BUILDDIRLNK); fi

distrib: linux_distrib clean win_distrib clean
