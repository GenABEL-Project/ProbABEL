VERSION = 0.1-9e
SRCDIR = src
BINDIR = bin
DOCDIR = doc
LINREG = $(BINDIR)/palinear
LOGREG = $(BINDIR)/palogist
COXREG = $(BINDIR)/pacoxph
EXECUTABLES = $(LOGREG) $(LINREG) $(COXREG)

REGFILES =  $(SRCDIR)/data.h $(SRCDIR)/mematrix.h $(SRCDIR)/reg1.h $(SRCDIR)/usage.h $(SRCDIR)/main.cpp
COXBASE = $(SRCDIR)/chinv2 $(SRCDIR)/cholesky2 $(SRCDIR)/chsolve2 $(SRCDIR)/dmatrix $(SRCDIR)/coxfit2
COXSRC = $(COXBASE:=.c)
COXOBJ = $(COXBASE:=.o)

CPP = g++
CFLAGS = -I $(SRCDIR)/include -O2

all: $(EXECUTABLES) doc
	cp $(SRCDIR)/extIDS.pl $(BINDIR)/.
	cp $(SRCDIR)/prepare_data.R $(BINDIR)/.
	cp $(SRCDIR)/probabel.pl $(BINDIR)/probabel.pl_example
	cp $(SRCDIR)/probabel_config.cfg $(BINDIR)/probabel_config.cfg_example

$(LINREG): $(REGFILES)
	$(CPP) $(CFLAGS) -DLINEAR $(SRCDIR)/main.cpp $(SRCDIR)/fvlib/*.cpp -o $(LINREG)

$(LOGREG): $(REGFILES)
	$(CPP) $(CFLAGS) -DLOGISTIC $(SRCDIR)/main.cpp $(SRCDIR)/fvlib/*.cpp -o $(LOGREG)

$(COXREG): $(COXSRC) $(REGFILES)
	$(CPP) $(CFLAGS) -DCOXPH $(COXSRC) $(SRCDIR)/main.cpp $(SRCDIR)/fvlib/*.cpp -o $(COXREG)

doc:
	$(MAKE) -C $(DOCDIR)/

clean: clean_doc
	rm -f $(BINDIR)/* $(SRCDIR)/*~ $(SRCDIR)/*.o *.zip *.tar.gz examples/*.out.txt examples/*out
	$(MAKE) clean -C $(DOCDIR)/

clean_doc:
	$(MAKE) clean_doc -C $(DOCDIR)/

linux_distrib: clean
	cd .. ; tar -czvf ProbABEL_$(VERSION).tar.gz ProbABEL

win_distrib: all
	cd .. ; zip -r9 ProbABEL_$(VERSION)_win.zip ProbABEL

distrib: linux_distrib clean win_distrib clean

.PHONY:	distrib linux_distrib win_distrib clean clean_doc doc
