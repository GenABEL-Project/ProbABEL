AM_CXXFLAGS = -I include -O2
AM_CFLAGS = -I include -O2
AM_LDFLAGS =

## Using wildcards in these lists doesn't work. Also GNU make's ($wildcard,) doesn't
## work. It gives warning message about portability, but in the end doesn't work,
## I tried :-).
REGFILES = data.h mematrix.h reg1.h usage.h main.cpp version.h mematri1.h

COXBASE = chinv2 cholesky2 chsolve2 dmatrix
COXSRC = $(COXBASE:=.c) coxfit2.c

FVSRC = fvlib/AbstractMatrix.cpp fvlib/CastUtils.cpp			\
 fvlib/convert_util.cpp fvlib/FileVector.cpp fvlib/FilteredMatrix.cpp	\
 fvlib/frutil.cpp fvlib/Logger.cpp fvlib/RealHandlerWrapper.cpp		\
 fvlib/ReusableFileHandle.cpp fvlib/Transposer.cpp
FVHEADERS = fvlib/AbstractMatrix.h fvlib/CastUtils.h fvlib/const.h	\
 fvlib/convert_util.h fvlib/FileVector.h fvlib/frutil.h			\
 fvlib/frversion.h fvlib/Logger.h fvlib/Transposer.h			\
 fvlib/ReusableFileHandle.h fvlib/RealHandlerWrapper.h			\
 fvlib/FilteredMatrix.h

bin_PROGRAMS = palinear palogist
#pacoxph

palinear_SOURCES = $(REGFILES) $(FVSRC) $(FVHEADERS)
palinear_CXXFLAGS = -DLINEAR

palogist_SOURCES = $(REGFILES) $(FVSRC) $(FVHEADERS)
palogist_CXXFLAGS = -DLOGISTIC

# pacoxph_SOURCES = $(REGFILES) $(FVSRC) $(FVHEADERS) $(COXSRC)
# pacoxph_CXXFLAGS = -DCOXPH -I include
# pacoxph_CFLAGS = -DCOXPH -I include

dist_bin_SCRIPTS = probabel.pl extIDS.pl

## Old Makefile rules
# palinear: $(REGFILES)
# 	$(CXX) $(CXXFLAGS) -DLINEAR main.cpp $(FVSRC) -o $@

# palogist: $(REGFILES)
# 	$(CPP) $(CFLAGS) -DLOGISTIC main.cpp fvlib/*.cpp -o $(LOGREG)

# pacoxph: $(COXSRC) $(REGFILES)
# 	$(CPP) $(CFLAGS) -DCOXPH $(COXSRC) main.cpp fvlib/*.cpp -o $(COXREG)
