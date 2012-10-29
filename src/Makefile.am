## Process this file with automake to produce Makefile.in

## Using wildcards in these lists doesn't work. Also GNU make's ($wildcard,) doesn't
## work. It gives warning message about portability, but in the end doesn't work,
## I tried :-).
REGFILES = data.h data.cpp gendata.h gendata.cpp mematrix.h mematri1.h	\
 comand_line_settings.h comand_line_settings.cpp reg1.h usage.h		\
 usage.cpp main.cpp utilities.h utilities.cpp phedata.h phedata.cpp	\
 cholesky.h cholesky.cpp regdata.h regdata.cpp coxph_data.h		\
 coxph_data.cpp maskedmatrix.cpp maskedmatrix.h eigen_mematrix.h 	\
 eigen_mematri1.h

## Filevector files:
FVSRC = fvlib/AbstractMatrix.cpp fvlib/CastUtils.cpp			\
 fvlib/convert_util.cpp fvlib/FileVector.cpp fvlib/FilteredMatrix.cpp	\
 fvlib/frutil.cpp fvlib/Logger.cpp fvlib/RealHandlerWrapper.cpp		\
 fvlib/ReusableFileHandle.cpp fvlib/Transposer.cpp
FVHEADERS = fvlib/AbstractMatrix.h fvlib/CastUtils.h fvlib/const.h	\
 fvlib/convert_util.h fvlib/FileVector.h fvlib/frutil.h			\
 fvlib/frversion.h fvlib/Logger.h fvlib/Transposer.h			\
 fvlib/ReusableFileHandle.h fvlib/RealHandlerWrapper.h			\
 fvlib/FilteredMatrix.h

## R headers, used by pacoxph
RHEADERS = include/R.h include/Rmath.h include/R_ext/Arith.h		\
 include/R_ext/Boolean.h include/Rconfig.h include/R_ext/Complex.h	\
 include/R_ext/Constants.h include/R_ext/Error.h			\
 include/R_ext/libextern.h include/R_ext/Memory.h			\
 include/R_ext/Print.h include/R_ext/Random.h include/R_ext/Utils.h	\
 include/R_ext/RS.h

AM_CXXFLAGS = -Wall

bin_PROGRAMS = palinear palogist
if BUILD_pacoxph
## Extra files for pacoxph
COXBASE = coxfit2 chinv2 cholesky2 chsolve2 dmatrix
COXSRC = $(COXBASE:=.c)

bin_PROGRAMS += pacoxph
endif

palinear_SOURCES = $(REGFILES) $(FVSRC) $(FVHEADERS)
palinear_CXXFLAGS = -DLINEAR $(AM_CXXFLAGS)

palogist_SOURCES = $(REGFILES) $(FVSRC) $(FVHEADERS)
palogist_CXXFLAGS = -DLOGISTIC $(AM_CXXFLAGS)

pacoxph_SOURCES = $(COXSRC) $(REGFILES) $(FVSRC) $(FVHEADERS)	\
 $(RHEADERS) survS.h survproto.h
pacoxph_CXXFLAGS = -DCOXPH -I $(top_srcdir)/src/include $(AM_CXXFLAGS)
pacoxph_CFLAGS = -DCOXPH -I $(top_srcdir)/src/include $(AM_CFLAGS)

## Install these scripts in the bin directory as well:
dist_bin_SCRIPTS = probabel.pl extIDS.pl

## Install the config file
dist_sysconf_DATA = probabel_config.cfg.example

## Insert the correct values of several variables in the installed probabel.pl
install-exec-hook:
	$(SED) -i "s;probabel_config.cfg;@sysconfdir@/probabel_config.cfg;" \
		 $(DESTDIR)$(bindir)/probabel.pl
	$(SED) -i "s;\./;@bindir@/;" $(DESTDIR)$(bindir)/probabel.pl
	$(SED) -i "s;PROBABEL_VERSION;@VERSION@;" $(DESTDIR)$(bindir)/probabel.pl
