## Process this file with automake to produce Makefile.in

## Using wildcards in these lists doesn't work. Also GNU make's ($wildcard,) doesn't
## work. It gives warning message about portability, but in the end doesn't work,
## I tried :-).
REGFILES = data.h data.cpp gendata.h gendata.cpp mematrix.h mematri1.h	\
 command_line_settings.h command_line_settings.cpp reg1.h usage.h		\
 usage.cpp main.cpp utilities.h utilities.cpp phedata.h phedata.cpp	\
 cholesky.h cholesky.cpp regdata.h regdata.cpp  \
 maskedmatrix.cpp maskedmatrix.h reg1.cpp main_functions_dump.h main_functions_dump.cpp

EIGENFILES = eigen_mematrix.h eigen_mematrix.cpp

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

bin_PROGRAMS =
if BUILD_palinear
bin_PROGRAMS += palinear
endif

if BUILD_palogist
bin_PROGRAMS += palogist
endif

if BUILD_pacoxph
## Extra files for pacoxph
COXBASE = coxfit2 chinv2 cholesky2 chsolve2 dmatrix
COXSRC = $(COXBASE:=.c)
bin_PROGRAMS += pacoxph
endif

if BUILD_extractsnp
bin_PROGRAMS += extract-snp
endif


palinear_SOURCES = $(REGFILES) $(FVSRC) $(FVHEADERS)
palinear_CPPFLAGS = -DLINEAR $(AM_CPPFLAGS)
if WITH_EIGEN
palinear_SOURCES += $(EIGENFILES)
palinear_CPPFLAGS += -DEIGEN
endif

palogist_SOURCES = $(REGFILES) $(FVSRC) $(FVHEADERS)
palogist_CPPFLAGS = -DLOGISTIC $(AM_CPPFLAGS)
if WITH_EIGEN
palogist_SOURCES += $(EIGENFILES)
palogist_CPPFLAGS += -DEIGEN
endif

pacoxph_SOURCES = $(COXSRC) $(REGFILES) $(FVSRC) $(FVHEADERS)	\
 $(RHEADERS) survS.h survproto.h coxph_data.h coxph_data.cpp
pacoxph_CXXFLAGS = -I $(top_srcdir)/src/include $(AM_CXXFLAGS)
pacoxph_CPPFLAGS = -DCOXPH $(AM_CPPFLAGS)
pacoxph_CFLAGS = -I $(top_srcdir)/src/include $(AM_CFLAGS)
if WITH_EIGEN
pacoxph_SOURCES += $(EIGENFILES)
pacoxph_CPPFLAGS += -DEIGEN
endif

extract_snp_SOURCES = extract-snp.cpp $(FVSRC) $(FVHEADERS)

## Install these scripts in the bin directory as well:
dist_bin_SCRIPTS = probabel

## Install this R script in the examples directory
scriptdir = $(pkgdatadir)/scripts
dist_script_DATA = prepare_data.R extIDS.pl

## Install the config file
dist_sysconf_DATA = probabel_config.cfg.example

## Insert the correct values of several variables in the installed
## probabel script
install-exec-hook:
	$(SED) -i -e "s;probabel_config.cfg;@sysconfdir@/probabel_config.cfg;" \
		$(DESTDIR)$(bindir)/probabel
	$(SED) -i -e "s;\./;@bindir@/;" $(DESTDIR)$(bindir)/probabel
	$(SED) -i -e "s;PROBABEL_VERSION;@VERSION@;" $(DESTDIR)$(bindir)/probabel
	ln -s $(DESTDIR)$(bindir)/probabel $(DESTDIR)$(bindir)/probabel.pl

## Make the script executable
install-data-hook:
	chmod a+rx $(DESTDIR)/$(scriptdir)/extIDS.pl
