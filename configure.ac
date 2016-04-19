# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT(ProbABEL, 0.5.0-beta2, genabel-devel@r-forge.wu-wien.ac.at)
AM_INIT_AUTOMAKE([silent-rules subdir-objects])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([src/mlinfo.h])
AC_CONFIG_HEADERS([src/config.h])

# Add the --disable-maintainer-mode option used by e.g. Debian and
# Gentoo to disable rebuilding the autotools-generated files after
# e.g. a patch has been applied to them. The default is 'disabled', we
# change it here to keep the old behaviour by default.
AM_MAINTAINER_MODE([enable])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_INSTALL
AC_PROG_LN_S


# Add an option to enable some default compiler options used during
# development
AC_ARG_ENABLE([dev-options],
    [AS_HELP_STRING([--enable-dev-options], [enable compile options
    for developers])],
    [devoptions=yes],
    [devoptions=no])

if test "x$devoptions" = "xyes"; then
   AC_MSG_NOTICE([Building with developer options enabled])
   if test -z "$CXXFLAGS"; then
      # User did not set CXXFLAGS, so we can put in our own defaults
      CXXFLAGS=""
   fi
   if test -z "$CPPFLAGS"; then
      # User did not set CPPFLAGS, so we can put in our own defaults
      CPPFLAGS="-Wall"
   fi
   CXXFLAGS+=" -g"
fi

if test "x$devoptions" = "xno"; then
   if test -z "$CXXFLAGS"; then
      # User did not set CXXFLAGS, so we can put in our own defaults
      CXXFLAGS="-O2"
   fi
   if test -z "$CPPFLAGS"; then
      # User did not set CPPFLAGS, so we can put in our own defaults
      CPPFLAGS+="-Wall -DNDEBUG"
   fi
fi

# If CXXFLAGS/CPPFLAGS are already set AC_PROG_CXX will not overwrite them
# with its own defaults
AC_PROG_CXX

# Tell the compiler to build not R version of filevector
AC_SUBST(CPPFLAGS, "$CPPFLAGS -D_NOT_R_FILEVECTOR")


# Since most of our code is in C++, set that language as the default
# for the subsequent checks
AC_LANG_PUSH([C++])

# Checks for libraries.

# Checks for header files.
AC_FUNC_ALLOCA
AC_CHECK_HEADERS([float.h inttypes.h libintl.h limits.h stddef.h \
                          stdint.h stdlib.h string.h sys/param.h \
                          wchar.h wctype.h])

# Option that lets the user specify a path to the EIGEN include files.
AC_ARG_WITH([eigen-include-path],
    [AS_HELP_STRING([--with-eigen-include-path],
      [location of the Eigen headers, defaults to /usr/include/eigen3])],
    [CPPFLAGS+=" -I${withval}"],
    [CPPFLAGS+=' -I/usr/include/eigen3'])

# Check for the EIGEN header files
AC_CHECK_HEADERS([Eigen/Dense Eigen/LU])

if test x$ac_cv_header_Eigen_Dense = xno; then
  AC_MSG_ERROR([Could not find the Eigen header files. Did you specify \
--with-eigen-include-path correctly?])
fi


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
# Disabled the next line to make cross compiling for Windows work.
#AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow putenv sqrt strdup strncasecmp floor])
AC_FUNC_STRTOD

# Check if we can use large (>2GB) files on 32 bit platforms
AC_SYS_LARGEFILE
if test $ac_cv_sys_file_offset_bits != 'no'; then
    CPPFLAGS="$CPPFLAGS -D_FILE_OFFSET_BITS=$ac_cv_sys_file_offset_bits"
fi
AC_FUNC_FSEEKO


# Check for presence of sed. Needed to replace path to /etc in probabel.pl
AC_PROG_SED

# Check for the presence of awk. Needed in the test suite
AC_PROG_AWK

# Check for presence of R
AC_CHECK_PROGS(R, R)
if test -z "$R"; then
  AC_MSG_NOTICE([R not found - Unable to do tests comparing ProbABEL output against R output])
fi
AM_CONDITIONAL([HAVE_R], test -n "$R")

# The user can disable building of the PDF of the manual, for example
# if the required LaTeX packages are not installed
AC_ARG_ENABLE([latex-doc],
    [AS_HELP_STRING([--disable-latex-doc], [disable building the PDF
    documentation from LaTeX source])],
    [latexdoc=no],
    [latexdoc=yes])

if test "x$latexdoc" = "xyes"; then
   AC_MSG_NOTICE([building of the PDF of the user manual from LaTeX source is enabled])
   # Check for presence of pdfLaTeX
   AC_CHECK_PROGS(PDFLATEX, pdflatex)
   if test -z "$PDFLATEX"; then
      AC_MSG_NOTICE([pdflatex not found - Unable to create PDF version of the user manual])
   fi
fi
AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")
AM_CONDITIONAL([BUILD_latexdoc], test "x$latexdoc" = "xyes")


# Check for the presence of Doxygen
AC_CHECK_PROGS([DOXYGEN], doxygen)
if test -z "$DOXYGEN"; then
   AC_MSG_NOTICE([Doxygen not found - Not generating code documentation])
fi
AM_CONDITIONAL([HAVE_DOXYGEN],
               [test -n "$DOXYGEN"]) AM_COND_IF([HAVE_DOXYGEN],
               [AC_CONFIG_FILES([doc/Doxyfile])])


# Building of each of the three modules can be disabled on request
AC_ARG_ENABLE([palinear],
    [AS_HELP_STRING([--disable-palinear], [disable building the palinear program])],
    [palinear=no],
    [palinear=yes])

if test "x$palinear" = "xyes"; then
   AC_MSG_NOTICE([building of palinear is enabled])
fi
AM_CONDITIONAL([BUILD_palinear], test "x$palinear" = "xyes")

AC_ARG_ENABLE([palogist],
    [AS_HELP_STRING([--disable-palogist], [disable building the palogist program])],
    [palogist=no],
    [palogist=yes])

if test "x$palogist" = "xyes"; then
   AC_MSG_NOTICE([building of palogist is enabled])
fi

AM_CONDITIONAL([BUILD_palogist], test "x$palogist" = "xyes")

AC_ARG_ENABLE([pacoxph],
    [AS_HELP_STRING([--disable-pacoxph], [disable building the pacoxph program])],
    [pacoxph=no],
    [pacoxph=yes])

if test "x$pacoxph" = "xyes"; then
   AC_MSG_NOTICE([building of pacoxph is enabled])
fi

AM_CONDITIONAL([BUILD_pacoxph], test "x$pacoxph" = "xyes")


AC_ARG_ENABLE([extract-snp],
    [AS_HELP_STRING([--enable-extract-snp], [enable building the
    extract-snp program. NOTE: this program is not finished yet!])],
    [extractsnp=yes],
    [extractsnp=no])

if test "x$extractsnp" = "xyes"; then
   AC_MSG_NOTICE([building of extract-snp is enabled])
fi

AM_CONDITIONAL([BUILD_extractsnp], test "x$extractsnp" = "xyes")


# Files to be generated by autotools
AC_CONFIG_FILES([
        Makefile
        src/Makefile
        doc/Makefile
        examples/Makefile
        checks/Makefile
        checks/R-tests/Makefile
])

# Create output files
AC_OUTPUT
