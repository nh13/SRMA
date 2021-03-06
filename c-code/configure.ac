#  See LICENSE for details
#
#  The version number in AC_INIT is for the srma package.
#  Developers - the version number should *NEVER* be changed except by, 
#  or with explicit permission from, Nils Homer.

# TODO: check that samtools sub-directory exists

AC_PREREQ(2.59);
AC_INIT([srma],[0.1.16],[srma-help@lists.sourceforge.net],[srma])
AC_COPYRIGHT([See LICENSE for copyright information.])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([dist-bzip2 subdir-objects])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])

# Sets build_cpu, build_vendor and build_os variables
AC_CANONICAL_BUILD

# Check for C and C++ compilers
AC_PROG_CC
AM_PROG_CC_C_O

AC_PROG_INSTALL
AC_GNU_SOURCE

# set CFLAGS and CXXFLAGS
default_CFLAGS="-Wall -g -O2 -pthread";
extended_CFLAGS="";# "-m64 -D_FILE_OFFSET_BITS=64";

# Define some variables
GITREV="Revision: undefined$";
gitrev="`git log | grep -m 1 commit | awk '{print $2}'`";
if test -n "${gitrev}"; then
	GITREV="${gitrev}";
fi
AC_REVISION(["${GITREV}"])
AC_DEFINE_UNQUOTED(GIT_REV, "${GITREV}", [Define to be the git revision number])

AC_ARG_ENABLE(intel64, [  --enable-intel64        optimize for Intel64 CPU such as Xeon and Core2], [extended_CFLAGS="${extended_CFLAGS} -mtune=nocona"], [])

CFLAGS="${default_CFLAGS} ${extended_CFLAGS}";

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
# Checks for library functions.
AC_FUNC_MALLOC 
AC_FUNC_REALLOC
AC_CHECK_LIB([m], [pow])
AC_CHECK_FUNC(pow, [AC_DEFINE(HAVE_POW, 1, [Define 1 if you have the function pow])],[])
AC_CHECK_FUNC(strdup, [AC_DEFINE(HAVE_STRDUP, 1, [Define 1 if you have the function strdup])],[])

# Check types
AC_CHECK_TYPE(int8_t)
AC_CHECK_TYPE(int32_t)
AC_CHECK_TYPE(int64_t)
AC_CHECK_TYPE(uint8_t)
AC_CHECK_TYPE(uint32_t)
AC_CHECK_TYPE(uint64_t)

# Set proprocessor SIZEOF_short_int, SIZEOF_int, SIZEOF_long_int
AC_CHECK_SIZEOF([short int])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long int])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h stdio.h stdint.h string.h sys/time.h unistd.h float.h])

AC_C_INLINE

# test for SAMtools
AC_MSG_CHECKING(for samtools)
if test -d "samtools" ; then
    AC_MSG_RESULT(yes)
    LDFLAGS="${LDFLAGS} -L./samtools -lbam";
    CFLAGS="${CFLAGS} -I./samtools";
    AC_DEFINE(HAVE_SAMTOOLS, 1, [Define 1 if the samtools library is present])
else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR("could not find the samtools source sub-directory")
fi

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
