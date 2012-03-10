AC_INIT([rcdd])
AC_CONFIG_SRCDIR([DESCRIPTION])

PKG_CPPFLAGS=$CPPFLAGS
PKG_CFLAGS=$CFLAGS
PKG_LDFLAGS=$LDFLAGS

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

AC_CHECK_LIB(gmp, __gmpz_init, [],
    [AC_MSG_ERROR([GNU MP not found, see README])])

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LDFLAGS)

AC_CONFIG_FILES(src/Makevars)
AC_OUTPUT
