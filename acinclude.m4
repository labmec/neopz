dnl
dnl acinclude.m4 for NeoPZ
dnl
dnl Process this file with GNU aclocal to produce a configure script.
dnl
dnl $Id: acinclude.m4,v 1.3 2003-06-03 13:52:01 erick Exp $
dnl

dnl
dnl Greetings!
dnl
AC_DEFUN(PZ_GREETINGS,
[
    echo
    echo "+-----------------------------------------------+"
    echo "             Welcome to NeoPZ project"
    echo "+-----------------------------------------------+"
    echo
    echo "Configuring PZ version:" $PZ_VERSION.$PZ_REV
    echo
])

dnl
dnl Checking g++ version
dnl
AC_DEFUN(PZ_PROG_CXX,
[
    AC_PROG_CXX
    case "$CXX" in
        c++ | g++)
           CXX_MAJOR=2
           CXX_MINOR=95
           AC_MSG_CHECKING(if $CXX version >= $CXX_MAJOR.$CXX_MINOR)
           AC_TRY_COMPILE([#include<features.h>],
             [
              #if !__GNUC_PREREQ($CXX_MAJOR, $CXX_MINOR)
              #error Bad version
              #endif
             ],
             AC_MSG_RESULT(ok),
             AC_MSG_ERROR($CXX invalid version! Must be >= $CXX_MAJOR.$CXX_MINOR))
             ;;
    esac
])

dnl
dnl Checking for ar
dnl
AC_DEFUN(PZ_PROG_AR,
[
    case "${AR-unset}" in
	unset) AC_CHECK_PROG(AR, ar, ar) ;;
	*) AC_CHECK_PROGS(AR, $AR ar, ar) ;;
    esac
    AC_SUBST(AR)
    AC_MSG_CHECKING(ar flags)
    case "${ARFLAGS-unset}" in
	unset) ARFLAGS="-rcsv" ;;
    esac
    AC_MSG_RESULT($ARFLAGS)
    AC_SUBST(ARFLAGS)
])

dnl
dnl Bye bye!
dnl
AC_DEFUN(PZ_BYEBYE,
[
    echo
    echo "Finished configuration for PZ version" $PZ_VERSION.$PZ_REV
    echo

    echo "+-----------------------------------------------+"
    echo
    echo "   You hopefully configured NeoPZ project"
    echo
    echo "   Options:"
    echo

    case "${metis_enabled}" in
      yes)
        echo "      -> MeTiS enabled."
      ;;
      no)
        echo "      -> MeTiS not enabled."
      ;;
    esac

    case "${sloan_enabled}" in
      yes)
        echo "      -> Sloan enabled."
      ;;
      no)
        echo "      -> Sloan not enabled."
      ;;
    esac

    case "${fad_enabled}" in
      yes)
        echo "      -> FAD enabled."
      ;;
      no)
        echo "      -> FAD not enabled."
      ;;
    esac

    echo
    echo "   type \"make\" to start compilation."
    echo "   type \"make install\" as root to install it."
    echo
    echo "+-----------------------------------------------+"
])

dnl --| NeoPZ |-----------------------------------------------------------------
