#
# NeoPZ Util makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@


utildir = $(prefix)/Util

util_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzutil.a

libpzutil_a_SOURCES = pzmap.cpp pzstring.cpp pzline.cpp pzplane.cpp pznumeric.cpp pzpolynomial.cpp TPZTimer.cpp TPZTimer.h

## --| NeoPZ |------------------------------------------------------------------
