#
# NeoPZ Util makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

CXXFLAGS = -I$(top_srcdir)/Common

utildir = $(prefix)/Util

util_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzutil.a

libpzutil_a_SOURCES = pzadmchunkgen.cpp pzavlmapgen.cpp pzchunkgen.cpp pzmanvectorgen.cpp pzmap.cpp pzstackgen.cpp pzstring.cpp pzvecgen.cpp pzline.cpp pzplane.cpp pznumeric.cpp pzpolynomial.cpp

## --| NeoPZ |------------------------------------------------------------------
