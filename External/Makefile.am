#
# NeoPZ External makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

externaldir = $(prefix)/External

external_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzexternal.a

libpzexternal_a_SOURCES = pzmetis.cpp 	pzrenumbering.cpp 	pzsloan.cpp tpznodesetcompute.cpp

## --| NeoPZ |------------------------------------------------------------------
