#
# NeoPZ Integral makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

integraldir = $(prefix)/Integral

integral_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzintegral.a

libpzintegral_a_SOURCES = pzquad.cpp

## --| NeoPZ |------------------------------------------------------------------
