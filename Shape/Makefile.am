#
# NeoPZ Shape makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

shapedir = $(prefix)/Shape

shape_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzshape.a

libpzshape_a_SOURCES = pzshapecube.cc \
	pzshapelinear.cc \
	pzshapepiram.cc \
	pzshapeprism.cc \
	pzshapequad.cc \
	pzshapetetra.cc \
	pzshapetriang.cc \
	TPZShapeDisc.cc 

## --| NeoPZ |------------------------------------------------------------------
