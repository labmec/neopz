#
# NeoPZ Shape makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

shapedir = $(prefix)/Shape

shape_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzshape.a

libpzshape_a_SOURCES = pzshapecube.cpp \
	pzshapelinear.cpp \
	pzshapepiram.cpp \
	pzshapeprism.cpp \
	pzshapequad.cpp \
	pzshapetetra.cpp \
	pzshapetriang.cpp \
	TPZShapeDisc.cpp 

## --| NeoPZ |------------------------------------------------------------------
