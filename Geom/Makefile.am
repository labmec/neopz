#
# NeoPZ Geom makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

geomdir = $(prefix)/Geom

geom_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzgeom.a

libpzgeom_a_SOURCES = pzgeopoint.cc \
	pzgeoprism.cc \
	pzgeopyramid.cc \
	pzgeoquad.cc \
	pzgeotetrahedra.cc \
	pzgeotriangle.cc \
	TPZGeoCube.cc \
	TPZGeoElement.cc \
	TPZGeoLinear.cc 

## --| NeoPZ |------------------------------------------------------------------
