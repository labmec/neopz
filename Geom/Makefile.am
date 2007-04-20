#
# NeoPZ Geom makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

geomdir = $(prefix)/Geom

geomincdir = $(includedir)/Geom
geominc_HEADERS= $(top_srcdir)/Geom/*.h


noinst_LIBRARIES = libpzgeom.a

libpzgeom_a_SOURCES = pzgeotriangle.cpp pzgeopoint.cpp pzgeoprism.cpp pzgeopyramid.cpp \
	pzgeoquad.cpp pzgeotetrahedra.cpp 	TPZGeoCube.cpp TPZGeoElement.cpp \
	TPZGeoLinear.cpp pznoderep.cpp

## --| NeoPZ |------------------------------------------------------------------

noinst_HEADERS = TPZGeoElement.h.h pznoderep.h
