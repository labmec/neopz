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

libpzgeom_a_SOURCES = pzgeopoint.cpp \
	pzgeoprism.cpp \
	pzgeopyramid.cpp \
	pzgeoquad.cpp \
	pzgeotetrahedra.cpp \
	pzgeotriangle.cpp \
	TPZGeoCube.cpp \
	TPZGeoElement.cpp \
	TPZGeoLinear.cpp 

## --| NeoPZ |------------------------------------------------------------------
INCLUDES = -I$(top_srcdir)/Analysis -I$(top_srcdir)/Common -I$(top_srcdir)/Integral -I$(top_srcdir)/Material -I$(top_srcdir)/Matrix -I$(top_srcdir)/Mesh -I$(top_srcdir)/Multigrid -I$(top_srcdir)/Refine -I$(top_srcdir)/Save -I$(top_srcdir)/Shape -I$(top_srcdir)/Util @FADINCLUDE_TFADET@ @FADINCLUDE_TFAD@ @FADINCLUDE_FAD@ @FADINCLUDE@
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@ @ATLASINCLUDE@ 
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@
