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

libpzshape_a_SOURCES = pzshapecube.cpp \
	pzshapelinear.cpp \
	pzshapepiram.cpp \
	pzshapeprism.cpp \
	pzshapequad.cpp \
	pzshapetetra.cpp \
	pzshapetriang.cpp \
	TPZShapeDisc.cpp 

## --| NeoPZ |------------------------------------------------------------------
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@ @ATLASINCLUDE@ 
INCLUDES = -I$(top_srcdir)/Common -I$(top_srcdir)/Matrix -I$(top_srcdir)/Multigrid -I$(top_srcdir)/Save -I$(top_srcdir)/Util @FADINCLUDE@ @FADINCLUDE_TFADET@ @FADINCLUDE_TFAD@ @FADINCLUDE_FAD@
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@
