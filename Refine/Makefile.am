#
# NeoPZ Refine makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id# Testing 2
#

-include $(top_srcdir)/@PZ_MK@


DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

refinedir = $(prefix)/Refine

refine_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzrefine.a

libpzrefine_a_SOURCES = pzrefpoint.cpp \
	pzrefprism.cpp \
	pzrefpyram.cpp \
	pzrefquad.cpp \
	pzreftetrahedra.cpp \
	pzreftriangle.cpp \
	TPZRefCube.cpp \
	TPZRefLinear.cpp \
	TPZRefPattern.cpp

## --| NeoPZ |------------------------------------------------------------------
INCLUDES = -I$(top_srcdir)/Analysis -I$(top_srcdir)/Common -I$(top_srcdir)/Geom -I$(top_srcdir)/Integral -I$(top_srcdir)/Material -I$(top_srcdir)/Matrix -I$(top_srcdir)/Mesh -I$(top_srcdir)/Multigrid -I$(top_srcdir)/Save -I$(top_srcdir)/Shape -I$(top_srcdir)/Util @FADINCLUDE_TFADET@ @FADINCLUDE_TFAD@ @FADINCLUDE_FAD@ @FADINCLUDE@
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@ @ATLASINCLUDE@ 
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@
