#
# NeoPZ Post makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

postdir = $(prefix)/Post

post_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzpost.a

libpzpost_a_SOURCES = pzdxmesh.cpp \
	pzgraphel1d.cpp \
	pzgraphel1dd.cpp \
	pzgraphel.cpp \
	pzgraphelq2d.cpp \
	pzgraphelq2dd.cpp \
	pzgraphelq3dd.cpp \
	pzgraphmesh.cpp \
	pzgraphnode.cpp \
	pzmvmesh.cpp \
	pztrigraph.cpp \
	pztrigraphd.cpp \
	pzv3dmesh.cpp \
	pzvisualmatrix.cpp

## --| NeoPZ |------------------------------------------------------------------
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@
INCLUDES = -I$(top_srcdir)/Common -I$(top_srcdir)/Geom -I$(top_srcdir)/Material -I$(top_srcdir)/Matrix -I$(top_srcdir)/Mesh -I$(top_srcdir)/Multigrid -I$(top_srcdir)/Refine -I$(top_srcdir)/Save -I$(top_srcdir)/Shape -I$(top_srcdir)/Util @FADINCLUDE_TFADET@ @FADINCLUDE_TFAD@ @FADINCLUDE_FAD@ @FADINCLUDE@
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@
