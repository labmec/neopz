#
# NeoPZ Multigrid makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

multigriddir = $(prefix)/Multigrid

multigrid_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzmultigrid.a

libpzmultigrid_a_SOURCES = pzmgsolver.cpp \
	pztransfer.cpp \
	pztrnsform.cpp 

## --| NeoPZ |------------------------------------------------------------------
INCLUDES = -I$(top_srcdir)/Common -I$(top_srcdir)/LinearSolvers -I$(top_srcdir)/Matrix -I$(top_srcdir)/Mesh -I$(top_srcdir)/Save -I$(top_srcdir)/Util
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@ @ATLASINCLUDE@ 
