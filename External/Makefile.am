#
# NeoPZ External makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

externaldir = $(prefix)/External

external_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzexternal.a

libpzexternal_a_SOURCES = pzmetis.cpp 	pzrenumbering.cpp 	pzsloan.cpp tpznodesetcompute.cpp

## --| NeoPZ |------------------------------------------------------------------
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@ @ATLASINCLUDE@ 
INCLUDES = -I$(top_srcdir)/Common -I$(top_srcdir)/Matrix -I$(top_srcdir)/Mesh -I$(top_srcdir)/Multigrid -I$(top_srcdir)/Save -I$(top_srcdir)/Util @METISINCLUDE@
noinst_HEADERS =  tpznodesetcompute.h
