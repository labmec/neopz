#
# NeoPZ Cosys makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

cosysdir = $(prefix)/Cosys

cosys_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzcosys.a

libpzcosys_a_SOURCES = pzcartsys.cpp \
	pzcosys.cpp \
	pzcylinsys.cpp \
	pzelgq2dcyl.cpp \
	pzesfersys.cpp 

## --| NeoPZ |------------------------------------------------------------------
INCLUDES = -I$(top_srcdir)/Common -I$(top_srcdir)/Geom -I$(top_srcdir)/Matrix -I$(top_srcdir)/Mesh -I$(top_srcdir)/Multigrid -I$(top_srcdir)/Save -I$(top_srcdir)/Shape -I$(top_srcdir)/Util @FADINCLUDE_TFADET@ @FADINCLUDE_TFAD@ @FADINCLUDE_FAD@ @FADINCLUDE@
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@ @ATLASINCLUDE@ 
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@
