#
# NeoPZ Pre makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

predir = $(prefix)/Pre

pre_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzpre.a

libpzpre_a_SOURCES = pzgengrid.cpp pzpargrid.cpp TPZGMSHReadMesh.cpp  TPZExtendGridDimension.cpp

## --| NeoPZ |------------------------------------------------------------------
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@ @ATLASINCLUDE@ 
INCLUDES = -I$(top_srcdir)/Common -I$(top_srcdir)/Geom -I$(top_srcdir)/Integral -I$(top_srcdir)/Matrix -I$(top_srcdir)/Mesh -I$(top_srcdir)/Multigrid -I$(top_srcdir)/Refine -I$(top_srcdir)/Save -I$(top_srcdir)/Shape -I$(top_srcdir)/Util @FADINCLUDE_TFADET@ @FADINCLUDE_TFAD@ @FADINCLUDE_FAD@ @FADINCLUDE@
