#
# NeoPZ LinearSolvers makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

linearsolversdir = $(prefix)/LinearSolvers

linearsolversincdir = $(includedir)/LinearSolvers
linearsolversinc_HEADERS= $(top_srcdir)/LinearSolvers/*.h

## --| NeoPZ |------------------------------------------------------------------
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@ @ATLASINCLUDE@ 
