#
# NeoPZ StrMatrix makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

strmatrixdir = $(prefix)/StrMatrix

strmatrix_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzstrmatrix.a

libpzstrmatrix_a_SOURCES = pzbdstrmatrix.cpp \
	pzbstrmatrix.cpp \
	pzfstrmatrix.cpp \
	pzsbstrmatrix.cpp \
	pzskylstrmatrix.cpp \
	pzstrmatrix.cpp \
	TPBSpStructMatrix.cpp \
	TPZFrontStructMatrix.cpp \
	TPZParFrontStructMatrix.cpp \
	TPZParSkylineStructMatrix.cpp \
	TPZSpStructMatrix.cpp

## --| NeoPZ |------------------------------------------------------------------
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@
INCLUDES = -I$(top_srcdir)/Analysis -I$(top_srcdir)/Common -I$(top_srcdir)/External -I$(top_srcdir)/Frontal -I$(top_srcdir)/Geom -I$(top_srcdir)/Integral -I$(top_srcdir)/Material -I$(top_srcdir)/Matrix -I$(top_srcdir)/Mesh -I$(top_srcdir)/Multigrid -I$(top_srcdir)/Post -I$(top_srcdir)/Refine -I$(top_srcdir)/Save -I$(top_srcdir)/Shape -I$(top_srcdir)/Util @METISINCLUDE@ @FADINCLUDE_TFADET@ @FADINCLUDE_TFAD@ @FADINCLUDE_FAD@ @FADINCLUDE@
