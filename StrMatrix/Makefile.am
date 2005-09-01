#
# NeoPZ StrMatrix makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

strmatrixdir = $(prefix)/StrMatrix

strmatrixincdir = $(includedir)/StrMatrix
strmatrixinc_HEADERS= $(top_srcdir)/StrMatrix/*.h

noinst_LIBRARIES = libpzstrmatrix.a

libpzstrmatrix_a_SOURCES = pzbdstrmatrix.cpp pzbstrmatrix.cpp pzfstrmatrix.cpp pzsbstrmatrix.cpp pzskylstrmatrix.cpp pzstrmatrix.cpp TPBSpStructMatrix.cpp TPZFrontStructMatrix.cpp TPZParFrontStructMatrix.cpp TPZParSkylineStructMatrix.cpp TPZSpStructMatrix.cpp tpzsparseblockdiagonalstructmatrix.cpp

## --| NeoPZ |------------------------------------------------------------------
