#
# NeoPZ Matrix makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

matrixdir = $(prefix)/Matrix

matrix_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzmatrix.a

libpzmatrix_a_SOURCES = pzblockdiag.cpp 	pzbndmat.cpp 	pzespmat.cpp 	pzfmatrix.cpp 	pzlink.cpp 	pzmatred.cpp 	pzmatrix.cpp 	pzsbndmat.cpp 	pzseqsolver.cpp 	pzsespmat.cpp 	pzsfulmat.cpp 	pzshtmat.cpp 	pzskylmat.cpp 	pzskylmatpar.cpp 	pzsolve.cpp 	pzstencil.cpp	pzblock.cpp 	pzstepsolver.cpp 	pzsysmp.cpp 	pztempmat.cpp 	pzysmp.cpp         pzmatrixid.cpp tpzsparseblockdiagonal.cpp

## --| NeoPZ |------------------------------------------------------------------
INCLUDES = -I$(top_srcdir)/Common -I$(top_srcdir)/LinearSolvers -I$(top_srcdir)/Save -I$(top_srcdir)/Util
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@
noinst_HEADERS =  tpzsparseblockdiagonal.h
