#
# NeoPZ Matrix makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

matrixdir = $(prefix)/Matrix

matrixincdir = $(includedir)/Matrix
matrixinc_HEADERS= $(top_srcdir)/Matrix/*.h

noinst_LIBRARIES = libpzmatrix.a

libpzmatrix_a_SOURCES = pzblockdiag.cpp pzbndmat.cpp pzespmat.cpp pzfmatrix.cpp pzlink.cpp pzmatred.cpp pzmatrix.cpp pzsbndmat.cpp pzseqsolver.cpp pzsespmat.cpp pzsfulmat.cpp pzshtmat.cpp pzskylmat.cpp pzskylmatpar.cpp pzsolve.cpp pzstencil.cpp pzblock.cpp pzstepsolver.cpp pzsysmp.cpp pzysmp.cpp pzmatrixid.cpp tpzsparseblockdiagonal.cpp pzspblockdiagpivot.cpp tpzverysparsematrix.cpp

## --| NeoPZ |------------------------------------------------------------------
noinst_HEADERS = pzspblockdiagpivot.h
