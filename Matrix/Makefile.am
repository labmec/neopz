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

libpzmatrix_a_SOURCES = pzblock.cc \
	pzblockdiag.cc \
	pzbndmat.cc \
	pzespmat.cc \
	pzfmatrix.cc \
	pzlink.cc \
	pzmatred.cc \
	pzmatrix.cc \
	pzsbndmat.cc \
	pzseqsolver.cc \
	pzsespmat.cc \
	pzsfulmat.cc \
	pzshtmat.cc \
	pzskylmat.cc \
	pzskylmatpar.cc \
	pzsolve.cc \
	pzstencil.cc \
	pzstepsolver.cc \
	pzsysmp.cc \
	pztempmat.cc \
	pzworkpool.cc \
	pzysmp.cc 

## --| NeoPZ |------------------------------------------------------------------
