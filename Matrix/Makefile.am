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

libpzmatrix_a_SOURCES = pzblock.cpp \
	pzblockdiag.cpp \
	pzbndmat.cpp \
	pzespmat.cpp \
	pzfmatrix.cpp \
	pzlink.cpp \
	pzmatred.cpp \
	pzmatrix.cpp \
	pzsbndmat.cpp \
	pzseqsolver.cpp \
	pzsespmat.cpp \
	pzsfulmat.cpp \
	pzshtmat.cpp \
	pzskylmat.cpp \
	pzskylmatpar.cpp \
	pzsolve.cpp \
	pzstencil.cpp \
	pzstepsolver.cpp \
	pzsysmp.cpp \
	pztempmat.cpp \
	pzworkpool.cpp \
	pzysmp.cpp 

## --| NeoPZ |------------------------------------------------------------------
