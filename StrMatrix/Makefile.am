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

libpzstrmatrix_a_SOURCES = pzbdstrmatrix.cc \
	pzbstrmatrix.cc \
	pzfstrmatrix.cc \
	pzsbstrmatrix.cc \
	pzskylstrmatrix.cc \
	pzstrmatrix.cc \
	TPBSpStructMatrix.cc \
	TPZFrontStructMatrix.cc \
	TPZParFrontStructMatrix.cc \
	TPZParSkylineStructMatrix.cc \
	TPZSpStructMatrix.cc

## --| NeoPZ |------------------------------------------------------------------
