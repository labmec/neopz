#
# NeoPZ Frontal makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

frontaldir = $(prefix)/Frontal

frontal_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzfrontal.a

libpzfrontal_a_SOURCES = pzvecfrontal.cc \
	tpzeqnarray.cc \
	TPZFileEqnStorage.cc \
	TPZFront.cc \
	TPZFrontMatrix.cc \
	TPZFrontNonSym.cc \
	TPZFrontSym.cc \
	TPZParFrontMatrix.cc \
	TPZStackEqnStorage.cc 

## --| NeoPZ |------------------------------------------------------------------
