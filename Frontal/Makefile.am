#
# NeoPZ Frontal makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

frontaldir = $(prefix)/Frontal

frontal_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzfrontal.a

libpzfrontal_a_SOURCES = pzvecfrontal.cpp \
	tpzeqnarray.cpp \
	TPZFileEqnStorage.cpp \
	TPZFront.cpp \
	TPZFrontMatrix.cpp \
	TPZFrontNonSym.cpp \
	TPZFrontSym.cpp \
	TPZParFrontMatrix.cpp \
	TPZStackEqnStorage.cpp 

## --| NeoPZ |------------------------------------------------------------------
