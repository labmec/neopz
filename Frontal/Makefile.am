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
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@ @ATLASINCLUDE@ 
INCLUDES = -I$(top_srcdir)/Common -I$(top_srcdir)/Matrix -I$(top_srcdir)/Save -I$(top_srcdir)/Util
