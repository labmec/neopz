#
# NeoPZ Analysis makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

analysisdir = $(prefix)/Analysis

analysis_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzanalysis.a

libpzanalysis_a_SOURCES = pzanalysis.cpp \
		pzanalysiserror.cpp \
		pzmganalysis.cpp \
		pznonlinanalysis.cpp \
		pzonedref.cpp \
		pzsmanal.cpp \
		pzsmfrontalanal.cpp \
		pzeuleranalysis.cpp 

## --| NeoPZ |------------------------------------------------------------------
