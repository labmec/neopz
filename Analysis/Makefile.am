#
# NeoPZ Analysis makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

analysisdir = $(prefix)/Analysis

analysis_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzanalysis.a

libpzanalysis_a_SOURCES = pzanalysis.cc \
		pzanalysiserror.cc \
		pzmganalysis.cc \
		pznonlinanalysis.cc \
		pzonedref.cc \
		pzsmanal.cc \
		pzsmfrontalanal.cc \
		TPZIterativeAnalysis.cc

## --| NeoPZ |------------------------------------------------------------------
