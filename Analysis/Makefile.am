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

analysisincdir = $(includedir)/Analysis
analysisinc_HEADERS= $(top_srcdir)/Analysis/*.h

noinst_LIBRARIES = libpzanalysis.a

libpzanalysis_a_SOURCES = pzanalysis.cpp pzanalysiserror.cpp pzmganalysis.cpp pznonlinanalysis.cpp pzonedref.cpp pzsmanal.cpp pzsmfrontalanal.cpp pzeuleranalysis.cpp pztransientanalysis.cpp

## --| NeoPZ |------------------------------------------------------------------
