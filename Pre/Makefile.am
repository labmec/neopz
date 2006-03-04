#
# NeoPZ Pre makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

predir = $(prefix)/Pre

preincdir = $(includedir)/Pre
preinc_HEADERS= $(top_srcdir)/Pre/*.h


noinst_LIBRARIES = libpzpre.a

libpzpre_a_SOURCES = pzgengrid.cpp pzpargrid.cpp TPZGMSHReadMesh.cpp  TPZExtendGridDimension.cpp pzreadtetgen.cpp

## --| NeoPZ |------------------------------------------------------------------
