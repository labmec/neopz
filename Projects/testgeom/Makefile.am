#
# NeoPZ adaptivity makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

bin_PROGRAMS = geom

geom_SOURCES = geom.cpp pzgbmesh.cpp

geom_LDADD = $(top_builddir)/lib/libpz.a @METISDIR@/lib/libmetis.a

## --| NeoPZ |------------------------------------------------------------------
