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

bin_PROGRAMS = main

main_SOURCES = main.cpp pzadaptmesh.cpp pzcclonemesh.cpp pzgclonemesh.cpp

main_LDADD = $(top_builddir)/lib/libpz.a @METISDIR@/Lib/libmetis.a

## --| NeoPZ |------------------------------------------------------------------
