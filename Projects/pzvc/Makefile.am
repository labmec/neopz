#
# NeoPZ pzvc makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

bin_PROGRAMS = main

main_SOURCES = main.cpp pzadaptmesh.cpp pzcclonemesh.cpp pzgclonemesh.cpp

main_LDADD = $(top_srcdir)/lib/libpz.a @METISDIR@/Lib/libmetis.a

## --| NeoPZ |------------------------------------------------------------------
