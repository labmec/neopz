#
# NeoPZ adaptivity makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

bin_PROGRAMS = main

main_SOURCES = main.cc pzadaptmesh.cc pzcclonemesh.cc pzgclonemesh.cc

main_LDADD = $(top_srcdir)/lib/libpz.a @METISDIR@/Lib/libmetis.a

## --| NeoPZ |------------------------------------------------------------------
