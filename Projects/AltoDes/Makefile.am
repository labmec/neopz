#
# NeoPZ High Performance Makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@ -DCONVERGENCE_OUTPUT

bin_PROGRAMS = main

main_SOURCES = main.cpp 

main_LDADD = $(top_builddir)/lib/libpz.a @METISDIR@/lib/libmetis.a

## --| NeoPZ |------------------------------------------------------------------
