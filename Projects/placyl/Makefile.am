#
# NeoPZ placyl makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

bin_PROGRAMS = main

main_SOURCES = main.cpp tpzmatplacaiso.cpp pzelgq2dcyl.cpp

main_LDADD = $(top_srcdir)/lib/libpz.a @METISDIR@/Lib/libmetis.a

## --| NeoPZ |------------------------------------------------------------------
