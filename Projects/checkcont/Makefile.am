#
# NeoPZ checkcont makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

bin_PROGRAMS = continue2

continue2_SOURCES = continue2.cpp

continue2_LDADD = $(top_srcdir)/lib/libpz.a @METISDIR@/lib/libmetis.a

## --| NeoPZ |------------------------------------------------------------------
