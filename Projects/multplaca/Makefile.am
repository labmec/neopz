#
# NeoPZ multplaca makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

bin_PROGRAMS = jornadas

jornadas_SOURCES = jornadas.cpp

jornadas_LDADD = $(top_srcdir)/lib/libpz.a @METISDIR@/Lib/libmetis.a

## --| NeoPZ |------------------------------------------------------------------
