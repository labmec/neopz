#
# NeoPZ materialhyperelastico makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

bin_PROGRAMS = testhyperelastic hyperelastic

testhyperelastic_SOURCES = testhyperelastic.cc
hyperelastic_SOURCES = hyperelastic.cc

testhyperelastic_LDADD = $(top_srcdir)/lib/libpz.a @METISDIR@/lib/libmetis.a
hyperelastic_LDADD = $(top_srcdir)/lib/libpz.a @METISDIR@/lib/libmetis.a

## --| NeoPZ |------------------------------------------------------------------
