#
# NeoPZ materialhyperelastico makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

bin_PROGRAMS = testhyperelastic hyperelastic

testhyperelastic_SOURCES = testhyperelastic.cpp
hyperelastic_SOURCES = hyperelastic.cpp

testhyperelastic_LDADD = $(top_builddir)/lib/libpz.a @METISDIR@/lib/libmetis.a
hyperelastic_LDADD = $(top_builddir)/lib/libpz.a @METISDIR@/lib/libmetis.a

## --| NeoPZ |------------------------------------------------------------------
