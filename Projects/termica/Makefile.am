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



main_LDADD = $(top_builddir)/lib/libpz.a \
	$(top_builddir)/External/sloan/libsloan.a @LIBATLAS@ -lpthread @METISDIR@/lib/libmetis.a

## --| NeoPZ |------------------------------------------------------------------

main_SOURCES = termi.cpp fast.cpp

