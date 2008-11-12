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



## --| NeoPZ |------------------------------------------------------------------

main_SOURCES = main.cpp malha.cpp

main_DEPENDENCIES = $(top_builddir)/lib/libpz.a

