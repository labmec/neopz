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

bin_PROGRAMS = multigrid

multigrid_SOURCES = main.cpp

multigrid_LDADD = $(PZ_LIBDIR) $(PZ_LIB) -lpthread 

## --| NeoPZ |------------------------------------------------------------------
