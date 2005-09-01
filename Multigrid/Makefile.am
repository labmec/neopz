#
# NeoPZ Multigrid makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

multigriddir = $(prefix)/Multigrid

multigridincdir = $(includedir)/Multigrid
multigridinc_HEADERS= $(top_srcdir)/Multigrid/*.h

noinst_LIBRARIES = libpzmultigrid.a

libpzmultigrid_a_SOURCES = pzmgsolver.cpp pztransfer.cpp pztrnsform.cpp 

## --| NeoPZ |------------------------------------------------------------------
