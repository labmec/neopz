#
# NeoPZ Multigrid makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

multigriddir = $(prefix)/Multigrid

multigrid_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzmultigrid.a

libpzmultigrid_a_SOURCES = pzmgsolver.cc \
	pztransfer.cc \
	pztrnsform.cc 

## --| NeoPZ |------------------------------------------------------------------
