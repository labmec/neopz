#
# NeoPZ Refine makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

refinedir = $(prefix)/Refine

refine_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzrefine.a

libpzrefine_a_SOURCES = pzrefpoint.cc \
	pzrefprism.cc \
	pzrefpyram.cc \
	pzrefquad.cc \
	pzreftetrahedra.cc \
	pzreftriangle.cc \
	TPZRefCube.cc \
	TPZRefLinear.cc 

## --| NeoPZ |------------------------------------------------------------------
