#
# NeoPZ Refine makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id# Testing 2
#

-include $(top_builddir)/@PZ_MK@


DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

refinedir = $(prefix)/Refine

refineincdir = $(includedir)/Refine
refineinc_HEADERS= $(top_srcdir)/Refine/*.h $(top_srcdir)/Refine/RefPatterns/*


noinst_LIBRARIES = libpzrefine.a

libpzrefine_a_SOURCES = pzrefpoint.cpp pzrefprism.cpp pzrefpyram.cpp pzrefquad.cpp pzreftetrahedra.cpp pzreftriangle.cpp TPZRefCube.cpp TPZRefLinear.cpp TPZRefPattern.cpp TPZRefPatternDataBase.cpp TPZRefPatternTools.cpp

## --| NeoPZ |------------------------------------------------------------------
