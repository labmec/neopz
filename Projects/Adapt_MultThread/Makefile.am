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

bin_PROGRAMS = Adapt_MultThread

Adapt_MultThread_SOURCES = main.cpp pzadaptmesh.cpp pzcclonemesh.cpp pzgclonemesh.cpp c0-simplequad.cpp \
		c01-lshape-quad.cpp c02-lshape-triang.cpp c03-lshape-quadtri.cpp \
		c04-simplehexa.cpp c05-cornercube.cpp c06-Tetra.cpp c07-Prism.cpp \
		c08-Mixed.cpp c09-Plate-Cedric.cpp c10-Prism-Pref.cpp c11-Piram-Tetra.cpp \
		c15-Exp-Hexa.cpp

Adapt_MultThread_LDADD = $(top_builddir)/lib/libpz.a $(top_builddir)/External/sloan/libsloan.a -lpthread @METISDIR@/lib/libmetis.a

## --| NeoPZ |------------------------------------------------------------------
