#
# NeoPZ adaptivity makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

bin_PROGRAMS = Adapt_MultThread

Adapt_MultThread_SOURCES = main.cc pzadaptmesh.cc pzcclonemesh.cc pzgclonemesh.cc c0-simplequad.cc \
		c01-lshape-quad.cc c02-lshape-triang.cc c03-lshape-quadtri.cc \
		c04-simplehexa.cc c05-cornercube.cc c06-Tetra.cc c07-Prism.cc \
		c08-Mixed.cc c09-Plate-Cedric.cc c10-Prism-Pref.cc c11-Piram-Tetra.cc \
		c15-Exp-Hexa.cc

Adapt_MultThread_LDADD =  -lpthread $(top_srcdir)/lib/libpz.a @METISDIR@/lib/libmetis.a 

## --| NeoPZ |------------------------------------------------------------------
