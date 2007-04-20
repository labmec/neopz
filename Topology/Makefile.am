#
# NeoPZ [ directory ] makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@


topologydir = $(prefix)/Topology

topologyincdir = $(includedir)/Topology
topologyinc_HEADERS= $(top_srcdir)/Topology/*.h

noinst_LIBRARIES = libpztopology.a


## --| NeoPZ |------------------------------------------------------------------
noinst_HEADERS = tpztriangle.h tpztetrahedron.h tpzquadrilateral.h tpzprism.h \
		tpzpoint.h tpzpyramid.h tpzline.h tpzdisc.h tpzcube.h pznoderep.h.h
libpztopology_a_SOURCES = tpztriangle.cpp tpztetrahedron.cpp \
		tpzquadrilateral.cpp tpzprism.cpp tpzpoint.cpp tpzpyramid.cpp tpzline.cpp tpzdisc.cpp \
	tpzcube.cpp
