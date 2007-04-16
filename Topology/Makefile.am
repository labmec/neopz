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

libpzstrmatrix_a_SOURCES = pztriangle.cpp

## --| NeoPZ |------------------------------------------------------------------
