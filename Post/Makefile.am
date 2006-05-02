#
# NeoPZ Post makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

postdir = $(prefix)/Post

postincdir = $(includedir)/Post
postinc_HEADERS= $(top_srcdir)/Post/*.h

noinst_LIBRARIES = libpzpost.a

libpzpost_a_SOURCES = pzdxmesh.cpp pzgraphel1d.cpp pzgraphel1dd.cpp \
	pzgraphel.cpp pzgraphelq2d.cpp pzgraphelq2dd.cpp pzgraphelq3dd.cpp pzgraphmesh.cpp \
	pzgraphnode.cpp pzmvmesh.cpp pztrigraph.cpp pztrigraphd.cpp pzv3dmesh.cpp \
	pzvisualmatrix.cpp tpzgraphelt3d.cpp

## --| NeoPZ |------------------------------------------------------------------

