#
# NeoPZ Post makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

postdir = $(prefix)/Post

post_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzpost.a

libpzpost_a_SOURCES = pzdxmesh.cc \
	pzgraphel1d.cc \
	pzgraphel1dd.cc \
	pzgraphel.cc \
	pzgraphelq2d.cc \
	pzgraphelq2dd.cc \
	pzgraphelq3dd.cc \
	pzgraphmesh.cc \
	pzgraphnode.cc \
	pzmvmesh.cc \
	pztrigraph.cc \
	pztrigraphd.cc \
	pzv3dmesh.cc 

## --| NeoPZ |------------------------------------------------------------------
