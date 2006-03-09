#
# NeoPZ Heman makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#
DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@
-include $(top_builddir)/@PZ_MK@
bin_PROGRAMS =  HeMan
HeMan_SOURCES = c0-simplequad.cpp c1-simplecube.cpp main.cc MeshReader.cpp \
	Quad_1_Rib_Marked.cpp Quad_2_Rib_Ad_Marked.cpp

## tpzdxfdraw.cpp

## noinst_HEADERS =  tpzdxfdraw.h
## HeMan_LDADD = $(top_builddir)/lib/libpz.a @METISLIBDIR@ @ATLASLIBDIR@ @LOG4CXXLIBDIR@ -llog4cxx -lpthread
## /home/ic/luis/Documents/dxflib-2.0.1.3-1.src/lib/libdxf.a

## --| NeoPZ |------------------------------------------------------------------
HeMan_LDFLAGS = -all-static
