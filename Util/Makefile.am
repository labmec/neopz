#
# NeoPZ Util makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

utildir = $(prefix)/Util

util_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzutil.a

libpzutil_a_SOURCES = pzadmchunkgen.cc		\
	pzavlmap.cc				\
	pzavlmapgen.cc				\
	pzchunkgen.cc				\
	pzmanvectorgen.cc			\
	pzmap.cc				\
	pzstackgen.cc				\
	pzstring.cc				\
	pzvecgen.cc				\
	pzvecmesh.cc				\
	pzline.cc				\
	pzplane.cc				\
	pznumeric.cc				\
	pzpolynomial.cc

## --| NeoPZ |------------------------------------------------------------------
