#
# NeoPZ Pre makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

predir = $(prefix)/Pre

pre_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzpre.a

libpzpre_a_SOURCES = pzgengrid.cc pzpargrid.cc TPZGMSHReadMesh.cc 

## --| NeoPZ |------------------------------------------------------------------
