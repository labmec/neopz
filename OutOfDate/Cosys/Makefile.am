#
# NeoPZ Cosys makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

cosysdir = $(prefix)/Cosys

cosys_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzcosys.a

libpzcosys_a_SOURCES = pzcartsys.cc \
	pzcosys.cc \
	pzcylinsys.cc \
	pzelgq2dcyl.cc \
	pzesfersys.cc 

## --| NeoPZ |------------------------------------------------------------------
