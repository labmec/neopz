#
# NeoPZ [ directory ] makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

noinst_LIBRARIES = libSave.a

libSave_a_SOURCES = \
	pzsave.cpp \
	pzfilebuffer.cpp \
	pzbfilestream.cpp

## --| NeoPZ |------------------------------------------------------------------
