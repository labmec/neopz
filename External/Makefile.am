#
# NeoPZ External makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

externaldir = $(prefix)/External

externalincdir = $(includedir)/External
externalinc_HEADERS= $(top_srcdir)/External/*.h


noinst_LIBRARIES = libpzexternal.a

libpzexternal_a_SOURCES = pzmetis.cpp pzrenumbering.cpp pzsloan.cpp tpznodesetcompute.cpp

## --| NeoPZ |------------------------------------------------------------------
SUBDIRS = sloan
