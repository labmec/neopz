#
# NeoPZ adaptivity makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

bin_PROGRAMS = main

main_SOURCES = main.cc 

main_LDADD = $(top_builddir)/lib/libpz.a @METISDIR@/lib/libmetis.a /home/ic/luis/Documents/dxflib-2.0.1.3-1.src/lib/libdxf.a

## --| NeoPZ |------------------------------------------------------------------
