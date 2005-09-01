#
# NeoPZ Cosys makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

cosysdir = $(prefix)/Cosys

cosysincdir = $(includedir)/Cosys
cosysinc_HEADERS= $(top_srcdir)/Cosys/*.h

#cosys_HEADERS = pzcartsys.h \
#	pzcosys.h \
#	pzcylinsys.h \
#	pzelgq2dcyl.h \
#	pzesfersys.h 


noinst_LIBRARIES = libpzcosys.a

libpzcosys_a_SOURCES = pzcartsys.cpp pzcosys.cpp pzcylinsys.cpp pzelgq2dcyl.cpp pzesfersys.cpp 

## --| NeoPZ |------------------------------------------------------------------

