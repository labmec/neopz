DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

bin_PROGRAMS = lara 



## --| NeoPZ |------------------------------------------------------------------

lara_SOURCES = adapt.cpp main.cpp malhas.cpp pzexplfinvolanal.cpp \
	pzexplfinvolanal.h tools.cpp

lara_DEPENDENCIES = $(top_builddir)/lib/libpz.a

# noinst_HEADERS = TPZFakeFunction.h
# _SOURCES = adapt.h
