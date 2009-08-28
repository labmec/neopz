DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

bin_PROGRAMS = main 



## --| NeoPZ |------------------------------------------------------------------

main_SOURCES = main.cpp malhas.cpp tools.cpp

main_DEPENDENCIES = $(top_builddir)/lib/libpz.a

noinst_HEADERS = TPZFakeFunction.h
