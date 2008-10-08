
DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@ 

INCLUDES = 
METASOURCES = AUTO
bin_PROGRAMS = cccc
cccc_SOURCES = convtest.cpp ctriangletests.cpp tpzarc3d.cpp tpzblendnaca.cpp \
	tpzchangeel.cpp tpzcurvedtriangle.cpp tpzellipse.cpp tpzgeoblend.cpp tpzmathtools.cpp \
	tpzquadraticquad.cpp tpzquadratictetra.cpp tpzquadratictrig.cpp

cccc_DEPENDENCIES = $(top_builddir)/lib/libpz.a
