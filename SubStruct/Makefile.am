-include /home/longhin/Projetos/PZInstall/PZ_user.mk

bin_PROGRAMS = substruct
substruct_SOURCES = converge substruct.cpp time.txt tpzdohrmatrix.cpp \
	tpzdohrprecond.cpp tpzdohrsubstruct.cpp tpzgensubstruct.cpp

# set the include path found by configure
INCLUDES= $(PZ_INCLUDES) $(all_includes)

# the library search path.
substruct_LDFLAGS = $(PZ_LIBDIR) $(all_libraries) 

AM_CXXFLAGS = $(PZ_CXXFLAGS)

substruct_LDADD = -L/hd/longhin/Projetos/metis-4.0/ $(PZ_LIB)
noinst_HEADERS = tpzdohrsubstruct.h tpzdohrmatrix.h tpzdohrprecond.h \
	tpzgensubstruct.h
