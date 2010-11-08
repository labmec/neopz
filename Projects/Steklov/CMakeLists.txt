
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(steklov_SRCS substructeigen.cpp tpzcurvedtriangle.cpp GenerateSquare.cpp)

kde4_add_executable(steklov ${steklov_SRCS})

target_link_libraries(steklov ${KDE4_KDECORE_LIBS})

install(TARGETS steklov ${INSTALL_TARGETS_DEFAULT_ARGS})


########### install files ###############




#original Makefile.am contents follow:

#
#DISTCLEANFILES = @PZ_DISTCLEANFILES@
#CLEANFILES     = @PZ_CLEANFILES@
#
#-include $(top_builddir)/@PZ_MK@ 
#
#INCLUDES = 
#METASOURCES = AUTO
#bin_PROGRAMS = steklov
#steklov_SOURCES = substructeigen.cpp tpzcurvedtriangle.cpp GenerateSquare.cpp
