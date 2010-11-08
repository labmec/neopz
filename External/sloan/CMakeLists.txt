
########### next target ###############

set(pzexternal_STAT_SRCS
    ${PROJECT_SOURCE_DIR}/External/sloan/diamtr.cpp
    ${PROJECT_SOURCE_DIR}/External/sloan/gegra.cpp
    ${PROJECT_SOURCE_DIR}/External/sloan/iorele.cpp
    ${PROJECT_SOURCE_DIR}/External/sloan/isorti.cpp
    ${PROJECT_SOURCE_DIR}/External/sloan/label.cpp
    ${PROJECT_SOURCE_DIR}/External/sloan/number.cpp
    ${PROJECT_SOURCE_DIR}/External/sloan/profi1.cpp
    ${PROJECT_SOURCE_DIR}/External/sloan/rootls.cpp
    ${PROJECT_SOURCE_DIR}/External/sloan/vsrtp1.cpp
	PARENT_SCOPE)


########### install files ###############




#original Makefile.am contents follow:

#INCLUDES = 
#METASOURCES = AUTO
#
#noinst_LIBRARIES = libsloan.a
#
#libsloan_a_SOURCES = diamtr.cpp gegra.cpp iorele.cpp isorti.cpp label.cpp\
#	number.cpp profi1.cpp rootls.cpp sloan.h vsrtp1.cpp
