
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(pzcosys_STAT_SRCS
    pzcartsys.cpp
    pzcosys.cpp
    pzcylinsys.cpp
    pzesfersys.cpp)

add_library(pzcosys STATIC ${pzcosys_STAT_SRCS})


