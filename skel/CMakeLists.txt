
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(]_STAT_SRCS ${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp)

add_library(] STATIC ${]_STAT_SRCS})


