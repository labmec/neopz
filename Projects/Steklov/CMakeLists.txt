
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(steklov_SRCS substructeigen.cpp tpzcurvedtriangle.cpp GenerateSquare.cpp)

add_executable(steklov ${steklov_SRCS})

target_link_libraries(steklov ${KDE4_KDECORE_LIBS})

install(TARGETS steklov ${INSTALL_TARGETS_DEFAULT_ARGS})

