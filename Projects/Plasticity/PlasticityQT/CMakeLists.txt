FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

IF (UNIX)
  include_directories(/usr/local/include)
  set (QWT_LIB /usr/local/lib/libqwt.so)
ENDIF()

IF (APPLE)
  include_directories(/usr/local/qwt-6.0.3-svn/lib/qwt.framework/Headers)
  set (QWT_LIB /usr/local/qwt-6.0.3-svn/lib/qwt.framework/qwt)
ENDIF()

# UI Files
QT4_WRAP_UI(UISrcs initialpointdock.ui mainwindow.ui globalconfig.ui)
QT4_WRAP_CPP(MOCSrcs initialpointdock.h mainwindow.h globalconfig.h canvaspicker.h plot.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS main.cpp ${sources} ${headers})

ADD_EXECUTABLE(PlasticityQT ${main_SRCS} ${MOCSrcs} ${UISrcs})
TARGET_LINK_LIBRARIES(PlasticityQT ${QT_QTMAIN_LIBRARY} ${QT_LIBRARIES} ${QWT_LIB} pz)

file(COPY ../ensaio_all_columns.txt DESTINATION .)
