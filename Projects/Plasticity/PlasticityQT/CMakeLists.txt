FIND_PROGRAM(QT_QMAKE_EXECUTABLE NAMES qmake PATHS "/usr/local/opt/qt5/bin")



FIND_PACKAGE(Qt REQUIRED)
INCLUDE(${QT_USE_FILE})

IF (UNIX)
  include_directories(/usr/local/include)
  set (QWT_LIB /usr/local/lib/libqwt.so)
ENDIF()

IF (APPLE)
  include_directories(/usr/local/Cellar/qwt/6.1.3_3/lib/qwt.framework/Headers)
  set (QWT_LIB /usr/local/Cellar/qwt/6.1.3_3/lib/qwt.framework)
ENDIF()

# UI Files
QT5_WRAP_UI(UISrcs mainwindow.ui globalconfig.ui)
QT5_WRAP_CPP(MOCSrcs mainwindow.h globalconfig.h canvaspicker.h plot.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

file(GLOB headers *.h)
file(GLOB sources *.cpp)

file(GLOB headersTP TP*.h)
file(GLOB sourcesTP TP*.cpp)

file (GLOB to_remove_from_sources main.cpp mainCLI.cpp)
list (REMOVE_ITEM sources ${to_remove_from_sources})

set(main_SRCS main.cpp ${sources} ${headers})
set(main_SRCS_CLI mainCLI.cpp ${sourcesTP} ${headersTP})

ADD_EXECUTABLE(PlasticityQT ${main_SRCS} ${MOCSrcs} ${UISrcs})
TARGET_LINK_LIBRARIES(PlasticityQT ${QT_QTMAIN_LIBRARY} ${QT_LIBRARIES} ${QWT_LIB} pz)

ADD_EXECUTABLE(PlasticityCLI ${main_SRCS_CLI})
TARGET_LINK_LIBRARIES(PlasticityCLI pz)


file(COPY ../ensaio_all_columns.txt DESTINATION .)
