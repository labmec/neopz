# @file neopz/External/FAD/TinyFad/Specializations/CMakeLists.txt

target_include_directories(pz PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/External/FAD/TinyFad/Specializations>
        )

set(headers
    tinyfadone.h
    tinyfadtwo.h
    tinyfadthree.h
    tinyfadfour.h
    tinyfadfive.h
    tinyfadsix.h
    tinyfadseven.h
    tinyfadeight.h
    tinyfadnine.h
    tinyfadten.h
    tinyfadeleven.h
    tinyfadtwelve.h
    tinyfadthirteen.h
    tinyfadfourteen.h
    tinyfadfifteen.h
    tinyfadsixteen.h
    tinyfadseventeen.h
    tinyfadeighteen.h
    tinyfadnineteen.h
    tinyfadtwenty.h
    )

set(sources generate.cc)

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/External/FAD/TinyFad/Specializations)

target_sources(pz PRIVATE ${headers} ${sources})
