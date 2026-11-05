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

## Do not compile generate.cc into the pz library (it contains a standalone main)
set(sources)

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/External/FAD/TinyFad/Specializations)

# Add header sources to the pz library. Keep generate.cc as a separate small tool.
target_sources(pz PRIVATE ${headers})

# Create a tiny executable from generate.cc so the file is not included in the library
# (generate.cc contains a main() used to produce header specializations).
add_executable(tinyfad_generate generate.cc)
target_include_directories(tinyfad_generate PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Optionally install the tool so users can run it after `make install`.
install(TARGETS tinyfad_generate
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
