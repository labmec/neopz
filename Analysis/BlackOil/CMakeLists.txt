# @file neopz/Analysis/BlackOil/CMakeLists.txt

set(public_headers
    pzblackoilanalysis.h
    )

set(headers
    pzblackoilanalysis.h
   )
set(sources
    pzblackoilanalysis.cpp
   )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Analysis/BlackOil)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)