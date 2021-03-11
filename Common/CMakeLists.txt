# @file neopz/Common/CMakeLists.txt  -- CMake file for the Common module

target_include_directories(pz PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Common>
                           )

set(headers
    TPZBlas.h
    TPZGuiInterface.h
    TPZLapack.h
    TPZReschedulableTask.h
    TPZTask.h
    TPZTaskGroup.h
    TPZThreadPool.h
    doxcommon.h
    fpo_exceptions.h
    pzeltype.h
    pzerror.h
    pzreal.h
    pzversion.h
    )

set(sources
    TPZGuiInterface.cpp
    TPZReschedulableTask.cpp
    TPZTask.cpp
    TPZTaskGroup.cpp
    TPZThreadPool.cpp
    pzerror.cpp
    pzreal.cpp
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Common)

target_sources(pz PRIVATE ${headers} ${sources})
