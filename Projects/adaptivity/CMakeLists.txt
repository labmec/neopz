
set(main_SRCS
    main.cpp
    pzadaptmesh.cpp
    pzcclonemesh.cpp
    pzgclonemesh.cpp)

add_executable(adaptivity ${main_SRCS})

if(USING_METIS)
	target_link_libraries(adaptivity pz metis)
else()
	target_link_libraries(adaptivity pz)
endif()

