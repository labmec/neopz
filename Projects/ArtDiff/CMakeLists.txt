
set(main_SRCS main.cpp)

add_executable(ArtDiff ${main_SRCS})

if(USING_METIS)
	target_link_libraries(ArtDiff pz metis)
else()
	target_link_libraries(ArtDiff pz)
endif()


