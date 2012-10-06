
set(main_SRCS multipleload.cpp)

add_executable(multipleload ${main_SRCS})

if(USING_METIS)
	target_link_libraries(multipleload pz pthread ${METIS_LIB})
else()
	target_link_libraries(multipleload pz pthread)
endif()


