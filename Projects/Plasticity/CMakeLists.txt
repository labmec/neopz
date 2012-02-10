file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

if(USING_FAD)

	add_executable(Plasticity ${main_SRCS})

	target_link_libraries(Plasticity pz pthread)

endif()

