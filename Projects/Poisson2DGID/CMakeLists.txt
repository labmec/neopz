########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(Poisson2DGID ${sources} ${header})

target_link_libraries(Poisson2DGID pz pthread)
	