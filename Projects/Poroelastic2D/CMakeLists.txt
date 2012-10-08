#file(GLOB headers *.h)
#file(GLOB sources *.cpp)

#file (GLOB to_remove_from_sources Validation1D.cpp)
#list (REMOVE_ITEM sources ${to_remove_from_sources})

#set(main_SRCS ${sources} ${headers})

#add_executable(Poroelastic2D ${main_SRCS} ${headers})
#target_link_libraries(Poroelastic2D pz)

add_subdirectory(Validations)
add_subdirectory(PetroleumApplications)
