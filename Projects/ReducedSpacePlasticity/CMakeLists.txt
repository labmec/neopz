file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

# Plasticity materials work only with double type.
if ((STATE_TYPE STREQUAL "double") AND (REAL_TYPE STREQUAL "double"))

  if(BUILD_PLASTICITY_MATERIALS)

    add_executable(ReducedSpacePlasticity ${main_SRCS} ${headers})
    target_link_libraries(ReducedSpacePlasticity pz)

  endif()

endif()