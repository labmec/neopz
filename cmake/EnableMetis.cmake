function(enable_metis target)
  # NeoPZ itself calls METIS directly (see External/pzmetis.cpp), so it always
  # needs a real METIS include dir + library, regardless of MUMPS.
  # When MUMPS was built with its own METIS (MUMPS_metis), reuse that exact
  # metis::metis target instead of doing a fresh find_package(METIS): searching
  # again can pick up a different METIS installed on the system (e.g. via
  # Homebrew) and mixing two METIS builds in the same binary leads to
  # ABI/version conflicts.
  if(USING_MUMPS AND MUMPS_metis AND TARGET metis::metis)
    message(STATUS "[NeoPZ] Reusing METIS already configured by MUMPS (metis::metis target)")
    target_link_libraries(${target} PRIVATE metis::metis)
  else()
    find_package(METIS REQUIRED)
    target_link_libraries(${target} PRIVATE ${METIS_LIBRARIES})
    target_include_directories(${target} PRIVATE ${METIS_INCLUDE_DIRS})
  endif()
  target_compile_definitions(${target} PRIVATE PZ_USING_METIS)
  target_compile_definitions(${target} INTERFACE PZ_USING_METIS)
endfunction()
