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

    # Another dependency can add a directory holding a different metis.h
    # (e.g. /opt/homebrew/include). BEFORE places the MUMPS METIS headers
    # ahead of every other include directory of the target.
    get_target_property(_metis_includes metis::metis INTERFACE_INCLUDE_DIRECTORIES)
    if(_metis_includes)
      target_include_directories(${target} BEFORE PRIVATE ${_metis_includes})
      message(STATUS "[NeoPZ] METIS headers taken first from: ${_metis_includes}")
    endif()

    # libmetis is loaded at run time through the rpath list. A directory such
    # as /opt/homebrew/lib may hold another libmetis with the same file name.
    # Entries set explicitly in BUILD_RPATH/INSTALL_RPATH come before the
    # entries CMake derives from the link directories.
    set(_metis_location "")
    get_target_property(_metis_configs metis::metis IMPORTED_CONFIGURATIONS)
    if(_metis_configs)
      foreach(_cfg IN LISTS _metis_configs)
        get_target_property(_metis_location metis::metis IMPORTED_LOCATION_${_cfg})
        if(_metis_location)
          break()
        endif()
      endforeach()
    endif()
    if(NOT _metis_location)
      get_target_property(_metis_location metis::metis IMPORTED_LOCATION)
    endif()
    if(_metis_location)
      get_filename_component(_metis_libdir "${_metis_location}" DIRECTORY)
      get_target_property(_build_rpath ${target} BUILD_RPATH)
      if(NOT _build_rpath)
        set(_build_rpath "")
      endif()
      set_property(TARGET ${target} PROPERTY BUILD_RPATH "${_metis_libdir}" ${_build_rpath})
      set_property(TARGET ${target} APPEND PROPERTY INSTALL_RPATH "${_metis_libdir}")
      message(STATUS "[NeoPZ] METIS library directory added first to rpath: ${_metis_libdir}")
    else()
      message(WARNING "[NeoPZ] Could not read the location of metis::metis; rpath left unchanged")
    endif()
  else()
    find_package(METIS REQUIRED)
    target_link_libraries(${target} PRIVATE ${METIS_LIBRARIES})
    target_include_directories(${target} PRIVATE ${METIS_INCLUDE_DIRS})
  endif()
  target_compile_definitions(${target} PRIVATE PZ_USING_METIS)
  target_compile_definitions(${target} INTERFACE PZ_USING_METIS)
endfunction()
