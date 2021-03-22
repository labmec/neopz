function(find_include_from_mkl incvar dirname incname)
  if(WIN32)
    set(SEARCH_ROOT_DIR "C:")
  else()
    set(SEARCH_ROOT_DIR "/")
  endif()
  while(NOT ${dirname} STREQUAL ${SEARCH_ROOT_DIR})
    find_path(${incvar}
      NAMES ${incname}
      PATHS ${dirname}
      PATH_SUFFIXES "include" "latest/include" "include/mkl"
      NO_DEFAULT_PATH #let us try not to mix installs
      )
    if(NOT tmp_incvar)
      #one level up...
      get_filename_component(tmp_var ${dirname} DIRECTORY)
      set(dirname ${tmp_var})
    else()
      break()
    endif()
  endwhile()
  set(${incvar} ${${incvar}} PARENT_SCOPE)
endfunction()
