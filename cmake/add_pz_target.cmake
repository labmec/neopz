#
# Create a target dependent on NeoPZ giving its name, its source files, files to be copied to the binary dir 
# and required NeoPZ CMake options.
# Usage:
#     add_pz_target(
#       NAME myTarget
#       SOURCES source.cpp header.h
#       FILES [optional] files that will be copied to target binary dir, e.g.: QuadMesh.msh, input.json, etc
#       REQUIRED [optional] compile definitions that NeoPZ in order to run 'myTarget', e.g.: PZ_USING_MKL, PZ_LOG, etc
#       )
#


include (CMakeParseArguments)

function(add_pz_target)
  cmake_parse_arguments(
    PARSED_ARGS # prefix of output variables
    "" # list of names of the boolean arguments (only defined ones will be true)
    "NAME" # list of names of mono-valued arguments
    "SOURCES;FILES;REQUIRED" # list of names of multi-valued arguments (output variables are lists)
    ${ARGN} # arguments of the function to parse, here we take the all original ones
    )
  # note: if it remains unparsed arguments, here, they can be found in variable PARSED_ARGS_UNPARSED_ARGUMENTS
  if(NOT PARSED_ARGS_NAME)
    message(FATAL_ERROR "You must provide a name for the target")
  endif(NOT PARSED_ARGS_NAME)
  list(LENGTH PARSED_ARGS_SOURCES N_SOURCES)
  if (N_SOURCES LESS 1)
    message(FATAL_ERROR "You must provide sources for the target")
  endif(N_SOURCES LESS 1)
  add_executable(${PARSED_ARGS_NAME} "")
  target_sources(${PARSED_ARGS_NAME} PRIVATE ${PARSED_ARGS_SOURCES})
  if(CMAKE_IS_PZ_BUILDTREE)
    target_link_libraries(${PARSED_ARGS_NAME} PRIVATE pz)
  else()
    target_link_libraries(${PARSED_ARGS_NAME} PRIVATE NeoPZ::pz)
    target_include_directories(${PARSED_ARGS_NAME} PRIVATE ${PZ_INCLUDE_DIRS})
  endif()
  foreach(file ${PARSED_ARGS_FILES})
    configure_file(${file} ${file} COPYONLY)
  endforeach(file)

  foreach(reqr ${PARSED_ARGS_REQUIRED})
    check_pz_opt(${reqr} res)
    if(NOT ${res})
      message(FATAL_ERROR "This target requires option ${reqr} in the NeoPZ library and it is not set")
    endif()
  endforeach(reqr)  
endfunction(add_pz_target)
