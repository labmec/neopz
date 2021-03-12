#
# Create a target dependent on NeoPZ with possible linkage to extra libraries, copied files to binary dir variable number of sources.
# This function is meant to be used for targets IN THIS PROJECT.
#
# Usage:
#     add_pz_internal_target(
#       NAME myTarget
#       SOURCES source.cpp header.h
#       FILES [optional] files that will be copied to target binary dir
#


include (CMakeParseArguments)

function(add_pz_internal_target)
  cmake_parse_arguments(
    PARSED_ARGS # prefix of output variables
    "" # list of names of the boolean arguments (only defined ones will be true)
    "NAME" # list of names of mono-valued arguments
    "SOURCES;FILES" # list of names of multi-valued arguments (output variables are lists)
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
  target_link_libraries(${PARSED_ARGS_NAME} PRIVATE pz)
  foreach(file ${PARSED_ARGS_FILES})
    configure_file(${file} ${file} COPYONLY)
  endforeach(file)  
endfunction(add_pz_target)
