file(GLOB scripts *.sh)

# TODO: Add dependency to PerfTests/SubStruct

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/substruct_perf01.sh"
                "${CMAKE_CURRENT_BINARY_DIR}/substruct_perf01.sh")	
