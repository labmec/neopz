file(GLOB scripts *.sh)

# TODO: Add dependency to PerfTests/SubStruct

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/run_cubo.sh"
                "${CMAKE_CURRENT_BINARY_DIR}/run_cubo.sh")
