file(GLOB scripts *.sh)

# TODO: Add dependency to PerfTests/SubStruct

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/run_cubo.sh"
                "${CMAKE_CURRENT_BINARY_DIR}/run_cubo.sh")

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/run_cubo_threads.sh"
                "${CMAKE_CURRENT_BINARY_DIR}/run_cubo_threads.sh")

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/run_cubo_threads_binckpt.sh"
                "${CMAKE_CURRENT_BINARY_DIR}/run_cubo_threads_binckpt.sh")

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/pthread_vs_tbb-cubo_binckpt.sh"
                "${CMAKE_CURRENT_BINARY_DIR}/pthread_vs_tbb-cubo_binckpt.sh")

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/pthread_vs_tbb-predio_md5.sh"
                "${CMAKE_CURRENT_BINARY_DIR}/pthread_vs_tbb-predio_md5.sh")

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/run_predio_threads.sh"
                "${CMAKE_CURRENT_BINARY_DIR}/run_predio_threads.sh")

