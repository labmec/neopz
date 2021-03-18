function(enable_threads target)
    #prefer pthread library when there are other options available
    set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
    #prefer -pthread compiler and linker flag
    set(THREADS_PREFER_PTHREAD_FLAG TRUE)
    find_package(Threads REQUIRED)
    #Public linkage due to ParallelFor function
    target_link_libraries(${target} PUBLIC ${CMAKE_THREAD_LIBS_INIT})
    #simple utility function for detecting if any linking is needed
    #for using functions from <atomic>
    include(cmake/check_for_atomic.cmake)
    check_for_atomic(${target})
endfunction()
