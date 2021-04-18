# Here we will search for Catch2. If eigen is not found, we will download it.
function(enable_Catch2)
    #perhaps catch was already downloaded when running cmake
    if(NOT TARGET Catch2::Catch2)
        find_package(Catch2 QUIET)
        if(NOT TARGET Catch2::Catch2)
            Include(FetchContent)
            if(NOT Catch2_SOURCE_DIR)
                message(STATUS "Downloading Catch2")
            endif()
            FetchContent_Declare(
                Catch2
                GIT_REPOSITORY https://github.com/catchorg/Catch2.git
                GIT_TAG        v2.13.6)
            FetchContent_MakeAvailable(Catch2)
        endif()
    endif()
    if(TARGET Catch2::Catch2)
        message(STATUS "Catch2 found at ${Catch2_SOURCE_DIR}")
    else()
        message(FATAL_ERROR "Could not satisfy dependency: Catch2")
    endif()
endfunction()