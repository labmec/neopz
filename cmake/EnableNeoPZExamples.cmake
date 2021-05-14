function(enable_neopzexamples)
    Include(FetchContent)
    FetchContent_Declare(
        NeoPZExamples
        GIT_REPOSITORY https://github.com/labmec/NeoPZExamples.git
        GIT_TAG        main)
    message(STATUS "Checking NeoPZExamples...")
    FetchContent_MakeAvailable(NeoPZExamples)
    message(STATUS "NeoPZExamples found at ${NeoPZExamples_SOURCE_DIR}")
    set(CMAKE_IS_PZ_BUILDTREE ON PARENT_SCOPE)
endfunction()