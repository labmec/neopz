function(enable_neopzexamples)
    Include(FetchContent)
    FetchContent_Declare(
        NeoPZExamples
        GIT_REPOSITORY https://github.com/labmec/NeoPZExamples.git
        GIT_TAG        main
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/Projects
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/Projects
        BINARY_DIR ${PROJECT_BINARY_DIR}/Projects)
    message(STATUS "Checking NeoPZExamples...")
    FetchContent_MakeAvailable(NeoPZExamples)
    message(STATUS "NeoPZExamples found at ${NeoPZExamples_SOURCE_DIR}")
endfunction()