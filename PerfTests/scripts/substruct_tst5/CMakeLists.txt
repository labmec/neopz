file(GLOB scripts *.py)

# TODO: Add dependency to PerfTests/prog/substruct

file(COPY __init__.py DESTINATION .)

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/test.py"
                "${CMAKE_CURRENT_BINARY_DIR}/test.py")	

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/test6.py"
                "${CMAKE_CURRENT_BINARY_DIR}/test6.py")

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/test12.py"
                "${CMAKE_CURRENT_BINARY_DIR}/test12.py")

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/test18.py"
                "${CMAKE_CURRENT_BINARY_DIR}/test18.py")

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/test18.py"
                "${CMAKE_CURRENT_BINARY_DIR}/test18.py")