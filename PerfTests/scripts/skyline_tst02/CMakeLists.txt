file(GLOB scripts *.py)

# TODO: Add dependency to PerfTests/prog/skyline

file(COPY __init__.py DESTINATION .)

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/test.py"
                "${CMAKE_CURRENT_BINARY_DIR}/test.py")	

