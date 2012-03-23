# @file neopz/Refine/RefPatterns/CMakeLists.txt  -- CMake file for the refpatterns files

file(GLOB refpatterns *.rpt)
file(GLOB text *.txt)
install(FILES ${refpatterns} ${text} DESTINATION ${NEOPZ_INCLUDE}/refpatterns)

