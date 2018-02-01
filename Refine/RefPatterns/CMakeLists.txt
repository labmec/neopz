# @file neopz/Refine/RefPatterns/CMakeLists.txt  -- CMake file for the refpatterns files

file(GLOB refpatterns *.rpt)
file(GLOB text *.txt)
install(FILES ${refpatterns} ${text} DESTINATION ${NEOPZ_INCLUDE}/Refine/RefPatterns)

#source_group(RefPatterns FILES ${refpatternsVar})

if(BUILD_TUTORIAL)

	add_custom_target(REFPATTERNS_Files SOURCES ${refpatterns})

endif(BUILD_TUTORIAL)
