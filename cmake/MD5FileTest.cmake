#message(STATUS "Filename='${filename}'")
#message(STATUS "ExpectedMD5=${expectedMD5}###")

#calculating MD5
execute_process(
  COMMAND ${CMAKE_COMMAND}
    #-D "dir:STRING=${dir}"
    #-D "testname:STRING=${testname}"
    #-P "${scriptname}"
    -E md5sum ${filename}
  OUTPUT_VARIABLE out
  ERROR_VARIABLE err
  RESULT_VARIABLE result
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_STRIP_TRAILING_WHITESPACE
  )

#getting MD5 from the string
STRING(SUBSTRING "${out}" 0 32 calculatedMD5)
MESSAGE ("Calculated MD5:${calculatedMD5}###")

#checking
if(NOT "${calculatedMD5}" STREQUAL "${expectedMD5}")
    message(FATAL_ERROR "Wrong MD5 for '${filename}', expected: '${expectedMD5}', got: '${calculatedMD5}' !!")
endif()


