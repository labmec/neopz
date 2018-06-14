if(BUILD_PYTHON_BINDING)
    message (STATUS "SWIG 2 and Python found. Python module will built")
	target_include_directories(pz PUBLIC 
							   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						       #$<INSTALL_INTERFACE:.> 
						      )

    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(SWSIZE SWIGWORDSIZE64)
    else(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(SWSIZE SWIGWORDSIZE32)
    endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

    add_custom_command(OUTPUT TPZPython_wrap.cxx ${CMAKE_CURRENT_BINARY_DIR}/TPZPython.py
	COMMAND ${SWIG_EXECUTABLE} -c++ -Wall -python -D${SWSIZE} -o ${CMAKE_CURRENT_BINARY_DIR}/TPZPython_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} TPZPython.i
	COMMENT "Building SWIG source in ${SWSIZE} mode: See swigpython.out for errors and warning"
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
        
    set(TPZPython_src ${CMAKE_CURRENT_BINARY_DIR}/TPZPython_wrap.cxx)

    include_directories(${PYTHON_INCLUDE_PATH})

    add_library(TPZPython SHARED ${TPZPython_src} ${pzsources})

    IF (APPLE OR WIN32)
        target_link_libraries(TPZPython pz ${PYTHON_LIBRARY}) #python)
    ELSE()
        target_link_libraries(TPZPython pz)
    ENDIF()

    IF (NOT WIN32)
		SET_TARGET_PROPERTIES(TPZPython PROPERTIES PREFIX "_")
		SET_TARGET_PROPERTIES(TPZPython PROPERTIES SUFFIX ".so")
	ELSE()
		SET_TARGET_PROPERTIES(TPZPython PROPERTIES PREFIX "_")
		SET_TARGET_PROPERTIES(TPZPython PROPERTIES SUFFIX ".pyd")
    ENDIF()
    
endif()
