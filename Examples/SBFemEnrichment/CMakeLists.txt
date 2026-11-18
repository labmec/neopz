cmake_minimum_required(VERSION 3.14)

# The same examples work in the NeoPZ tree or as an installed-package consumer.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(SBFemEnrichmentExamples LANGUAGES CXX)
    find_package(NeoPZ CONFIG REQUIRED)
    if(NOT NeoPZ_SBFEM_ENRICHMENT)
        message(FATAL_ERROR "Install NeoPZ with USING_LAPACK=ON and USING_SBFEM_ENRICHMENT=ON")
    endif()
    set(enrichment_library NeoPZ::pz)
    enable_testing()
else()
    set(enrichment_library pz)
endif()

add_executable(sbfem_eigenvalues eigenvalues.cpp)
add_executable(sbfem_evaluate evaluate.cpp)
foreach(example sbfem_eigenvalues sbfem_evaluate)
    target_compile_features(${example} PRIVATE cxx_std_17)
    target_link_libraries(${example} PRIVATE ${enrichment_library})
endforeach()

foreach(config isotropic orthotropic bimaterial)
    configure_file(${config}.json ${config}.json COPYONLY)
endforeach()

add_test(NAME SBFemExampleIsotropicEigenvalues COMMAND sbfem_eigenvalues
    "${CMAKE_CURRENT_BINARY_DIR}/isotropic.json" --verify-isotropic)
add_test(NAME SBFemExampleValuesAndGradients COMMAND sbfem_evaluate
    "${CMAKE_CURRENT_BINARY_DIR}/isotropic.json" 0.35 0.20 --verify-gradients)
add_test(NAME SBFemExampleOrthotropic COMMAND sbfem_evaluate
    "${CMAKE_CURRENT_BINARY_DIR}/orthotropic.json" 0.35 0.20 --verify-gradients)
add_test(NAME SBFemExampleBimaterial COMMAND sbfem_evaluate
    "${CMAKE_CURRENT_BINARY_DIR}/bimaterial.json" 0.35 0.20 --verify-gradients)
add_test(NAME SBFemExampleRejectTip COMMAND sbfem_evaluate
    "${CMAKE_CURRENT_BINARY_DIR}/isotropic.json" 0 0)
set_tests_properties(SBFemExampleRejectTip PROPERTIES WILL_FAIL TRUE)
set_tests_properties(SBFemExampleIsotropicEigenvalues SBFemExampleValuesAndGradients
    SBFemExampleOrthotropic SBFemExampleBimaterial SBFemExampleRejectTip
    PROPERTIES LABELS "sbfem;enrichment;example" TIMEOUT 180)
