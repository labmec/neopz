# Reuse NeoPZ's Catch2 runner; no second test_library or GFEM checkout needed.
add_unit_test(Test_TPZGFEMEnrichment test_TPZGFemEnrichment.cpp)
add_unit_test(Test_TPZGFEMMaterials test_TPZGFemMaterials.cpp)
add_unit_test(Test_TPZGFEMRegression test_TPZGFemRegression.cpp)

set_tests_properties(
    Test_TPZGFEMEnrichment Test_TPZGFEMMaterials Test_TPZGFEMRegression
    PROPERTIES LABELS "sbfem;enrichment" TIMEOUT 180)
