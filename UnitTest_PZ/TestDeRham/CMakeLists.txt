add_unit_test(TestDeRham
    TestDeRham.cpp
    #materials for checking dimensions of span/kernel
    TPZMatDeRhamH1.cpp
    TPZMatDeRhamHCurl.cpp
    TPZMatDeRhamHDiv.cpp
    TPZMatDeRhamL2.cpp
    #materials for testing inclusion
    TPZMatDeRhamH1HCurl.cpp
    TPZMatDeRhamHDivL2.cpp
    TPZMatDeRhamHCurlHDiv.cpp
    )