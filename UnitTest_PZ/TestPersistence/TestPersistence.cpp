/**
 * @file BlendUnitTest.cpp
 * @brief Define a Unit Test for testing features related to writing to and reading
 from files.
 *
 */

#include <catch2/catch_test_macros.hpp>

#include "pzfmatrix.h"
#include "tpzautopointer.h"

TEST_CASE("MatrixWriteTest","[persistence_tests]") {
  //create any dummy variable to be written on disk
  using namespace std::complex_literals;
  TPZFMatrix<CSTATE> mymat = {{
      {1.0,1.0+1i,1.0+2i},
      {2.0,2.0+1i,2.0+2i},
      {3.0,3.0+1i,3.0+2i}}};

  const std::string mat_file("testfile.pz");
  /*
    TPZPersistenceManager::PopulateClassIdMap is a fundamental function
    for writing/reading from files in NeoPZ, as the identifiers of NeoPZ
    classes are registered in a map, allowing for identifying instances
    when reading from a file. Therefore, we just test opening and closing a file.
  */
  SECTION("PopulateMap"){
    TPZPersistenceManager::OpenWrite(mat_file);
    TPZPersistenceManager::CloseWrite();
  }

  SECTION("WriteToFile"){
    TPZPersistenceManager::OpenWrite(mat_file);
    TPZPersistenceManager::WriteToFile(&mymat);
    TPZPersistenceManager::CloseWrite();
  }

  SECTION("ReadFromFile"){
    TPZPersistenceManager::OpenRead(mat_file);
    TPZAutoPointer<TPZFMatrix<CSTATE>> mymat_file =
      dynamic_cast<TPZFMatrix<CSTATE>*>(TPZPersistenceManager::ReadFromFile());
    TPZPersistenceManager::CloseRead();

    SECTION("MatData"){
      REQUIRE(mymat_file);
      REQUIRE(mymat_file->Rows() == mymat.Rows());
      REQUIRE(mymat_file->Cols() == mymat.Cols());
    }
    SECTION("MatVals"){
      const int nrows = mymat.Rows();
      const int ncols = mymat.Cols();
      for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
          REQUIRE(mymat_file->GetVal(i,j) == mymat.GetVal(i,j));
        }
      }
    }
  }
}