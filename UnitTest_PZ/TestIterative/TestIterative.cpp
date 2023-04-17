/**
 * @file TestIterative.cpp
 * @brief Define a Unit Test using Catch2 for usage of preconditioners in NeoPZ
 *
 */

#include "pzcmesh.h"
#include "TPZGeoMeshTools.h"
#include "TPZLinearAnalysis.h"
#include "pzstepsolver.h"
#include "TPZIdentitySolver.h"
#include "Poisson/TPZMatPoisson.h"
#include "TPZSSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "TPZVTKGenerator.h"
#include <iostream>


#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

template<>
struct Catch::StringMaker<long double> {
    static std::string convert(long double ref);
    static int precision;
};

int Catch::StringMaker<long double>::precision = 10;

std::string Catch::StringMaker<long double>::convert(long double value) {
    std::ostringstream out;
    out.precision(precision);
    out << std::fixed << value;
    return out.str();
}

namespace testiterative{
  //! Creates computational mesh for 2d poisson problem
  TPZCompMesh *CreateCMesh();
}

TEST_CASE("Poisson equation","[iterative_testss]") {
  auto oldPrecision = Catch::StringMaker<STATE>::precision;
  Catch::StringMaker<STATE>::precision = std::numeric_limits<STATE>::max_digits10;
  

  TPZAutoPointer<TPZCompMesh> cmesh = testiterative::CreateCMesh();
  constexpr bool optimizeBandwidth{false};
  TPZLinearAnalysis an(cmesh,optimizeBandwidth);
  {
    TPZStepSolver<STATE> solv;
    solv.SetDirect(ECholesky);
    an.SetSolver(solv);
  }
  TPZSSpStructMatrix<STATE> str(cmesh);
  an.SetStructuralMatrix(str);
  an.Assemble();

  auto *solver = dynamic_cast<TPZStepSolver<STATE>*>(an.Solver());

  solver->Matrix()->SetSymmetry(SymProp::Herm);
  solver->Matrix()->SetDefPositive(true);



  constexpr int niter{40};
  constexpr REAL tol{1e-12};
  constexpr int64_t fromCurrent{0};
  SECTION("Identity matrix as precond"){
    std::cout<<"testing identity matrix as a precond"<<std::endl;
    bool success = false;
    TPZIdentitySolver<STATE> pre;
    solver->SetCG(niter, pre, tol, fromCurrent);
    an.Solve();
    success = true;
    REQUIRE(success);
  }
  SECTION("Exact inverse as precond"){
    std::cout<<"testing exact inverse as precond"<<std::endl;
    //we copy the original matrix
    TPZAutoPointer<TPZMatrix<STATE>> mtrxcp = solver->Matrix()->Clone();
    //we use a direct solver as a precond
    TPZStepSolver<STATE> pre(mtrxcp);
    pre.SetDirect(ECholesky);
    solver->SetCG(niter, pre, tol, fromCurrent);
    an.Solve();
    const auto niter_res = solver->NumIterations();
    REQUIRE(niter_res == 1);
    const auto tol_res = solver->GetTolerance();
    REQUIRE(tol_res <= tol);
  }

  auto pre_type = GENERATE(Precond::Jacobi,
                           Precond::BlockJacobi,
                           Precond::Element,
                           Precond::NodeCentered);

  const bool overlap = GENERATE(true,false);
  const auto precondname = std::string(Precond::Name(pre_type)) +
                                       (overlap ? " with overlap" : " without overlap");
  SECTION(precondname){
    std::cout<<"testing "<<precondname<<std::endl;
    TPZAutoPointer<TPZMatrixSolver<STATE>> pre =
      an.BuildPreconditioner<STATE>(pre_type, overlap);
    const bool created_precond = pre;
    REQUIRE(created_precond);
    solver->SetCG(niter, *pre, tol, fromCurrent);
    an.Solve();
    const auto tol_res = solver->GetTolerance();
    REQUIRE(tol_res <= tol);
  }
  Catch::StringMaker<STATE>::precision = oldPrecision;
}

namespace testiterative{

  TPZCompMesh *CreateCMesh()
  {
    constexpr int rhsPOrder{2};
    const auto rhs = [](const TPZVec<REAL>&loc, TPZVec<STATE> &u){
      const REAL &x = loc[0];
      const REAL &y = loc[1];
      u[0] = 2*y*y+2*x*x-4;
      u[0] *= -1;
    };
    
    //dimension of the problem
    constexpr int dim{2};
    //n divisions in x direction
    constexpr int nDivX{4};
    //n divisions in y direction
    constexpr int nDivY{4};
  
    //TPZManVector<Type,N> is a vector container with static + dynamic storage. one can also use TPZVec<Type> for dynamic storage
    TPZManVector<int,2> nDivs={nDivX,nDivY};

    //all geometric coordinates in NeoPZ are in the 3D space

    //lower left corner of the domain
    TPZManVector<REAL,3> minX={-1,-1,0};
    //upper right corner of the domain
    TPZManVector<REAL,3> maxX={ 1, 1,0};

    /*vector for storing materials(regions) identifiers
     * for the TPZGeoMeshTools::CreateGeoMeshOnGrid function.
     * In this example, we have different materials in each
     * section of the boundary*/
    TPZManVector<int,5> matIdVec={1,-1,-2,-3,-4};
    //whether to create boundary elements
    constexpr bool genBoundEls{true};
    //type of elements
    constexpr MMeshType meshType{MMeshType::ETriangular};


    //polynomial order of elements
    constexpr int pOrder{2};
    //defining the geometry of the problem
    //TPZAutoPointer is a smart pointer from the NeoPZ library
    TPZAutoPointer<TPZGeoMesh> gmesh =
      TPZGeoMeshTools::CreateGeoMeshOnGrid(dim,minX,maxX,matIdVec,nDivs,meshType,genBoundEls);

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

    auto *mat = new TPZMatPoisson<STATE>(matIdVec[0],dim);
    //TPZMatLaplacian solves div(k grad(u)) = -f
    mat->SetForcingFunction(rhs,rhsPOrder);
    cmesh->InsertMaterialObject(mat);

    //now we insert the boundary conditions
    for(auto i = 0; i < 4; i++)
    {
      //TPZFMatrix<T> implements a full matrix of type T
      /* val1 and val2 are used for calculating the boundary
       * conditions. val1 goes in the matrix and val2 in the rhs.
       * for dirichlet boundary conditions, only the value of 
       * val2 is used.*/
      TPZFMatrix<STATE> val1(1,1,0.);
      TPZManVector<STATE,1> val2={0};
      //dirichlet=0,neumann=1,robin=2
      constexpr int boundType{0};
      //TPZBndCond is a material type for boundary conditions
      TPZBndCond * bnd = mat->CreateBC(mat,matIdVec[i+1],boundType,val1,val2);
      cmesh->InsertMaterialObject(bnd);
    }
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(pOrder);
    cmesh->AutoBuild();
    return cmesh;
  }
}