/**
    \file poisson.cpp
    How to solve the poisson equation in a bidimensional domain using NeoPZ
*/
#include "pzgmesh.h" //for TPZGeoMesh
#include "pzcmesh.h" //for TPZCompMesh
#include "TPZGeoMeshTools.h" //for TPZGeoMeshTools::CreateGeoMeshOnGrid
#include "MMeshType.h" //for MMeshType
#include "pzmanvector.h"//for TPZManVector
#include "Poisson/TPZMatPoisson.h" //for TPZMatPoisson
#include "TPZBndCond.h" //for TPZBndCond
#include "TPZLinearAnalysis.h" //for TPZLinearAnalysis
#include "TPZSSpStructMatrix.h" //symmetric sparse matrix storage
#include "pzskylstrmatrix.h" //symmetric skyline matrix storage
#include "pzstepsolver.h" //for TPZStepSolver
#include <catch2/catch.hpp>

void CheckErrorPoisson(const int pOrder, TPZVec<REAL>&error);


TEST_CASE("error_poisson","[error_tests]")
{
  TPZManVector<REAL,3> lastError(3,10000);
  for(int p = 1; p < 5; p++){
    TPZManVector<REAL,3> error(3,0.);
    CheckErrorPoisson(p,error);
    bool checkConv = true;
    for(auto i = 0; i < 3; i++){
      if(error[i]>lastError[i]) checkConv = false;
    }
    REQUIRE(checkConv);
    lastError=error;
  }
  auto oldPrecision = Catch::StringMaker<REAL>::precision;
  Catch::StringMaker<REAL>::precision = std::numeric_limits<REAL>::max_digits10;
  //the solution is contained in the p=4 approx space
  for(auto ier : lastError){
    CAPTURE(ier);
    REQUIRE(ier == Approx(0.0).margin(std::numeric_limits<REAL>::epsilon()*10));
  }
  Catch::StringMaker<REAL>::precision = oldPrecision;
}

void CheckErrorPoisson(const int pOrder, TPZVec<REAL>&error)
{
  
  /* We will solve div(grad(u)) = 2y^2+2x^2-4
   * in the domain Omega=[-1,1]x[-1,1] embedded in a 3D space
   * with homogeneous boundary conditions u=0 on dOmega*/
  constexpr int solOrder{4};
  auto exactSol = [](const TPZVec<REAL> &loc,
         TPZVec<STATE>&u,
         TPZFMatrix<STATE>&gradU){
        const auto &x=loc[0];
        const auto &y=loc[1];
        u[0]=(x*x-1)*(y*y-1);
        gradU(0,0) = 2*x*(y*y-1);
        gradU(1,0) = 2*y*(x*x-1);
        gradU(2,0) = 0;//optional
  };
  constexpr int rhsPOrder{2};
  const auto rhs = [](const TPZVec<REAL>&loc, TPZVec<STATE> &u){
        const REAL &x = loc[0];
        const REAL &y = loc[1];
        u[0] = 2*y*y+2*x*x-4;
        //see comment on TPZMatPoisson
        u[0] *= -1;
  };
  
  //dimension of the problem
  constexpr int dim{2};
  //n divisions in x direction
  constexpr int nDivX{2};
  //n divisions in y direction
  constexpr int nDivY{2};
  
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

  //defining the geometry of the problem
  //TPZAutoPointer is a smart pointer from the NeoPZ library
  TPZAutoPointer<TPZGeoMesh> gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim,minX,maxX,matIdVec,nDivs,meshType,genBoundEls);
  ///Defines the computational mesh based on the geometric mesh
  TPZAutoPointer<TPZCompMesh>  cmesh = new TPZCompMesh(gmesh);
  //using traditional H1 elements
  cmesh->SetAllCreateFunctionsContinuous();

  /* The TPZMaterial class is used for implementing the weak formulation.
   * Each instance has an associated material id, which should correspond to the
   * material ids used when creating the geometric mesh. In this way, you could
   * have different materials on different mesh regions */

  auto *mat = new TPZMatPoisson(matIdVec[0],dim);

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
      TPZVec<STATE> val2(1,0.);
      //dirichlet=0,neumann=1,robin=2
      constexpr int boundType{0};
      //TPZBndCond is a material type for boundary conditions
      TPZBndCond * bnd = mat->CreateBC(mat,matIdVec[i+1],boundType,val1,val2);
      cmesh->InsertMaterialObject(bnd);
    }
  
  //seting the polynomial order in the computational mesh
  cmesh->SetDefaultOrder(pOrder);
  //actually creates the computational elements
  cmesh->AutoBuild();

  /*The TPZLinearAnalysis class manages the creation of the algebric
  * problem and the matrix inversion*/
  TPZLinearAnalysis an(cmesh);

  //sets number of threads to be used by the solver
  constexpr int nThreads{4};
  //defines storage scheme to be used for the FEM matrices
  //in this case, a symmetric skyline matrix is used
  TPZSkylineStructMatrix<STATE> matskl(cmesh);
  matskl.SetNumThreads(nThreads);
  an.SetStructuralMatrix(matskl);
  	
  ///Setting a direct solver
  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt);
  an.SetSolver(step);

  //assembles the system
  an.Assemble();
	
  ///solves the system
  an.Solve();
  //let us set the exact solution and suggest an integration rule
  //for calculating the error
  an.SetExact(exactSol,solOrder);
  an.PostProcess(error);
}

