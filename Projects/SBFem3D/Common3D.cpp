#include "Common3D.h"

#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include "TPZMatElasticity2D.h"
#include "TPZMatLaplacian.h"
#include "pzbndcond.h"

#include "TPZAcademicGeoMesh.h"
#include "pzgengrid.h"
#include "TPZBuildSBFem.h"

#include "TPZVTKGeoMesh.h"

#include "TPZSBFemElementGroup.h"
#include "pzinterpolationspace.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"

TLaplaceExampleSmooth ExactLaplace;


void SolveSist(TPZAnalysis *an, TPZCompMesh *Cmesh)
{
    //    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(Cmesh);
    TPZSkylineStructMatrix strmat(Cmesh);
    //    TPZSymetricSpStructMatrix strmat(Cmesh);
    strmat.SetNumThreads(0);
    an->SetStructuralMatrix(strmat);
    
    long neq = Cmesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an->SetSolver(step);
    
    an->Assemble();
    
    //    std::ofstream andrade("../Andrade.mtx");
    //    andrade.precision(16);
    //    an->Solver().Matrix()->Print("Andrade",andrade,EMatrixMarket);
    //    std::cout << "Leaving Assemble\n";
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
    
    if(neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
    
    an->Solve();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2-t1 << " Time for solving " << t3-t2 << std::endl;
#endif
    
    
}

void HarmonicNeumannLeft(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = -M_PI*exp(M_PI*x[0])*sin(M_PI*x[1]);
}

void HarmonicNeumannRight(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    val[0] = M_PI*exp(M_PI*x[0])*sin(M_PI*x[1]);
}

void Harmonic_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    val[0] = exp(M_PI*xv[0])*sin(M_PI*xv[1]);
    deriv(0,0) = M_PI*val[0];
    deriv(1,0) = M_PI*exp(M_PI*xv[0])*cos(M_PI*xv[1]);
    
}


void InsertMaterialObjects3D(TPZCompMesh *cmesh, int problemtype)
{
    // Plane strain assumption
    int planestress = 0;
    
    if(problemtype != 0 && problemtype != 1) DebugStop();
    // Getting mesh dimension
    int dim = 2;
    int matId1 = Emat1;
    
    TPZMaterial *material;
    int nstate = 1;
    bool elasticity = false;
    if (problemtype == 0) {
        elasticity = true;
    }
    if (elasticity)
    {
        TPZMatElasticity2D *matloc = new TPZMatElasticity2D(matId1);
        material = matloc;
        nstate = 2;
        // Setting up paremeters
        matloc->SetfPlaneProblem(planestress);
        //        REAL lamelambda = 1.0e9,lamemu = 0.5e3, fx= 0, fy = 0;
        REAL lamelambda = 0.,lamemu = 0.5e3, fx= 0, fy = 0;
        matloc->SetParameters(lamelambda,lamemu, fx, fy);
        //material->SetElasticParameters(40.0,0.0);
        REAL Sigmaxx = 0.0, Sigmayx = 0.0, Sigmayy = 0.0, Sigmazz = 0.0;
        matloc->SetPreStress(Sigmaxx,Sigmayx,Sigmayy,Sigmazz);
    }
    else
    {
        TPZMatLaplacian *matloc = new TPZMatLaplacian(matId1);
        matloc->SetForcingFunction(ExactLaplace.ForcingFunction());
        matloc->SetDimension(3);
        matloc->SetSymmetric();
        material = matloc;
        nstate = 1;
    }
    //material->SetBiotAlpha(Alpha);cade o metodo?
    
    
    TPZFMatrix<STATE> val1(nstate,nstate,0.), val2(nstate,1,0.);
    val2(0,0) = 0.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BCond1;
    if(elasticity==0)
    {
        BCond1 = material->CreateBC(material,Ebc1,0, val1, val2);
        BCond1->SetForcingFunction(ExactLaplace.Exact());
    }
    else
    {
        BCond1 = material->CreateBC(material,Ebc1,1, val1, val2);
    }
    
    
    val2(0,0) = 0.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BSkeleton = material->CreateBC(material,ESkeleton,1, val1, val2);
    
    
    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BSkeleton);
    
}

TPZCompMesh *SetupSquareMesh3D(int nelx, int nrefskeleton, int porder, bool elasticityproblem)
{
    
    TPZAcademicGeoMesh acadgmesh(nelx,TPZAcademicGeoMesh::EHexa);
    TPZManVector<int,6> bcids(6,-1);
    acadgmesh.SetBCIDVector(bcids);
    acadgmesh.SetMaterialId(EGroup);
    
    TPZAutoPointer<TPZGeoMesh> gmesh = acadgmesh.CreateGeoMesh();
    
    
    std::map<int,int> matmap;
    matmap[EGroup] = 1;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    
    build.StandardConfiguration();
    build.DivideSkeleton(nrefskeleton);
    //        AddSkeletonElements(gmesh);
    /// generate the SBFem elementgroups
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 0;
    if (elasticityproblem) {
        problemtype = 0;
    }
    else
    {
        problemtype = 1;
    }
    InsertMaterialObjects3D(SBFem,problemtype);
    
    
    build.BuildComputationMesh(*SBFem);
    
    if(1)
    {
        std::ofstream outg("GMesh3D.txt");
        gmesh->Print(outg);
        std::ofstream out("Geometry3D.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return SBFem;
}



