//
//  PerkinsGonzalez.cpp
//  PZ
//
//  Created by Philippe Devloo on 3/25/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "PerkinsGonzalez.h"
#include "tpzautopointer.h"
#include "TPZGenGrid2D.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzelasAXImat.h"
#include "pzelasmat.h"
#include "pzbndcond.h"

#include "TPZVTKGeoMesh.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

TPZAutoPointer<TPZCompMesh> BuildCompMesh(TPZAutoPointer<TPZGeoMesh> gmesh);

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.project.perkinsgonzalez");
#endif

void forceX(TPZVec<REAL> &x, TPZVec<REAL> &force)
{
    if(x[0] > 2 && x[0] < 12) force[0] = 1.;
    else force[0] = 0.;
}

void TemperatureFunction(const TPZVec<REAL> &rz, REAL &temp)
{
    if(rz[0] < 20) temp = 100.;
    else temp = 0.;
    
}

int main()
{
//    TPZGenGrid2D(TPZVec<int> &nx, TPZVec<REAL> &x0, TPZVec<REAL> &x1, int numl = 1, REAL rot = 0.5);

    const int nel=2;
    TPZVec<int> nx(2,nel);
    nx[1] = 1;
    TPZVec<REAL> x0(3,0.),x1(3,2.);
    x0[0] = 1.;
    x1[1] = 1.;
    TPZGenGrid2D gengrid(nx,x0,x1);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gengrid.Read(gmesh);
    gengrid.SetBC(gmesh,4,-1);
    gengrid.SetBC(gmesh,6,-2);
    TPZAutoPointer<TPZCompMesh> cmesh = BuildCompMesh(gmesh);
    
    std::ofstream gmeshfile("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh.operator->(), gmeshfile,true);
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        gmesh->Print(sout);
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZSkylineStructMatrix skylstruct(cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    TPZAnalysis an(cmesh);
    an.SetStructuralMatrix(skylstruct);
    an.SetSolver(step);
    an.Run();
    
    //an.Solution().Print("Solution",std::cout);
    
    TPZStack<std::string> scalnames,vecnames;
    /*
    if(!strcmp("Sigmarr",name.c_str()))       return 3;
    if(!strcmp("Sigmazz",name.c_str()))       return 4;
    if(!strcmp("Sigmatt",name.c_str()))       return 5;
    if(!strcmp("Taurz",name.c_str()))         return 6;
     */
    scalnames.Push("SigmaX");
    scalnames.Push("SigmaY");
    scalnames.Push("TauXY");
    scalnames.Push("Pressure");
    vecnames.Push("displacement");
    an.DefineGraphMesh(2, scalnames, vecnames, "perkins.vtk");
    int postprocessresolution = 2;
    an.PostProcess(postprocessresolution);
    std::cout << "Finished\n";
    return 0;

}

TPZAutoPointer<TPZCompMesh> BuildCompMesh(TPZAutoPointer<TPZGeoMesh> gmesh)
{

    TPZManVector<int64_t> nodeindexes(1,0);
    int64_t index;
    int pointbc(-3);
    gmesh->CreateGeoElement(EPoint, nodeindexes, pointbc, index);
    gmesh->BuildConnectivity();
    TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
    REAL Elast = 1000.;
    REAL nu = 0.3;
    REAL fx(0.),fy(0.);
    TPZElasticityMaterial *aximat = new TPZElasticityMaterial(1,Elast,nu,fx,fy);
    //aximat->SetForcingFunction(forceX);
    //aximat->SetTemperature(100.);
    //aximat->SetTemperatureFunction(TemperatureFunction);
    TPZMaterial * autoaximat(aximat);
    cmesh->InsertMaterialObject(autoaximat);
    int mixed = 2;
    TPZFMatrix<STATE> val1(2,2,0.),val2(2,1,0.);
    val1(1,1) = 100.;
    val1(0,0) = 100.;
    TPZBndCond *bnd = aximat->CreateBC(autoaximat, pointbc, mixed, val1, val2);
    TPZMaterial * autobnd(bnd);
    cmesh->InsertMaterialObject(autobnd);
    val1.Zero();
    val2(0,0) = 10.;
    int neumann = 1;
    int rightbc = -2;
    TPZBndCond *bnd2 = aximat->CreateBC(autoaximat, rightbc, neumann, val1, val2);
    TPZMaterial * autobnd2(bnd2);
    cmesh->InsertMaterialObject(autobnd2);
    val2(0,0) = -10.;
    int leftbc = -1;
    TPZBndCond *bnd3 = aximat->CreateBC(autoaximat, leftbc, neumann, val1, val2);
    TPZMaterial * autobnd3(bnd3);
    cmesh->InsertMaterialObject(autobnd3);
    int porder = 3;
    cmesh->SetDefaultOrder(porder);
    cmesh->AutoBuild();
    return cmesh;
}
