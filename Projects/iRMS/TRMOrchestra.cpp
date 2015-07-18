//
//  TRMOrchestra.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMOrchestra.h"
#include "pzmaterial.h"
#include "TPZVecL2.h"
#include "pzl2projection.h"
#include "pzbndcond.h"
#include "TRMFlowConstants.h"


/** @brief Default constructor */
TRMOrchestra::TRMOrchestra(){
    
    /** @brief Define the global geometry being used */
    fgmesh = NULL;
    
//    /** @brief Define the mixed system analysis */
//    fFluxPressureAnalysis = NULL;
//    
//    /** @brief Define the transpor equation analysis */
//    fTransportAnalysis = NULL;
    
    /** @brief Define simulation data */
    fSimulationData = NULL;
    
}

/** @brief Default desconstructor */
TRMOrchestra::~TRMOrchestra(){
    
}

/** @brief Create a primal analysis using space odissey */
void TRMOrchestra::CreateAnalysisPrimal(TRMSpaceOdissey spacegenerator){
    
    TPZManVector<int,2> dx(2,1), dy(2,1), dz(2,1);
    dx[0] = 50.0;
    dy[0] = 50.0;
    dz[0] = 50.0;
    
    spacegenerator.CreateGeometricBoxMesh(dx, dy, dz);
//    spacegenerator.CreateGeometricReservoirMesh();
    spacegenerator.PrintGeometry();
    fgmesh = spacegenerator.GetGmesh();
    spacegenerator.CreateH1Cmesh();
    
    TPZAutoPointer<TPZCompMesh > Cmesh = spacegenerator.GetH1Cmesh();
    
    // Analysis
    bool mustOptimizeBandwidth = true;
    TPZAnalysis * AnalysisPrimal = new TPZAnalysis(Cmesh.operator->(),mustOptimizeBandwidth);
    int numofThreads = 8;
    
    TPZSkylineNSymStructMatrix skylnsym(Cmesh.operator->());
    TPZStepSolver<STATE> step;
    skylnsym.SetNumThreads(numofThreads);
    step.SetDirect(ELU);
    AnalysisPrimal->SetStructuralMatrix(skylnsym);
    AnalysisPrimal->SetSolver(step);;
    AnalysisPrimal->Run();
    std::cout << "Primal dof: " << AnalysisPrimal->Rhs().Rows() << std::endl;
    
    const int dim = 3;
    int div = 1;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "PrimalDarcy.vtk";
    scalnames.Push("Pressure");
    vecnames.Push("MinusKGradU");
    AnalysisPrimal->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    AnalysisPrimal->PostProcess(div);
    
}

/** @brief Create a dual analysis using space odissey */
void TRMOrchestra::CreateAnalysisDual(){
    
    TPZManVector<int,2> dx(2,1), dy(2,1), dz(2,1);
    dx[0] = 50.0;
    dy[0] = 50.0;
    dz[0] = 50.0;
    
//    spacegenerator.CreateGeometricBoxMesh(dx, dy, dz);
    fSpaceGenerator.CreateGeometricReservoirMesh();
#ifdef DEBUG
    fSpaceGenerator.PrintGeometry();
#endif
    
    fSpaceGenerator.CreateMixedCmesh();
    
    ProjectExactSolution();
    
    
    TPZAutoPointer<TPZCompMesh > Cmesh = fSpaceGenerator.GetMixedCmesh();
    
    // transfer the solution from the meshes to the multiphysics mesh
    TPZManVector<TPZAutoPointer<TPZCompMesh>,3 > meshvec(2);
    meshvec[0] = fSpaceGenerator.GetFluxCmesh();
    meshvec[1] = fSpaceGenerator.GetPressureMesh();
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, Cmesh);
    
    
    // Analysis
    bool mustOptimizeBandwidth = true;
    fFluxPressureAnalysis.SetCompMesh(Cmesh.operator->(), mustOptimizeBandwidth);
//    TPZAnalysis * AnalysisDual = new TPZAnalysis(Cmesh.operator->(),mustOptimizeBandwidth);
    int numofThreads = 0;
#ifdef DEBUG
    {
        std::ofstream out("../MFCompMesh.txt");
        Cmesh->Print(out);
    }
#endif
    
    TPZSkylineStructMatrix strmat(Cmesh.operator->());
//    TPZSkylineNSymStructMatrix strmat(Cmesh.operator->());
    TPZStepSolver<STATE> step;
    strmat.SetNumThreads(numofThreads);
    step.SetDirect(ELDLt);
    fFluxPressureAnalysis.SetStructuralMatrix(strmat);
    fFluxPressureAnalysis.SetSolver(step);
    fFluxPressureAnalysis.Assemble();
    fFluxPressureAnalysis.Rhs() *= -1.0;
    std::cout << "Dual dof: " << fFluxPressureAnalysis.Rhs().Rows() << std::endl;
    std::cout << "Rhs norm " << Norm(fFluxPressureAnalysis.Rhs()) << std::endl;
    fFluxPressureAnalysis.Solve();
    
    const int dim = 3;
    int div = 1;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "../DualDarcy.vtk";
    scalnames.Push("WeightedPressure");
    scalnames.Push("DivOfBulkVeclocity");
    vecnames.Push("BulkVelocity");
    fFluxPressureAnalysis.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    fFluxPressureAnalysis.PostProcess(div);
    
}

/** @brief Create computational meshes using space odissey */
void TRMOrchestra::CreateCompMeshes(TRMRawData &rawdata){
    
}

/// Transfer the flux solution to the saturation mesh
void TRMOrchestra::TransferToSaturationMesh()
{
    //
}

/** @brief Project an exact solution */
void TRMOrchestra::ProjectExactSolution()
{
    TPZAutoPointer<TPZCompMesh> mesh = fSpaceGenerator.GetFluxCmesh();
    if (!mesh) {
        DebugStop();
    }
    // copiar os materiais
    std::map<int, TPZMaterial *> matmap;// = mesh->MaterialVec();
    std::map<int, TPZMaterial *>::iterator it;
    for (it= mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        it->second->Clone(matmap);
    }
    
    // put L2 projection material
    TPZVecL2 *vecmat = new TPZVecL2(_ReservMatId);
    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(TRMOrchestra::ExactFlux);
    vecmat->SetForcingFunction(force);
    for (it = mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        delete it->second;
    }
    mesh->MaterialVec().clear();
    mesh->InsertMaterialObject(vecmat);
    {
        TPZAnalysis an(mesh,false);
        TPZSkylineStructMatrix strmat(mesh);
        std::set<int> matids;
        matids.insert(_ReservMatId);
        strmat.SetMaterialIds(matids);
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Run();
        TPZStack<std::string> scalnames,vecnames;
        vecnames.Push("state");
        an.DefineGraphMesh(3, scalnames, vecnames, "../ProjectFlux.vtk");
        an.PostProcess(0);
    }
    delete vecmat;
    mesh->MaterialVec() = matmap;
    
    mesh = fSpaceGenerator.GetPressureMesh();
    matmap.clear();// = mesh->MaterialVec();
    for (it= mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        it->second->Clone(matmap);
    }
    int nstate = 1;
    TPZManVector<STATE,1> sol(1,1.);
    int dim = 3;
    TPZL2Projection *l2proj = new TPZL2Projection(_ReservMatId,dim,nstate,sol);
    TPZAutoPointer<TPZFunction<STATE> > force2 = new TPZDummyFunction<STATE>(TRMOrchestra::ExactPressure);
    l2proj->SetForcingFunction(force2);

    for (it = mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        delete it->second;
    }
    mesh->MaterialVec().clear();
    mesh->InsertMaterialObject(l2proj);
    {
        TPZAnalysis an(mesh,false);
        TPZSkylineStructMatrix strmat(mesh);
        std::set<int> matids;
        matids.insert(_ReservMatId);
        strmat.SetMaterialIds(matids);
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Run();
        TPZStack<std::string> scalnames,vecnames;
        scalnames.Push("state");
        an.DefineGraphMesh(3, scalnames, vecnames, "../ProjectPressure.vtk");
        an.PostProcess(0);
    }
    delete l2proj;
    mesh->MaterialVec() = matmap;
    
    mesh = fSpaceGenerator.GetMixedCmesh();
    for (it = mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        TPZBndCond * bnd = dynamic_cast<TPZBndCond *>(it->second);
        if (bnd)
        {
            bnd->SetForcingFunction(0,force2);
        }
    }

}

/** @brief exact pressure */
void TRMOrchestra::ExactPressure(const TPZVec<REAL> &x, TPZVec<STATE> &pressure)
{
    pressure[0] = x[2];
}

/** @brief exact pressure */
void TRMOrchestra::ExactFlux(const TPZVec<REAL> &x, TPZVec<STATE> &flux)
{
    flux[0] = 0.;
    flux[1] = 0.;
    flux[2] = -1.;
    
}

