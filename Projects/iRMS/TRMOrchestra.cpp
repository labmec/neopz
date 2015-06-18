//
//  TRMOrchestra.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMOrchestra.h"


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
void TRMOrchestra::CreateAnalPrimal(TRMSpaceOdissey spacegenerator){
    
    TPZManVector<int,2> dx(2,10), dy(2,1), dz(2,1);
    dx[0] = 100;
    dy[0] = 100;
    dz[0] = 100;
    
    spacegenerator.CreateGeometricBoxMesh(dx, dy, dz);
//    spacegenerator.CreateGeometricReservoirMesh();
    spacegenerator.PrintGeometry();
    fgmesh = spacegenerator.GetGmesh();
    spacegenerator.CreateH1Cmesh();
    
    TPZAutoPointer<TPZCompMesh > Cmesh = spacegenerator.GetH1Cmesh();
    
    // Analysis
    bool mustOptimizeBandwidth = true;
    TPZAnalysis * AnalPrimal = new TPZAnalysis(Cmesh.operator->(),mustOptimizeBandwidth);
    int numofThreads = 8;
    
    TPZSkylineNSymStructMatrix skylnsym(Cmesh.operator->());
    TPZStepSolver<STATE> step;
    skylnsym.SetNumThreads(numofThreads);
    step.SetDirect(ELU);
    AnalPrimal->SetStructuralMatrix(skylnsym);
    AnalPrimal->SetSolver(step);;
    AnalPrimal->Run();
    
    const int dim = 3;
    int div = 2;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "PrimalDarcy.vtk";
    scalnames.Push("Pressure");
    vecnames.Push("MinusKGradU");
    AnalPrimal->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    AnalPrimal->PostProcess(div);
    
}

/** @brief Create a dual analysis using space odissey */
void TRMOrchestra::CreateAnalDual(TRMSpaceOdissey spacegenerator){
    
    TPZManVector<int,2> dx(2,10), dy(2,1), dz(2,1);
    dx[0] = 100;
    dy[0] = 100;
    dz[0] = 100;
    
    spacegenerator.CreateGeometricBoxMesh(dx, dy, dz);
    spacegenerator.PrintGeometry();
    fgmesh = spacegenerator.GetGmesh();
    spacegenerator.CreateMixedCmesh();
    
    TPZAutoPointer<TPZCompMesh > Cmesh = spacegenerator.GetMixedCmesh();
    
    // Analysis
    bool mustOptimizeBandwidth = true;
    TPZAnalysis * AnalDual = new TPZAnalysis(Cmesh.operator->(),mustOptimizeBandwidth);
    int numofThreads = 8;
    
    TPZSkylineNSymStructMatrix skylnsym(Cmesh.operator->());
    TPZStepSolver<STATE> step;
    skylnsym.SetNumThreads(numofThreads);
    step.SetDirect(ELU);
    AnalDual->SetStructuralMatrix(skylnsym);
    AnalDual->SetSolver(step);
    AnalDual->Assemble();
    AnalDual->Rhs() *= -1.0;
    AnalDual->Solve();
    
    const int dim = 3;
    int div = 2;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "DualDarcy.vtk";
    scalnames.Push("WeightedPressure");
    scalnames.Push("DivOfBulkVeclocity");
    vecnames.Push("BulkVelocity");
    AnalDual->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    AnalDual->PostProcess(div);
    
}

/** @brief Create computational meshes using space odissey */
void TRMOrchestra::CreateCompMeshes(TRMRawData &rawdata){
    
}
