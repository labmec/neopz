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
void TRMOrchestra::CreateAnalysisPrimal(TRMSpaceOdissey spacegenerator){
    
    TPZManVector<int,2> dx(2,1), dy(2,1), dz(2,1);
    dx[0] = 50.0;
    dy[0] = 50.0;
    dz[0] = 50.0;
    
//    spacegenerator.CreateGeometricBoxMesh(dx, dy, dz);
    spacegenerator.CreateGeometricReservoirMesh();
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
void TRMOrchestra::CreateAnalysisDual(TRMSpaceOdissey spacegenerator){
    
    TPZManVector<int,2> dx(2,1), dy(2,1), dz(2,1);
    dx[0] = 50.0;
    dy[0] = 50.0;
    dz[0] = 50.0;
    
//    spacegenerator.CreateGeometricBoxMesh(dx, dy, dz);
    spacegenerator.CreateGeometricReservoirMesh();
    spacegenerator.PrintGeometry();
    fgmesh = spacegenerator.GetGmesh();
    spacegenerator.CreateMixedCmesh();
    
    TPZAutoPointer<TPZCompMesh > Cmesh = spacegenerator.GetMixedCmesh();
    
#ifdef DEBUG
    {
        std::ofstream out("../MFCompMesh.txt");
        Cmesh->Print(out);
    }
#endif
    
    // Analysis
    bool mustOptimizeBandwidth = true;
    TPZAnalysis * AnalysisDual = new TPZAnalysis(Cmesh.operator->(),mustOptimizeBandwidth);
    int numofThreads = 0;
    
    TPZSkylineStructMatrix strmat(Cmesh.operator->());
//    TPZSkylineNSymStructMatrix strmat(Cmesh.operator->());
    TPZStepSolver<STATE> step;
    strmat.SetNumThreads(numofThreads);
    step.SetDirect(ELU);
    AnalysisDual->SetStructuralMatrix(strmat);
    AnalysisDual->SetSolver(step);
    AnalysisDual->Assemble();
    AnalysisDual->Rhs() *= -1.0;
    AnalysisDual->Solve();
    std::cout << "Dual dof: " << AnalysisDual->Rhs().Rows() << std::endl;
    
    const int dim = 3;
    int div = 1;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "DualDarcy.vtk";
    scalnames.Push("WeightedPressure");
    scalnames.Push("DivOfBulkVeclocity");
    vecnames.Push("BulkVelocity");
    AnalysisDual->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    AnalysisDual->PostProcess(div);
    
}

/** @brief Create computational meshes using space odissey */
void TRMOrchestra::CreateCompMeshes(TRMRawData &rawdata){
    
}
