//
//  TRMMonolithicMultiphaseAnalysis.cpp
//  PZ
//
//  Created by Omar on 5/11/16.
//
//

#include "TRMMonolithicMultiphaseAnalysis.h"


TRMMonolithicMultiphaseAnalysis::TRMMonolithicMultiphaseAnalysis() : TPZAnalysis() {
    
    fSimulationData = NULL;
    fmeshvec.Resize(2); // Start with monophasic approach
    ferror = 1.0;
    fdx_norm = 1.0;
    
}

TRMMonolithicMultiphaseAnalysis::~TRMMonolithicMultiphaseAnalysis(){
    
}

/** @brief Resize and fill residue and solution vectors */
void TRMMonolithicMultiphaseAnalysis::AdjustVectors(){
    
    if(fSolution.Rows() == 0){
        DebugStop();
    }
    
    TPZBuildMultiphysicsMesh::AddElements(fmeshvec, this->Mesh());
    TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, this->Mesh());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, this->Mesh());
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());    
    
    fX.Resize(fSolution.Rows(),1);
    fX.Zero();
    fX_n.Resize(fSolution.Rows(),1);
    fX_n.Zero();
    fR.Resize(fSolution.Rows(),1);
    fR.Zero();
    fR_n.Resize(fSolution.Rows(),1);
    fR_n.Zero();
}



void TRMMonolithicMultiphaseAnalysis::NewtonIteration(){
 
    this->Assemble();
    this->Rhs() += fR; // total residue
    this->Rhs() *= -1.0;
    
    this->Solve(); // update correction
    fdx_norm = Norm(this->Solution()); // correction variation
    
    fX_n += this->Solution(); // update state
    
    this->Mesh()->LoadSolution(fX_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());

    this->AssembleResidual();
    fR_n = this->Rhs();
    fR_n += fR; // total residue
    ferror =  Norm(fR_n); // residue error

    
}

void TRMMonolithicMultiphaseAnalysis::ExcecuteOneStep(){
       
    this->SimulationData()->SetCurrentStateQ(false);
    this->LoadSolution(fX);
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    this->AssembleResidual();
    fR = this->Rhs();
    
    
    STATE v = -0.25;
    fX_n(0,0) = v;
    fX_n(1,0) = v;
    fX_n(2,0) = v;
    fX_n(3,0) = v;

    fX_n(4,0) = -v;
    fX_n(5,0) = -v;
    fX_n(6,0) = -v;
    fX_n(7,0) = -v;
    
    STATE  pr = 1.0e+7;
    STATE  pl = 1.1e+7;
    fX_n(36,0) = pl;
    fX_n(37,0) = pr;
    fX_n(38,0) = pr;
    fX_n(39,0) = pl;
    fX_n(40,0) = pl;
    fX_n(41,0) = pr;
    fX_n(42,0) = pr;
    fX_n(43,0) = pl;
    fX_n(44,0) = 1.0;
    
    this->SimulationData()->SetCurrentStateQ(true);
    this->LoadSolution(fX_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());

    ferror = 1.0;
    
    STATE epsilon_res = this->SimulationData()->epsilon_res();
    STATE epsilon_cor = this->SimulationData()->epsilon_cor();
    int n  =   this->SimulationData()->n_corrections();


    
    for (int k = 1; k <= n; k++) {

        this->Assemble();
        fR_n = this->Rhs();
//        this->NewtonIteration();
        
#ifdef PZDEBUG
        fR.Print("R = ", std::cout,EMathematicaInput);
        fX.Print("X = ", std::cout,EMathematicaInput);        
        fR_n.Print("Rn = ", std::cout,EMathematicaInput);
        fX_n.Print("Xn = ", std::cout,EMathematicaInput);
#endif
        
        if(ferror < epsilon_res || fdx_norm < epsilon_cor)
        {
            std::cout << "Converged with iterations:  " << k << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
            fX = fX_n;
            return;
        }
        
    }
    
    std::cout << "Warning:: Exit with iterations:  " << n << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
    
    
}

void TRMMonolithicMultiphaseAnalysis::PostProcessStep(){
    
    const int dim = 3;
    int div = 0;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile;
    if (fSimulationData->IsInitialStateQ()) {
        plotfile =  "DualMonolithicDarcyOnBox_I.vtk";
    }
    else{
        plotfile =  "DualMonolithicDarcyOnBox.vtk";
    }

    if (fSimulationData->IsOnePhaseQ()) {
        scalnames.Push("p");
        scalnames.Push("div_u");
        vecnames.Push("u");
    }

    if (fSimulationData->IsTwoPhaseQ()) {
        scalnames.Push("p");
        scalnames.Push("s_a");        
        scalnames.Push("div_u");
        vecnames.Push("u");
    }

    this->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    this->PostProcess(div);
    
}