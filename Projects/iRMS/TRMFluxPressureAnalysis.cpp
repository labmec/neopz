//
//  TRMFluxPressureAnalysis.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMFluxPressureAnalysis.h"


TRMFluxPressureAnalysis::TRMFluxPressureAnalysis() : TPZAnalysis() {
    
    fSimulationData = NULL;
    fmeshvec.Resize(2); // Start with monophasic approach
    ferror = 1.0;
    fdx_norm = 1.0;
    
}

TRMFluxPressureAnalysis::~TRMFluxPressureAnalysis(){
    
}

/** @brief Resize and fill residue and solution vectors */
void TRMFluxPressureAnalysis::AdjustVectors(){
    
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
    fR_n.Resize(fSolution.Rows(),1);
    fR_n.Zero();
}



void TRMFluxPressureAnalysis::NewtonIteration(){
    
    this->Assemble();
    this->Rhs() += fR_n; // total residue
    this->Rhs() *= -1.0;
    
    this->Solve(); // update correction
    fdx_norm = Norm(this->Solution()); // correction variation
    
    fX_n += this->Solution(); // update state
    
    this->Mesh()->LoadSolution(fX_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    this->AssembleResidual();
    fR_n += this->Rhs(); // total residue
    ferror =  Norm(fR_n); // residue error
    
    
}

void TRMFluxPressureAnalysis::ExcecuteOneStep(){
    
    //    this->SimulationData()->SetCurrentStateQ(false);
    //    this->LoadSolution(fX);
    //
    //    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    //    this->AssembleResidual();
    //    fR = this->Rhs();
    
    //    this->SimulationData()->SetCurrentStateQ(true);
    //    this->LoadSolution(fX_n);
    //    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    //    this->AssembleResidual();
    //    fR_n = this->Rhs();
    //
    //    fR_n += fR;
    
    ferror = 1.0;
    
    
    STATE epsilon_res = this->SimulationData()->epsilon_res();
    STATE epsilon_cor = this->SimulationData()->epsilon_cor();
    int n  =   this->SimulationData()->n_corrections();
    
    
    
    for (int k = 1; k <= n; k++) {
        
        this->NewtonIteration();
        
        //#ifdef PZDEBUG
        //        fR.Print("R = ", std::cout,EMathematicaInput);
        //        fR_n.Print("Rn = ", std::cout,EMathematicaInput);
        //        fX_n.Print("X = ", std::cout,EMathematicaInput);
        //#endif
        
        
        if(ferror < epsilon_res || fdx_norm < epsilon_cor)
        {
            std::cout << "Converged with iterations:  " << k << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
            fX = fX_n;
            return;
        }
        
    }
    
    std::cout << "Exit with iterations:  " << n << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
    
    
}

/** @brief Update memory using the Transfer object at state n */
void TRMFluxPressureAnalysis::UpdateMemory_at_n(){

    Mesh()->LoadSolution(fX_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, Mesh());
    
    // Volumetric update    
    fTransfer->Transfer_u_To_Mixed_Memory(fmeshvec[0], Mesh());
    fTransfer->Transfer_p_To_Mixed_Memory(fmeshvec[1], Mesh());
    
}

/** @brief Update memory using the Transfer object */
void TRMFluxPressureAnalysis::UpdateMemory(){
    
    Mesh()->LoadSolution(fX);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, Mesh());
    
    // Volumetric update
    fTransfer->Transfer_u_To_Mixed_Memory(fmeshvec[0], Mesh());
    fTransfer->Transfer_p_To_Mixed_Memory(fmeshvec[1], Mesh());
    
}

void TRMFluxPressureAnalysis::PostProcessStep(){
    
    const int dim = 3;
    int div = 0;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "DualSegregatedDarcyOnBox.vtk";
    scalnames.Push("p");
    scalnames.Push("div_u");
    vecnames.Push("u");
    this->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    this->PostProcess(div);
    
}