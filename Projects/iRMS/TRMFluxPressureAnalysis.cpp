//
//  TRMFluxPressureAnalysis.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMFluxPressureAnalysis.h"


TRMFluxPressureAnalysis::TRMFluxPressureAnalysis() : TPZAnalysis() {
    
    /** @brief define the simulation data */
    fSimulationData = NULL;
    
    /** @brief define the transfer matrices */
    fTransfer = NULL;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    fmeshvec.Resize(2);
    
    /** @brief Part of residue at n state  */
    fR_n.Resize(0,0);
    
    /** @brief Solution ate n state */
    fX_n.Resize(0,0);
    
    /** @brief Solution at past state */
    fX.Resize(0, 0);
    
    /** @brief Residue error */
    ferror    = 1.0;
    
    /** @brief Correction variation */
    fdx_norm  = 1.0;
    
}

TRMFluxPressureAnalysis::~TRMFluxPressureAnalysis(){
    
}

/** @brief Copy constructor $ */
TRMFluxPressureAnalysis::TRMFluxPressureAnalysis(const TRMFluxPressureAnalysis &copy)
{    
    fSimulationData = copy.fSimulationData;
    fTransfer       = copy.fTransfer;
    fmeshvec        = copy.fmeshvec;
    fR_n            = copy.fR_n;
    fX_n            = copy.fX_n;
    fX              = copy.fX;
    ferror          = copy.ferror;
    fdx_norm        = copy.fdx_norm;

}

/** @brief Copy assignemnt operator $ */
TRMFluxPressureAnalysis & TRMFluxPressureAnalysis::operator=(const TRMFluxPressureAnalysis &other)
{
    if (this != & other) {  // prevent self-assignment
        
        fSimulationData = other.fSimulationData;
        fTransfer       = other.fTransfer;
        fmeshvec        = other.fmeshvec;
        fR_n            = other.fR_n;
        fX_n            = other.fX_n;
        fX              = other.fX;
        ferror          = other.ferror;
        fdx_norm        = other.fdx_norm;
    }
    return *this;
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
    
    this->UpdateMemory_at_n();
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    this->AssembleResidual();
    fR_n += this->Rhs();    // total residue
    ferror =  Norm(fR_n);   // residue error
    
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