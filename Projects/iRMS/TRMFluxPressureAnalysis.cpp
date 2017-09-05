//
//  TRMFluxPressureAnalysis.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMFluxPressureAnalysis.h"
#include "pzcheckmesh.h"


TRMFluxPressureAnalysis::TRMFluxPressureAnalysis() : TPZAnalysis() {
    
    /** @brief define the simulation data */
    fSimulationData = NULL;
    
    /** @brief define the transfer matrices */
    fTransfer = NULL;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    fmeshvec.Resize(2);
    
    /** @brief Part of residue at n state  */
    fR_n.Resize(0,0);
    
    /** @brief Part of residue at last state  */
    fR.Resize(0,0);
    
    /** @brief Solution ate n state */
    fX_n.Resize(0,0);
    
    /** @brief Solution at past state */
    fX.Resize(0, 0);
    
    /** @brief Residue error */
    ferror    = 0.0;
    
    /** @brief Correction variation */
    fdx_norm  = 0.0;
    
    /** @brief number of newton corrections */
    fk_iterations = 0;
    
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
    fR              = copy.fR;
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
        fR              = other.fR;
        fX_n            = other.fX_n;
        fX              = other.fX;
        ferror          = other.ferror;
        fdx_norm        = other.fdx_norm;
    }
    return *this;
}

/** @brief Resize and fill residue and solution vectors */
void TRMFluxPressureAnalysis::AdjustVectors(){
    
    if(fSolution.Rows() == 0 /* || fRhs.Rows() == 0 */){
        DebugStop();
    }
    
    fX.Resize(fSolution.Rows(),1);
    fX.Zero();
    fX_n.Resize(fSolution.Rows(),1);
    fX_n.Zero();
    fR_n.Resize(fSolution.Rows(),1);
    fR_n.Zero();
    fR.Resize(fSolution.Rows(),1);
    fR.Zero();
}

void TRMFluxPressureAnalysis::NewtonIteration(){
    
    this->Assemble();
    this->Rhs() += fR; // total residue
    this->Rhs() *= -1.0;
    
    this->Solve(); // update correction
    fdx_norm = Norm(this->Solution()); // correction variation
    
    fX_n += this->Solution(); // update state

    this->UpdateMemory_at_n();
    
    this->AssembleResidual();
    fR_n = this->Rhs();
    fR_n += fR; // total residue
    ferror =  Norm(fR_n); // residue error
    
}

void TRMFluxPressureAnalysis::QuasiNewtonIteration(){
    
    if (k_ietrarions() <= 4) {
        this->Assemble();
    }
    else{
        this->AssembleResidual();
    }
    
    this->Rhs() += fR; // total residue
    this->Rhs() *= -1.0;
    
    this->Solve(); // update correction

    fdx_norm = Norm(this->Solution()); // correction variation
    
    fX_n += this->Solution(); // update state
    
    this->UpdateMemory_at_n();
    
    this->AssembleResidual();
    fR_n = this->Rhs();
    fR_n += fR; // total residue
    ferror =  Norm(fR_n); // residue error
    
}

void TRMFluxPressureAnalysis::ExcecuteOneStep(){
    
    this->SimulationData()->SetCurrentStateQ(false);
    this->LoadSolution(fX);
    this->UpdateMemory();
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    
    this->AssembleResidual();
    fR = this->Rhs();
    
    this->SimulationData()->SetCurrentStateQ(true);
    this->LoadSolution(fX_n);
    this->UpdateMemory_at_n();
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());

    this->AssembleResidual();
    fR_n = this->Rhs();
    
    ferror = Norm(fR+fR_n);
    this->Set_k_ietrarions(0);
    
    STATE epsilon_res = this->SimulationData()->epsilon_res();
    STATE epsilon_cor = this->SimulationData()->epsilon_cor();
    int n  =   this->SimulationData()->n_corrections();
    
    // check if the new state is almost stationary
    if(ferror < epsilon_res){
        return;
    }
    
    for (int k = 1; k <= n; k++) {
        
        this->Set_k_ietrarions(k);
        
        if (fSimulationData->IsQuasiNewtonQ()) {
            this->QuasiNewtonIteration();
        }
        else{
            this->NewtonIteration();
        }
        

#ifdef PZDEBUG
//        fR.Print("R = ", std::cout,EMathematicaInput);
//        fX.Print("X = ", std::cout,EMathematicaInput);
//        fR_n.Print("Rn = ", std::cout,EMathematicaInput);
//        fX_n.Print("Xn = ", std::cout,EMathematicaInput);
#endif
        
        if(ferror < epsilon_res || fdx_norm < epsilon_cor)
        {
            std::cout << "Parabolic:: Converged with iterations:  " << k << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
            return;
        }
    
    }
    
    std::cout << "Parabolic:: Exit max iterations with min dt:  " << fSimulationData->dt()/86400.0 << "; (day) " << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
    
    
}

/** @brief Update memory using the Transfer object at state n */
void TRMFluxPressureAnalysis::UpdateMemory_at_n(){

    Mesh()->LoadSolution(fX_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, Mesh());
    
    // Volumetric update    
    fTransfer->u_To_Mixed_Memory(fmeshvec[0], Mesh());
    fTransfer->p_To_Mixed_Memory(fmeshvec[1], Mesh());
    fTransfer->p_avg_Memory_Transfer(Mesh());
    
}

/** @brief Update memory using the Transfer object */
void TRMFluxPressureAnalysis::UpdateMemory(){
    
    Mesh()->LoadSolution(fX);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, Mesh());
    
    // Volumetric update
    fTransfer->u_To_Mixed_Memory(fmeshvec[0], Mesh());
    fTransfer->p_To_Mixed_Memory(fmeshvec[1], Mesh());
    fTransfer->p_avg_Memory_Transfer(Mesh());
    
}

void TRMFluxPressureAnalysis::PostProcessStep(){
    
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    const int dim = this->Mesh()->Dimension();
    int div = 0;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile;
    if (fSimulationData->IsInitialStateQ()) {
        plotfile =  "DualSegregatedDarcyOnBox_I.vtk";
        return;
    }
    else{
        plotfile =  "DualSegregatedDarcyOnBox.vtk";
    }
    
    scalnames.Push("p");
    scalnames.Push("div_u");
    scalnames.Push("cfl");
    scalnames.Push("phi");
    vecnames.Push("u");
    vecnames.Push("kappa");
    
    this->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    this->PostProcess(div);
    
}
