//
//  TRMTransportAnalysis.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMTransportAnalysis.h"


TRMTransportAnalysis::TRMTransportAnalysis() : TPZAnalysis() {
    
    /** @brief define the simulation data */
    fSimulationData = NULL;
    
    /** @brief define the transfer matrices */
    fTransfer = NULL;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    fmeshvec.Resize(1);
    
    /** @brief Part of residue at n state  */
    fR_n.Resize(0,0);
    
    /** @brief Part of residue at past state  */
    fR.Resize(0,0);
    
    /** @brief Solution ate n state */
    fX_n.Resize(0,0);
    
    /** @brief Solution at past state */
    fX.Resize(0, 0);
    
    /** @brief Residue error */
    ferror    = 1.0;
    
    /** @brief Correction variation */
    fdx_norm  = 1.0;
    
    /** @brief number of newton corrections */
    fk_iterations = 0;
    
}

TRMTransportAnalysis::~TRMTransportAnalysis(){
    
}

/** @brief Copy constructor $ */
TRMTransportAnalysis::TRMTransportAnalysis(const TRMTransportAnalysis &copy)
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

/** @brief Assignemnt operator $ */
TRMTransportAnalysis & TRMTransportAnalysis::operator=(const TRMTransportAnalysis &other)
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
void TRMTransportAnalysis::AdjustVectors(){
    
    if(fSolution.Rows() == 0  /*|| fRhs.Rows() == 0 */ ){
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
    fR.Resize(fSolution.Rows(),1);
    fR.Zero();
}

void TRMTransportAnalysis::NewtonIteration(){
    
    this->Assemble();
    this->Rhs() += fR; // total residue
    this->Rhs() *= -1.0;
    
    this->Solve(); // update correction
    fdx_norm = Norm(this->Solution()); // correction variation
    
    fX_n += this->Solution(); // update state
    
//    this->Rhs().Print("rhs = ");
//    fX_n.Print("x_n = ");
    
    this->Mesh()->LoadSolution(fX_n);
    this->UpdateMemory_at_n();
    
    this->Assemble();
    fR_n = this->Rhs();
    fR_n += fR; // total residue
    ferror =  Norm(fR_n); // residue error
    
}

void TRMTransportAnalysis::ExcecuteOneStep(){
    
    
    this->SimulationData()->SetCurrentStateQ(false);
    this->LoadSolution(fX);
    this->UpdateMemory();
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    this->AssembleResidual();
    fR = this->Rhs();
    
    this->SimulationData()->SetCurrentStateQ(true);
    this->LoadSolution(fX_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    this->UpdateMemory_at_n();
    
    ferror = 1.0;
    
    STATE dt_min    = fSimulationData->dt_min();
//    STATE dt_max    = fSimulationData->dt_max();
//    STATE dt_up     = fSimulationData->dt_up();
    STATE dt_down   = fSimulationData->dt_down();
    STATE dt        = fSimulationData->dt();
    
    STATE epsilon_res = this->SimulationData()->epsilon_res();
    STATE epsilon_cor = this->SimulationData()->epsilon_cor();
    int n  =   this->SimulationData()->n_corrections();
    
    for (int k = 1; k <= n; k++) {
        
        this->fk_iterations = k;
        this->NewtonIteration();
        
//#ifdef PZDEBUG
//        fR.Print("R = ", std::cout,EMathematicaInput);
//        fX.Print("X = ", std::cout,EMathematicaInput);
//        fR_n.Print("Rn = ", std::cout,EMathematicaInput);
//        fX_n.Print("Xn = ", std::cout,EMathematicaInput);
//        this->Solver().Matrix()->Print("K = ",std::cout,EMathematicaInput);
//#endif
        
        if(ferror < epsilon_res || fdx_norm < epsilon_cor)
        {
            std::cout << "Hyperbolic:: Converged with iterations:  " << k << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
//            if (k == 1 && dt_max > dt && dt_up > 1.0) {
//                dt *= dt_up;
//                if(dt_max < dt ){
//                    fSimulationData->Setdt(dt_max);
//                }
//                else{
//                    fSimulationData->Setdt(dt);
//                }
//                std::cout << "Hyperbolic:: Increasing time step to " << fSimulationData->dt()/86400.0 << "; (day): " << std::endl;
//            }
            
            fX = fX_n;
            return;
        }
        
        if(k == n  && dt > dt_min && dt_down < 1.0){
            dt *= dt_down;
            if(dt_min > dt ){
                fSimulationData->Setdt(dt_min);
            }
            else{
                fSimulationData->Setdt(dt);
            }
            std::cout << "Hyperbolic:: Decreasing time step to " << fSimulationData->dt()/86400.0 << "; (day): " << std::endl;
            std::cout << "Hyperbolic:: Restarting current time step correction " << std::endl;
            
            this->SimulationData()->SetCurrentStateQ(false);
            this->LoadSolution(fX);
            
            this->UpdateMemory();
            this->AssembleResidual();
            fR = this->Rhs();
            
            this->SimulationData()->SetCurrentStateQ(true);
            this->UpdateMemory_at_n();
            fX_n = fX;
            k = 1;
        }
        
        
    }
    
    std::cout << "Hyperbolic:: Exit max iterations with min dt:  " << fSimulationData->dt()/86400.0 << "; (day) " << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
    
    
}

/** @brief Update memory using the Transfer object at state n */
void TRMTransportAnalysis::UpdateMemory_at_n(){
    
    Mesh()->LoadSolution(fX_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, Mesh());
    
    // Volumetric update
    if (fSimulationData->IsTwoPhaseQ()) {
        fTransfer->s_To_Transport_Memory(fmeshvec[0], Mesh(),0);
    }
    
    // Volumetric update
    if (fSimulationData->IsThreePhaseQ()) {
        fTransfer->s_To_Transport_Memory(fmeshvec[0], Mesh(),0);
        fTransfer->s_To_Transport_Memory(fmeshvec[1], Mesh(),1);
    }
    
}

/** @brief Update memory using the Transfer object */
void TRMTransportAnalysis::UpdateMemory(){
    
    Mesh()->LoadSolution(fX);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, Mesh());
    
    // Volumetric update
    if (fSimulationData->IsTwoPhaseQ()) {
        fTransfer->s_To_Transport_Memory(fmeshvec[0], Mesh(),0);
    }

    // Volumetric update
    if (fSimulationData->IsThreePhaseQ()) {
        fTransfer->s_To_Transport_Memory(fmeshvec[0], Mesh(),0);        
        fTransfer->s_To_Transport_Memory(fmeshvec[1], Mesh(),1);
    }
    
}

void TRMTransportAnalysis::PostProcessStep(){
    
    const int dim = this->Mesh()->Reference()->Dimension();
    int div = 0;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "DualSegregatedDarcyOnBox_Saturations.vtk";
    scalnames.Push("s_a");
    scalnames.Push("s_b");
    
    if (fSimulationData->IsThreePhaseQ()) {
        scalnames.Push("s_c");
    }
    
    this->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    this->PostProcess(div);
    
}