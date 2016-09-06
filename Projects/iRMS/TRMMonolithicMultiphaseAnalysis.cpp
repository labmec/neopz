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
    fmeshvec.Resize(3);
    ferror = 1.0;
    fdx_norm = 1.0;
    fk_iterations = 0;
    
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

    this->Assemble(); // This is stupid! i cannot easily compute rhs without "fake" ek
    fR_n = this->Rhs();
    fR_n += fR; // total residue
    ferror =  Norm(fR_n); // residue error

    
}

void TRMMonolithicMultiphaseAnalysis::QuasiNewtonIteration(){
    
   
    if (k_ietrarions() == 1) {
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
    
    this->Mesh()->LoadSolution(fX_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    
    this->Assemble(); // This is stupid! i cannot easily compute rhs without "fake" ek
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
    
    this->SimulationData()->SetCurrentStateQ(true);
    this->LoadSolution(fX_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());

    ferror = 1.0;
    
    STATE dt_min    = fSimulationData->dt_min();
    STATE dt_max    = fSimulationData->dt_max();
    STATE dt_up     = fSimulationData->dt_up();
    STATE dt_down   = fSimulationData->dt_down();
    STATE dt        = fSimulationData->dt();
    
    STATE epsilon_res = this->SimulationData()->epsilon_res();
    STATE epsilon_cor = this->SimulationData()->epsilon_cor();
    int n  =   this->SimulationData()->n_corrections();
    
    for (int k = 1; k <= n; k++) {

        this->Set_k_ietrarions(k);
        
        if (fSimulationData->IsQuasiNewtonQ()) {
            this->QuasiNewtonIteration();
        }
        else{
            this->NewtonIteration();
        }

#ifdef PZDEBUG
//        this->Solver().Matrix()->Print("K = ", std::cout,EMathematicaInput);
//        fR.Print("R = ", std::cout,EMathematicaInput);
//        fX.Print("X = ", std::cout,EMathematicaInput);        
//        fR_n.Print("Rn = ", std::cout,EMathematicaInput);
//        fX_n.Print("Xn = ", std::cout,EMathematicaInput);
#endif
        
        if(ferror < epsilon_res || fdx_norm < epsilon_cor)
        {
            std::cout << "Monolithic:: Converged with iterations:  " << k << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
            if (k == 1 && dt_max > dt && dt_up > 1.0) {
                dt *= dt_up;
                if(dt_max < dt ){
                    fSimulationData->Setdt(dt_max);
                }
                else{
                    fSimulationData->Setdt(dt);
                }
                std::cout << "Monolithic:: Increasing time step to " << fSimulationData->dt()/86400.0 << "; (day): " << std::endl;
            }
            
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
            std::cout << "Monolithic:: Decreasing time step to " << fSimulationData->dt()/86400.0 << "; (day): " << std::endl;
            std::cout << "Monolithic:: Restarting current time step correction " << std::endl;
            
            this->SimulationData()->SetCurrentStateQ(false);
            this->LoadSolution(fX);
            
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
            this->AssembleResidual();
            fR = this->Rhs();
            
            this->SimulationData()->SetCurrentStateQ(true);
            fX_n = fX;
            k = 1;
        }

        
    }
    
    std::cout << "Warning:: Exit max iterations with min dt:  " << fSimulationData->dt()/86400.0 << "; (day) " << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
    
    
}

void TRMMonolithicMultiphaseAnalysis::PostProcessStep(){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    int dim = this->Mesh()->Reference()->Dimension();
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
        scalnames.Push("div_q");
        vecnames.Push("u");
        vecnames.Push("q");
    }

    if (fSimulationData->IsTwoPhaseQ()) {
        scalnames.Push("p");
        scalnames.Push("div_u");
        scalnames.Push("div_q");
        scalnames.Push("s_a");
        scalnames.Push("s_b");        
        vecnames.Push("u");
        vecnames.Push("q");
    }
    
    if (fSimulationData->IsThreePhaseQ()) {
        scalnames.Push("p");
        scalnames.Push("div_u");
        scalnames.Push("div_q");
        scalnames.Push("s_a");
        scalnames.Push("s_b");
        scalnames.Push("s_c");        
        vecnames.Push("u");
        vecnames.Push("q");
    }

    this->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    this->PostProcess(div);
    
}