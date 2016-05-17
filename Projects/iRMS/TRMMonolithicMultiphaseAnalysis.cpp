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

// set the type of the system
//void TRMMonolithicMultiphaseAnalysis::SetFluidData(TPZVec< TPZAutoPointer<Phase> > PVTData){
//    
//    TPZStack<std::string> System =  fSimulationData->GetsystemType();
//    int nphases = System.size();
//    
//    if (fSimulationData->IsOnePhaseQ()) {
//        
//        for(int iphase = 0; iphase < nphases; iphase++){
//            
//            if (!strcmp("Water", System[iphase].c_str())){
//                falpha_fluid = PVTData[0];
//                fbeta_fluid = PVTData[1];
//                fgamma_fluid = PVTData[2];
//            }
//            
//            if (!strcmp("Oil", System[iphase].c_str())){
//                falpha_fluid = PVTData[1];
//                fbeta_fluid  = PVTData[0];
//                fgamma_fluid = PVTData[2];
//            }
//            
//            if (!strcmp("Gas", System[iphase].c_str())){
//                falpha_fluid = PVTData[2];
//                fbeta_fluid = PVTData[0];
//                fgamma_fluid = PVTData[1];
//            }
//            
//        }
//        
//    }
//    
//    if(fSimulationData->IsTwoPhaseQ()){
//        
//        fmeshvecini.Resize(3);
//        fmeshvec.Resize(3);
//        
//        for(int iphase = 0; iphase < nphases; iphase++){
//            
//            switch (iphase) {
//                case 0:
//                {
//                    if (!strcmp("Water", System[iphase].c_str())){
//                        falpha_fluid = PVTData[0];
//                    }
//                    
//                    if (!strcmp("Oil", System[iphase].c_str())){
//                        falpha_fluid = PVTData[1];
//                    }
//                    
//                    if (!strcmp("Gas", System[iphase].c_str())){
//                        falpha_fluid = PVTData[2];
//                    }
//                    
//                    fgamma_fluid = PVTData[2];
//                    
//                }
//                    break;
//                    
//                case 1:
//                {
//                    if (!strcmp("Water", System[iphase].c_str())){
//                        fbeta_fluid = PVTData[0];
//                    }
//                    
//                    if (!strcmp("Oil", System[iphase].c_str())){
//                        fbeta_fluid = PVTData[1];
//                    }
//                    
//                    if (!strcmp("Gas", System[iphase].c_str())){
//                        fbeta_fluid = PVTData[2];
//                    }
//                    
//                }
//                    fgamma_fluid = PVTData[2];
//                    break;
//                default:
//                {
//                    DebugStop();
//                }
//                    break;
//            }
//            
//        }
//        
//    }
//    
//    if(fSimulationData->IsThreePhaseQ()){
//        
//        fmeshvecini.Resize(4);
//        fmeshvec.Resize(4);
//        
//        std::cout << "System not impelmented " << System << std::endl;
//        DebugStop();
//    }
//    
//}


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

void TRMMonolithicMultiphaseAnalysis::PostProcessStep(){
    
    const int dim = 3;
    int div = 0;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "DualMonolithicDarcyOnBox.vtk";
    scalnames.Push("p");
    scalnames.Push("div_u");
    vecnames.Push("u");
    this->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    this->PostProcess(div);
    
}