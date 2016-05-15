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
    
    this->Rhs() = -1.0*fResidue_n;
    this->Assemble();
    this->Solve();
    fdx_norm = Norm(this->Solution()); // correction variation
    
    fSolution_n += this->Solution(); // update
    
    this->Mesh()->LoadSolution(fSolution_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    this->AssembleResidual();
    
    fResidue_n  = this->Rhs() + fResidue;
    ferror =  Norm(fResidue_n); // residue error

    
}

void TRMMonolithicMultiphaseAnalysis::ExcecuteOneStep(){
    
    this->SimulationData()->SetCurrentStateQ(false);
    this->LoadSolution(fSolution);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    this->AssembleResidual();
    fResidue = this->Rhs();
    
    this->SimulationData()->SetCurrentStateQ(true);
    this->LoadSolution(fSolution);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());    
    this->AssembleResidual();
    fResidue_n = this->Rhs();
    
    fResidue_n += fResidue;
    ferror = 1.0;
    
    
    STATE epsilon_res = this->SimulationData()->epsilon_res();
    STATE epsilon_cor = this->SimulationData()->epsilon_cor();
    int n  =   this->SimulationData()->n_corrections();


    
    for (int k = 0; k < n; k++) {
        this->NewtonIteration();
        
#ifdef PZDEBUG
                std::stringstream sout;
                fResidue.Print("R = ", sout,EMathematicaInput);
                fResidue_n.Print("Rn = ", sout,EMathematicaInput);
                fSolution_n.Print("X = ", sout,EMathematicaInput);
                std::cout << sout << std::endl;
#endif
        
        if(ferror < epsilon_res || fdx_norm < epsilon_cor)
        {
            std::cout << "Converged with iterations:  " << k << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
            fSolution = fSolution_n;
            return;
        }
        
    }
    
    std::cout << "Exit with iterations:  " << n << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
    
    
}

void TRMMonolithicMultiphaseAnalysis::PostProcessStep(){
    
    const int dim = 3;
    int div = 1;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "DualMonolithicDarcyOnBox.vtk";
    scalnames.Push("p");
    scalnames.Push("div_u");
    vecnames.Push("u");
    this->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    this->PostProcess(div);
    
}