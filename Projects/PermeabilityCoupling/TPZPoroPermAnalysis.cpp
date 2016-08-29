//
//  TPZPoroPermAnalysis.cpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#include "TPZPoroPermAnalysis.h"


TPZPoroPermAnalysis::TPZPoroPermAnalysis() : TPZAnalysis() {
    
    /** @brief define the simulation data */
    fSimulationData = NULL;
    
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
    ferror    = 1.0;
    
    /** @brief Correction variation */
    fdx_norm  = 1.0;
    
    /** @brief number of newton corrections */
    fk_iterations = 0;
    
}

TPZPoroPermAnalysis::~TPZPoroPermAnalysis(){
    
}

/** @brief Copy constructor $ */
TPZPoroPermAnalysis::TPZPoroPermAnalysis(const TPZPoroPermAnalysis &copy)
{
    fSimulationData = copy.fSimulationData;
    fmeshvec        = copy.fmeshvec;
    fR_n            = copy.fR_n;
    fR              = copy.fR;
    fX_n            = copy.fX_n;
    fX              = copy.fX;
    ferror          = copy.ferror;
    fdx_norm        = copy.fdx_norm;
    
}

/** @brief Copy assignemnt operator $ */
TPZPoroPermAnalysis & TPZPoroPermAnalysis::operator=(const TPZPoroPermAnalysis &other)
{
    if (this != & other) {  // prevent self-assignment
        
        fSimulationData = other.fSimulationData;
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
void TPZPoroPermAnalysis::AdjustVectors(){
    
    if(fSolution.Rows() == 0 /* || fRhs.Rows() == 0 */){
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

void TPZPoroPermAnalysis::QuasiNewtonIteration(){
    
//    fX_n(3,0) = -0.01;
//    fX_n(7,0) = -0.01;
//    fX_n(8,0) = M_PI;
//    fX_n(9,0) = M_PI;
//    fX_n(10,0) = M_PI;
//    fX_n(11,0) = M_PI;
//    
//    this->Update_at_n_State();
//    this->Solution().Print("x = ");
    this->Assemble();
    this->Rhs() += fR; // total residue
    this->Rhs() *= -1.0;
//    this->Rhs().Print("rhs = ");
    this->Solve(); // update correction
    fdx_norm = Norm(this->Solution()); // correction variation
    
    fX_n += this->Solution(); // update state
    
    this->Update_at_n_State();
    
    this->AssembleResidual();
    fR_n = this->Rhs();
    fR_n += fR; // total residue
    ferror =  Norm(fR_n); // residue error
    
}

void TPZPoroPermAnalysis::ExcecuteOneStep(){
    
    this->SimulationData()->SetCurrentStateQ(false);
    this->UpdateState();

    this->AssembleResidual();
    fR = this->Rhs();
    
    this->SimulationData()->SetCurrentStateQ(true);
    this->Update_at_n_State();
    
    ferror = 1.0;
    
    STATE epsilon_res = this->SimulationData()->epsilon_res();
    STATE epsilon_cor = this->SimulationData()->epsilon_cor();
    int n  =   this->SimulationData()->n_corrections();
    
    for (int k = 1; k <= n; k++) {
        
        this->Set_k_ietrarions(k);
        this->QuasiNewtonIteration();

        
//#ifdef PZDEBUG
//        fR.Print("R = ", std::cout,EMathematicaInput);
//        fX.Print("X = ", std::cout,EMathematicaInput);
//        fR_n.Print("Rn = ", std::cout,EMathematicaInput);
//        fX_n.Print("Xn = ", std::cout,EMathematicaInput);
//#endif
        
        if(ferror < epsilon_res || fdx_norm < epsilon_cor)
        {
            std::cout << "PermeabilityCoupling:: Converged with iterations:  " << k << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
            fX = fX_n;
            return;
        }
        
    }
    
    std::cout << "PermeabilityCoupling:: Exit max iterations with min dt:  " << fSimulationData->dt() << "; (secs) " << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
    
    
}

/** @brief update last state solution */
void TPZPoroPermAnalysis::UpdateState(){
    this->LoadSolution(fX);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
}

/** @brief update current state solution */
void TPZPoroPermAnalysis::Update_at_n_State(){
    this->LoadSolution(fX_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
}

void TPZPoroPermAnalysis::PostProcessStep(){
    
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    const int dim = this->Mesh()->Dimension();
    int div = 2;
    TPZStack<std::string>scalnames, vecnames;
    scalnames.Push("s_x");
    scalnames.Push("s_y");
    scalnames.Push("t_xy");    
//    scalnames.Push("k_x");
//    scalnames.Push("k_y");
    scalnames.Push("p_ex");
    vecnames.Push("u");
    
    std::string plotfile = "Poro_Permeability_2D.vtk";

    this->DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    this->PostProcess(div,dim);
    
}

/** @brief execute the evolutionary problem */
void TPZPoroPermAnalysis::Run_Evolution(){
    
    int n = fSimulationData->n_steps();
    for (int i = 0; i < n; i++) {
        
        this->ExcecuteOneStep();
        this->PostProcessStep();
        
    }
    
}