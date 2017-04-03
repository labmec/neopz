//
//  TPZElasticAnalysis.cpp
//  PZ
//
//  Created by Omar on 3/5/17.
//
//

#include "TPZElasticAnalysis.h"

TPZElasticAnalysis::TPZElasticAnalysis() : TPZAnalysis() {
    
    /** @brief define the simulation data */
    fSimulationData = NULL;
    
    /** @brief define the tranfer object data */
    ftransfer = NULL;
    
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

TPZElasticAnalysis::~TPZElasticAnalysis(){
    
}

/** @brief Copy constructor $ */
TPZElasticAnalysis::TPZElasticAnalysis(const TPZElasticAnalysis &copy)
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
TPZElasticAnalysis & TPZElasticAnalysis::operator=(const TPZElasticAnalysis &other)
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
void TPZElasticAnalysis::AdjustVectors(){
    
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

void TPZElasticAnalysis::QuasiNewtonIteration(){
    
    if(fk_iterations == 1){
        this->Assemble();
    }
    else{
        this->AssembleResidual();
    }
    this->Rhs() += fR; // total residue
    this->Rhs() *= -1.0;
    
#ifdef PZDEBUG
//            this->Solver().Matrix()->Print("K = ", std::cout,EMathematicaInput);
//            this->Rhs().Print("Rn = ", std::cout,EMathematicaInput);
#endif
    
    this->Solve(); // update correction
    fdx_norm = Norm(this->Solution()); // correction variation
    
    fX_n += this->Solution(); // update state
    
    this->Update_at_n_State();
    
    this->AssembleResidual();
    fR_n = this->Rhs();
    fR_n += fR; // total residue
    ferror =  Norm(fR_n); // residue error
    
}

void TPZElasticAnalysis::NewtonIteration(){
    
    this->Assemble();
    this->Rhs() += fR; // total residue
    this->Rhs() *= -1.0;
    
#ifdef PZDEBUG
//        this->Solver().Matrix()->Print("K = ", std::cout,EMathematicaInput);
//        this->Rhs().Print("Rn = ", std::cout,EMathematicaInput);
#endif
    
    this->Solve(); // update correction
    fdx_norm = Norm(this->Solution()); // correction variation
    
    fX_n += this->Solution(); // update state
    
    this->Update_at_n_State();
    
    this->AssembleResidual();
    fR_n = this->Rhs();
    fR_n += fR; // total residue
    ferror =  Norm(fR_n); // residue error
    
}

void TPZElasticAnalysis::ExcecuteOneStep(){
    
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
        //        this->NewtonIteration();
        if(ferror < epsilon_res || fdx_norm < epsilon_cor)
        {
            std::cout << "elliptic:: Converged with iterations:  " << k << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
//            fX = fX_n;
            return;
        }
        
    }
    
    std::cout << "elliptic:: Exit max iterations with min dt:  " << fSimulationData->dt() << "; (secs) " << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
    
}

/** @brief update last state solution */
void TPZElasticAnalysis::UpdateState(){
    this->LoadSolution(fX);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    if (fSimulationData->IsRBApproxQ() && !fSimulationData->IsInitialStateQ()) {
        ftransfer->rb_elliptic_To_rb_elliptic(this->Mesh());
    }
    else{
       ftransfer->elliptic_To_elliptic(this->Mesh());
    }
    
}

/** @brief update current state solution */
void TPZElasticAnalysis::Update_at_n_State(){
    this->LoadSolution(fX_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    if (fSimulationData->IsRBApproxQ() && !fSimulationData->IsInitialStateQ()) {
        ftransfer->rb_elliptic_To_rb_elliptic(this->Mesh());
    }
    else{
        ftransfer->elliptic_To_elliptic(this->Mesh());
    }
}


void TPZElasticAnalysis::PostProcessStep(std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    const int dim = this->Mesh()->Dimension();
    int div = 5;
    TPZStack<std::string>scalnames, vecnames;
    scalnames.Push("sx");
    scalnames.Push("sy");
    scalnames.Push("sxy");
    vecnames.Push("u");
    vecnames.Push("ue");
    
    this->DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    this->PostProcess(div,dim);
}
