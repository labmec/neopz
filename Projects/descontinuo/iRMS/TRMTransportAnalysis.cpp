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
    ferror    = 0.0;
    
    /** @brief Correction variation */
    fdx_norm  = 0.0;
    
    /** @brief number of newton corrections */
    fk_iterations = 0;
    
    /** @brief active equations of alpha saturation */
    factive_sa.Resize(0);
    
    /** @brief no active equations of alpha saturation */
    fno_active_sa.Resize(0);
    
    /** @brief active equations of beta saturation */
    factive_sb.Resize(0);
    
    /** @brief no active equations of beta saturation */
    fno_active_sb.Resize(0);
    
    /** @brief Gradient reconstruction object */
    fgradreconst = new TPZGradientReconstruction(false,1.);
    
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
    
    if (fSimulationData->UseGradientR())
    {
        //        CleanUpGradients();
        this->SaturationReconstruction();
        fX_n = Solution();
    }
    
    this->UpdateMemory_at_n();
    
    this->AssembleResidual();
    fR_n = this->Rhs();
    fR_n += fR; // total residue
    ferror =  Norm(fR_n); // residue error
    
}

void TRMTransportAnalysis::QuasiNewtonIteration(){
    
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
    if (fSimulationData->UseGradientR())
    {
        CleanUpGradients();
        this->SaturationReconstruction();
        fX_n = Solution();
    }
    this->UpdateMemory_at_n();
    
    this->AssembleResidual();
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
    
    this->UpdateMemory_at_n();    
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());

    this->AssembleResidual();
    fR_n = this->Rhs();
    
    ferror = Norm(fR+fR_n);
    this->Set_k_ietrarions(0);
    
    REAL epsilon_res = this->SimulationData()->epsilon_res() * 0.01;
    REAL epsilon_cor = this->SimulationData()->epsilon_cor() * 0.01;
    int n  =   this->SimulationData()->n_corrections();
    
    for (int k = 1; k <= n; k++) {

        this->Set_k_ietrarions(k);
        this->NewtonIteration();// @omar:: I prefer no linearize this matrix
        
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
            return;
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
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, Mesh());    
    const int dim = this->Mesh()->Reference()->Dimension();
    int div = 0;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "DualSegregatedDarcyOnBox_Saturations.vtk";
    scalnames.Push("s_a");
    scalnames.Push("s_b");
    
    if (fSimulationData->IsInitialStateQ()) {
        plotfile =  "DualSegregatedDarcyOnBox_I.vtk";
        return;
    }
    
    if (fSimulationData->IsThreePhaseQ()) {
        scalnames.Push("s_c");
    }
    
    this->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    this->PostProcess(div);
    
}

void TRMTransportAnalysis::FilterSaturationGradients()
{
    
    int ncon_sa = fmeshvec[0]->NConnects();
    int ncon = Mesh()->NConnects();
    
    
    // DOF related with Constant
    for(int i = 0; i < ncon_sa; i++)
    {
        TPZConnect &con = Mesh()->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = Mesh()->Block().Position(seqnum);
        int blocksize = Mesh()->Block().Size(seqnum);
        int vs = factive_sa.size();
        factive_sa.Resize(vs+1);
        
        int ieq = blocksize-1;
        factive_sa[vs] = pos+ieq;
    }
    
    // DOF Related with S grandients
    for(int i = 0; i< ncon_sa; i++)
    {
        TPZConnect &con = Mesh()->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = Mesh()->Block().Position(seqnum);
        int blocksize = Mesh()->Block().Size(seqnum);
        int vs = fno_active_sa.size();
        fno_active_sa.Resize(vs+blocksize-1);
        for(int ieq = 0; ieq<blocksize-1; ieq++)
        {
            fno_active_sa[vs+ieq] = pos+ieq;
        }
        
    }
    
    if(fSimulationData->IsTwoPhaseQ()){
        return;
    }
    
    // DOF related with Constant
    for(int i = ncon - ncon_sa; i < ncon; i++)
    {
        TPZConnect &con = Mesh()->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = Mesh()->Block().Position(seqnum);
        int blocksize = Mesh()->Block().Size(seqnum);
        int vs = factive_sa.size();
        factive_sa.Resize(vs+1);
        
        int ieq = blocksize-1;
        factive_sa[vs] = pos+ieq;
    }
    
    // DOF Related with S grandients
    for(int i = ncon - ncon_sa; i < ncon; i++)
    {
        TPZConnect &con = Mesh()->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = Mesh()->Block().Position(seqnum);
        int blocksize = Mesh()->Block().Size(seqnum);
        int vs = fno_active_sa.size();
        fno_active_sa.Resize(vs+blocksize-1);
        for(int ieq = 0; ieq<blocksize-1; ieq++)
        {
            fno_active_sa[vs+ieq] = pos+ieq;
        }
        
    }
    
}

void TRMTransportAnalysis::CleanUpGradients(){
    
    int64_t numofdof_sa = fno_active_sa.size();
    int64_t numofdof_sb = fno_active_sb.size();
    TPZFMatrix<REAL> SolToLoad = Solution();
    for(int64_t i=0; i < numofdof_sa; i++)
    {
        Solution()(fno_active_sa[i],0) = 0.0;
    }
    for(int64_t i=0; i < numofdof_sb; i++)
    {
        Solution()(fno_active_sb[i],0) = 0.0;
    }
    
}

void TRMTransportAnalysis::SaturationReconstruction()
{
    Mesh()->LoadSolution(Solution());
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, Mesh());
    
    if(fSimulationData->IsTwoPhaseQ()){
        fmeshvec[0]->Reference()->ResetReference();
        fmeshvec[0]->LoadReferences();
        fgradreconst->ProjectionL2GradientReconstructed(fmeshvec[0], fSimulationData->L2_Projection_material_Id());
    }
    
    if(fSimulationData->IsThreePhaseQ()){
        fmeshvec[0]->Reference()->ResetReference();
        fmeshvec[0]->LoadReferences();
        fgradreconst->ProjectionL2GradientReconstructed(fmeshvec[0], fSimulationData->L2_Projection_material_Id());
        
        fmeshvec[1]->Reference()->ResetReference();
        fmeshvec[1]->LoadReferences();
        fgradreconst->ProjectionL2GradientReconstructed(fmeshvec[1], fSimulationData->L2_Projection_material_Id());
    }
    
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, Mesh());
    LoadSolution(Mesh()->Solution());
    
}


void TRMTransportAnalysis::FilterEquations(){
    FilterSaturationGradients();
    this->StructMatrix()->EquationFilter().Reset();
    this->StructMatrix()->EquationFilter().SetActiveEquations(factive_sa);
}
