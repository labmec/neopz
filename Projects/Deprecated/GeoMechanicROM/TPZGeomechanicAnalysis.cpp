//
//  TPZGeomechanicAnalysis.cpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#include "TPZGeomechanicAnalysis.h"


TPZGeomechanicAnalysis::TPZGeomechanicAnalysis() : TPZAnalysis() {
    
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
    
    /** @brief Strain-Stress solution data */
    fstrain_stress_duplets.Resize(0);
    
    /** @brief Residue error */
    ferror    = 1.0;
    
    /** @brief Correction variation */
    fdx_norm  = 1.0;
    
    /** @brief number of newton corrections */
    fk_iterations = 0;
    
}

TPZGeomechanicAnalysis::~TPZGeomechanicAnalysis(){
    
}

/** @brief Copy constructor $ */
TPZGeomechanicAnalysis::TPZGeomechanicAnalysis(const TPZGeomechanicAnalysis &copy)
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
TPZGeomechanicAnalysis & TPZGeomechanicAnalysis::operator=(const TPZGeomechanicAnalysis &other)
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
void TPZGeomechanicAnalysis::AdjustVectors(){
    
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

void TPZGeomechanicAnalysis::QuasiNewtonIteration(){
    
    if(fk_iterations == 1){
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
    
    this->Update_at_n_State();
    
    this->AssembleResidual();
    fR_n = this->Rhs();
    fR_n += fR; // total residue
    ferror =  Norm(fR_n); // residue error
    
}

void TPZGeomechanicAnalysis::NewtonIteration(){
    
    this->Assemble();
    this->Rhs() += fR; // total residue
    this->Rhs() *= -1.0;
    
    this->Solve(); // update correction
    fdx_norm = Norm(this->Solution()); // correction variation
    
    fX_n += this->Solution(); // update state
    
    this->Update_at_n_State();
    
    this->AssembleResidual();
    fR_n = this->Rhs();
    fR_n += fR; // total residue
    ferror =  Norm(fR_n); // residue error
    
}

void TPZGeomechanicAnalysis::ExcecuteOneStep(){
    
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
//        this->QuasiNewtonIteration();
        this->NewtonIteration();

        
#ifdef PZDEBUG
//        this->Solver().Matrix()->Print("K = ", std::cout,EMathematicaInput);
//        fX.Print("X = ", std::cout,EMathematicaInput);
//        fR.Print("R = ", std::cout,EMathematicaInput);
//        fX_n.Print("Xn = ", std::cout,EMathematicaInput);
//        fR_n.Print("Rn = ", std::cout,EMathematicaInput);
#endif
        
        if(ferror < epsilon_res || fdx_norm < epsilon_cor)
        {
            std::cout << "Geomechanic Coupling:: Converged with iterations:  " << k << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
            fX = fX_n;
            return;
        }
        
    }
    
    std::cout << "Geomechanic Coupling:: Exit max iterations with min dt:  " << fSimulationData->dt() << "; (secs) " << "; error: " << ferror <<  "; dx: " << fdx_norm << std::endl;
    
    
}

/** @brief update last state solution */
void TPZGeomechanicAnalysis::UpdateState(){
    this->LoadSolution(fX);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
}

/** @brief update current state solution */
void TPZGeomechanicAnalysis::Update_at_n_State(){
    this->LoadSolution(fX_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
}

/** @brief PostProcess results for nonlinear case */
void TPZGeomechanicAnalysis::PostNonlinearProcessStep(std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    const int dim = this->Mesh()->Dimension();
    int div = 0;
    TPZStack<std::string>scalnames, vecnames;
    scalnames.Push("u");
    vecnames.Push("sigma");
    
    this->DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    this->PostProcess(div,dim);
    
}

void TPZGeomechanicAnalysis::PostProcessStep(std::string plotfile){
    
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, this->Mesh());
    const int dim = this->Mesh()->Dimension();
    int div = 2;
    TPZStack<std::string>scalnames, vecnames;
//    scalnames.Push("s_x");
//    scalnames.Push("s_y");
//    scalnames.Push("t_xy");
//    scalnames.Push("e_x");
//    scalnames.Push("e_y");
//    scalnames.Push("e_xy");
//    scalnames.Push("ep_x");
//    scalnames.Push("ep_y");
//    scalnames.Push("ep_xy");
//    scalnames.Push("k_x");
//    scalnames.Push("k_y");
//    scalnames.Push("K_0");
//    scalnames.Push("phi");    
    scalnames.Push("p_ex");
    vecnames.Push("u");
//    vecnames.Push("v");

    this->DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    this->PostProcess(div,dim);
    
}

/** @brief execute the evolutionary problem */
void TPZGeomechanicAnalysis::Run_Evolution(TPZVec<REAL> &x, std::string plotfile){
    
    int n = fSimulationData->n_steps();
    REAL time = 0.0;
    REAL dt = this->SimulationData()->dt();
    
    for (int i = 0; i <= n; i++) {
        time = i * dt;
        std::cout<< "Geomechanic Coupling:: Current time (s) = " << time << std::endl;
        this->SimulationData()->SetTime(time);
        this->ExcecuteOneStep();
        this->PostProcessStep(plotfile);
//        this->AppendStrain_Stress(x);
//        this->AppendStrain_Pososity(x);
//        this->AppendStrain_Permeability(x);
//        this->AppendStrain_Pressure(x);
        
    }
    
}

/** @brief Compute the strain and the stress at x euclidean point for each time */
void TPZGeomechanicAnalysis::AppendStrain_Stress(TPZVec<REAL> & x){
    
    
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    int64_t n_elemenst = geometry->NElements();
    
    int sx_var = 2;
    int sy_var = 3;
    int eey_var = 12;
    int epy_var = 15;
    TPZVec<STATE> sx;
    TPZVec<STATE> sy;
    TPZVec<STATE> eey;
    TPZVec<STATE> epy;
    
    std::pair<REAL,REAL> duplet;
    
    TPZVec<REAL> parametric_space(dim,0.0);
    for (int64_t iel = 0; iel < n_elemenst; iel++) {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif

        if (gel->Dimension() != dim) {
            continue;
        }
        

        IsTargetElementQ = gel->ComputeXInverse(x, parametric_space, Tol);
        
        if(IsTargetElementQ){
            TPZCompEl * cel = gel->Reference();
            cel->Solution(parametric_space, sx_var, sx);
            cel->Solution(parametric_space, sy_var, sy);
            cel->Solution(parametric_space, eey_var, eey);
            cel->Solution(parametric_space, epy_var, epy);
            duplet.first = eey[0] + epy[0];
            duplet.second = sy[0] - sx[0];
            
            duplet.first = fabs(duplet.first);
            duplet.second = fabs(duplet.second);
        }
        
    }

    fstrain_stress_duplets.Push(duplet);
}

/** @brief Compute the strain and the Pososity at x euclidean point for each time */
void TPZGeomechanicAnalysis::AppendStrain_Pososity(TPZVec<REAL> & x){
    
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    int64_t n_elemenst = geometry->NElements();
    
    int phi_var = 10;
    int eex_var = 11;
    int epx_var = 14;
    int eey_var = 12;
    int epy_var = 15;
    TPZVec<STATE> phi;
    TPZVec<STATE> eex;
    TPZVec<STATE> epx;
    TPZVec<STATE> eey;
    TPZVec<STATE> epy;
    
    std::pair<REAL,REAL> duplet;
    
    TPZVec<REAL> parametric_space(dim,0.0);
    for (int64_t iel = 0; iel < n_elemenst; iel++) {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension() != dim) {
            continue;
        }
        
        
        IsTargetElementQ = gel->ComputeXInverse(x, parametric_space, Tol);
        
        if(IsTargetElementQ){
            TPZCompEl * cel = gel->Reference();
            cel->Solution(parametric_space, phi_var, phi);
            cel->Solution(parametric_space, eex_var, eex);
            cel->Solution(parametric_space, epx_var, epx);
            cel->Solution(parametric_space, eey_var, eey);
            cel->Solution(parametric_space, epy_var, epy);
            duplet.first = eex[0] + epx[0] + eey[0] + epy[0];
            duplet.second = phi[0];
            
            duplet.first = fabs(duplet.first);
            duplet.second = fabs(duplet.second);
        }
        
    }
    
    fstrain_porosity_duplets.Push(duplet);
    
}

/** @brief Compute the strain and the Permeability at x euclidean point for each time */
void TPZGeomechanicAnalysis::AppendStrain_Permeability(TPZVec<REAL> & x){
    
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    int64_t n_elemenst = geometry->NElements();
    
    int k_var   = 9;
    int eex_var = 11;
    int epx_var = 14;
    int eey_var = 12;
    int epy_var = 15;
    TPZVec<STATE> k;
    TPZVec<STATE> eex;
    TPZVec<STATE> epx;
    TPZVec<STATE> eey;
    TPZVec<STATE> epy;
    
    std::pair<REAL,REAL> duplet;
    
    TPZVec<REAL> parametric_space(dim,0.0);
    for (int64_t iel = 0; iel < n_elemenst; iel++) {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension() != dim) {
            continue;
        }
        
        
        IsTargetElementQ = gel->ComputeXInverse(x, parametric_space, Tol);
        
        if(IsTargetElementQ){
            TPZCompEl * cel = gel->Reference();
            cel->Solution(parametric_space, k_var, k);
            cel->Solution(parametric_space, eex_var, eex);
            cel->Solution(parametric_space, epx_var, epx);
            cel->Solution(parametric_space, eey_var, eey);
            cel->Solution(parametric_space, epy_var, epy);
            duplet.first = eex[0] + epx[0] + eey[0] + epy[0];
            duplet.second = k[0];
            
            duplet.first = fabs(duplet.first);
            duplet.second = fabs(duplet.second);
        }
        
    }
    
    fstrain_permeability_duplets.Push(duplet);
    
}

/** @brief Compute the strain and the Permeability at x euclidean point for each time */
void TPZGeomechanicAnalysis::AppendStrain_Pressure(TPZVec<REAL> & x){
    
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    int64_t n_elemenst = geometry->NElements();
    
    int p_var   = 6;
    int eex_var = 11;
    int epx_var = 14;
    int eey_var = 12;
    int epy_var = 15;
    TPZVec<STATE> p;
    TPZVec<STATE> eex;
    TPZVec<STATE> epx;
    TPZVec<STATE> eey;
    TPZVec<STATE> epy;
    
    std::pair<REAL,REAL> duplet;
    
    TPZVec<REAL> parametric_space(dim,0.0);
    for (int64_t iel = 0; iel < n_elemenst; iel++) {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->Dimension() != dim) {
            continue;
        }
        
        
        IsTargetElementQ = gel->ComputeXInverse(x, parametric_space, Tol);
        
        if(IsTargetElementQ){
            TPZCompEl * cel = gel->Reference();
            cel->Solution(parametric_space, p_var, p);
            cel->Solution(parametric_space, eex_var, eex);
            cel->Solution(parametric_space, epx_var, epx);
            cel->Solution(parametric_space, eey_var, eey);
            cel->Solution(parametric_space, epy_var, epy);
            duplet.first = eex[0] + epx[0] + eey[0] + epy[0];
            duplet.second = p[0];
            
            duplet.first = fabs(duplet.first);
            duplet.second = fabs(duplet.second);
        }
        
    }
    
    fstrain_pressure_duplets.Push(duplet);
    
}

/** @brief Compute the strain and the stress at x euclidean point for each time */
void TPZGeomechanicAnalysis::PlotStrainStress(std::string file_name){
    
#ifdef PZDEBUG
    if (fstrain_stress_duplets.size() == 0) {
        DebugStop();
    }
#endif
    
    int n_data = fstrain_stress_duplets.size();
    TPZFMatrix<REAL> points(n_data,2,0.0);
    for(int i = 0; i < n_data; i++){
        points(i,0) = fstrain_stress_duplets[i].first;
        points(i,1) = fstrain_stress_duplets[i].second;
    }
    
    {
        std::ofstream out(file_name.c_str());
        points.Print("data = ",out,EMathematicaInput);
    }
}

/** @brief Compute the strain and the Porosity at x euclidean point for each time */
void TPZGeomechanicAnalysis::PlotStrainPorosity(std::string file_name){
    
#ifdef PZDEBUG
    if (fstrain_porosity_duplets.size() == 0) {
        DebugStop();
    }
#endif
    
    int n_data = fstrain_porosity_duplets.size();
    TPZFMatrix<REAL> points(n_data,2,0.0);
    for(int i = 0; i < n_data; i++){
        points(i,0) = fstrain_porosity_duplets[i].first;
        points(i,1) = fstrain_porosity_duplets[i].second;
    }
    
    {
        std::ofstream out(file_name.c_str());
        points.Print("phi = ",out,EMathematicaInput);
    }
    
}

/** @brief Compute the strain and the Porosity at x euclidean point for each time */
void TPZGeomechanicAnalysis::PlotStrainPermeability(std::string file_name){
    
#ifdef PZDEBUG
    if (fstrain_permeability_duplets.size() == 0) {
        DebugStop();
    }
#endif
    
    int n_data = fstrain_permeability_duplets.size();
    TPZFMatrix<REAL> points(n_data,2,0.0);
    for(int i = 0; i < n_data; i++){
        points(i,0) = fstrain_permeability_duplets[i].first;
        points(i,1) = fstrain_permeability_duplets[i].second;
    }
    
    {
        std::ofstream out(file_name.c_str());
        points.Print("k = ",out,EMathematicaInput);
    }
    
}

/** @brief Compute the strain and the Porosity at x euclidean point for each time */
void TPZGeomechanicAnalysis::PlotStrainPressure(std::string file_name){
    
#ifdef PZDEBUG
    if (fstrain_pressure_duplets.size() == 0) {
        DebugStop();
    }
#endif
    
    int n_data = fstrain_pressure_duplets.size();
    TPZFMatrix<REAL> points(n_data,2,0.0);
    for(int i = 0; i < n_data; i++){
        points(i,0) = fstrain_pressure_duplets[i].first;
        points(i,1) = fstrain_pressure_duplets[i].second;
    }
    
    {
        std::ofstream out(file_name.c_str());
        points.Print("p = ",out,EMathematicaInput);
    }
    
}