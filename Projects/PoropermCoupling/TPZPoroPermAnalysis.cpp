//
//  TPZPoroPermAnalysis.cpp
//  PZ
//
//  Created by Omar and Manouchehr on 8/28/16.
//
//

#include "TPZPoroPermAnalysis.h"
//#define PZDEBUG

/** @brief default costructor */
TPZPoroPermAnalysis::TPZPoroPermAnalysis() : TPZAnalysis()
{
    
    /** @brief define the simulation data */
    m_SimulationData = NULL;
    
    /** @brief Vector of compmesh pointers. m_meshvec[0] = flowHdiv, m_meshvec[1] = PressureL2 */
    m_meshvec.Resize(2);
    
    /** @brief Part of residue at n+1 (current) state */
    m_R_n.Resize(0,0);
    
    /** @brief Part of residue at n (last) state  */
    m_R.Resize(0,0);
    
    /** @brief Solution at n+1 (current) state */
    m_X_n.Resize(0,0);
    
    /** @brief Solution  at n (last) state */
    m_X.Resize(0,0);
    
    /** @brief Strain-Stress solution data */
    m_strain_stress_duplets.Resize(0);
    
    /** @brief Residue error */
    m_error    = 1.0;
    
    /** @brief Correction variation */
    m_dx_norm  = 1.0;
    
    /** @brief number of newton corrections */
    m_k_iterations = 0;
    
}

/** @brief default destructor $ */
TPZPoroPermAnalysis::~TPZPoroPermAnalysis(){
    
}

/** @brief copy constructor $ */
TPZPoroPermAnalysis::TPZPoroPermAnalysis(const TPZPoroPermAnalysis &copy)
{
    m_SimulationData = copy.m_SimulationData;
    m_meshvec        = copy.m_meshvec;
    m_R_n            = copy.m_R_n;
    m_R              = copy.m_R;
    m_X_n            = copy.m_X_n;
    m_X              = copy.m_X;
    m_error          = copy.m_error;
    m_dx_norm        = copy.m_dx_norm;
    
}

/** @brief Copy assignemnt operator $ */
TPZPoroPermAnalysis & TPZPoroPermAnalysis::operator=(const TPZPoroPermAnalysis &other)
{
    if (this != & other) {  // prevent self-assignment
        
        m_SimulationData = other.m_SimulationData;
        m_meshvec        = other.m_meshvec;
        m_R_n            = other.m_R_n;
        m_R              = other.m_R;
        m_X_n            = other.m_X_n;
        m_X              = other.m_X;
        m_error          = other.m_error;
        m_dx_norm        = other.m_dx_norm;
    }
    return *this;
}

/** @brief Resize and fill residue and solution vectors */
void TPZPoroPermAnalysis::AdjustVectors(){
    
    if(fSolution.Rows() == 0 /* || fRhs.Rows() == 0 */)
    {
        DebugStop();
    }
    
    TPZBuildMultiphysicsMesh::AddElements(m_meshvec, this->Mesh());
    TPZBuildMultiphysicsMesh::AddConnects(m_meshvec, this->Mesh());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(m_meshvec, this->Mesh());
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
    
    m_X.Resize(fSolution.Rows(),1);
    m_X.Zero();
    m_X_n.Resize(fSolution.Rows(),1);
    m_X_n.Zero();
    m_R_n.Resize(fSolution.Rows(),1);
    m_R_n.Zero();
    m_R.Resize(fSolution.Rows(),1);
    m_R.Zero();
}

void TPZPoroPermAnalysis::QuasiNewtonIteration()
{
    
    if(m_k_iterations == 1)
    {
        this->Assemble();
    }
    else
    {
        this->AssembleResidual();
    }
    this->Rhs() += m_R; // total residue
    this->Rhs() *= -1.0;

#ifdef PZDEBUG
//        this->Solver().Matrix()->Print("J = ", std::cout,EMathematicaInput);
//        this->Rhs().Print("R = ", std::cout,EMathematicaInput);
#endif
    
    this->Solve(); // update correction
    m_dx_norm = Norm(this->Solution()); // correction variation
    
#ifdef PZDEBUG
//    this->Solution().Print("dx = ", std::cout,EMathematicaInput);
#endif
    
    m_X_n += this->Solution(); // update state

    this->Update_at_n_State();
    
    this->AssembleResidual();
    m_R_n = this->Rhs();
    
#ifdef PZDEBUG
//    m_X.Print("X = ", std::cout,EMathematicaInput);
//    m_X_n.Print("Xn = ", std::cout,EMathematicaInput);
//    m_R.Print("R = ", std::cout,EMathematicaInput);
//    m_R_n.Print("Rn = ", std::cout,EMathematicaInput);
#endif
    
    m_R_n += m_R; // total residue
#ifdef PZDEBUG
//    m_R_n.Print("Rt = ", std::cout,EMathematicaInput);
#endif
    m_error =  Norm(m_R_n); // residue error
    
}

void TPZPoroPermAnalysis::ExcecuteOneStep(){
    
    this->SimulationData()->SetCurrentStateQ(false);
    this->UpdateState();

    this->AssembleResidual();
    m_R = this->Rhs();
    
    this->SimulationData()->SetCurrentStateQ(true);
    
    this->Update_at_n_State();
    
    m_error = 1.0;
    
    STATE epsilon_res = this->SimulationData()->epsilon_res();
    STATE epsilon_cor = this->SimulationData()->epsilon_cor();
    int n  =   this->SimulationData()->n_iterations();
    
    for (int k = 1; k <= n; k++)
    {
                
        this->Set_k_ietrarions(k);
        this->QuasiNewtonIteration();
        
        if(m_error < epsilon_res || (m_dx_norm < epsilon_cor && k > 3 ) )
        {
            std::cout << "Permeability Coupling:: Converged with iterations:  " << k << "; error: " << m_error <<  "; dx: " << m_dx_norm << std::endl;
            m_X = m_X_n;
            return;
        }
        
    }
    
    std::cout << "Permeability Coupling:: Exit max iterations with min dt:  " << m_SimulationData->dt() << "; (secs) " << "; error: " << m_error <<  "; dx: " << m_dx_norm << std::endl;
    
    
}

/** @brief update last state (at n state) solution */
void TPZPoroPermAnalysis::UpdateState()
{
    this->LoadSolution(m_X);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
}

/** @brief update current state (at n+1 state) solution */
void TPZPoroPermAnalysis::Update_at_n_State()
{
    this->LoadSolution(m_X_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
}

void TPZPoroPermAnalysis::PostProcessStep()
{
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_meshvec, this->Mesh());
    const int dim = this->Mesh()->Dimension();
    int div = m_SimulationData->n_div();
 
    TPZManVector<std::string,50> scalnames = m_SimulationData->scalar_names();
    TPZManVector<std::string,50> vecnames = m_SimulationData->vector_names();
    
    std::string plotfile = m_SimulationData->name_vtk_file();

    this->DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    this->PostProcess(div,dim);
    
}

/** @brief execute the evolutionary problem */
void TPZPoroPermAnalysis::Run_Evolution(TPZVec<REAL> &x)
{
    
    int n = m_SimulationData->n_steps();
    REAL time = 0.0;
    REAL dt = this->SimulationData()->dt();

    for (int i = 0; i < n; i++)
    {
        this->ExcecuteOneStep();
        this->PostProcessStep();
//        this->AppendStrain_Stress(x);
//        this->AppendStrain_Pososity(x);
//        this->AppendStrain_Permeability(x);
//        this->AppendStrain_Pressure(x);
        time = (i+1)* dt;
        std::cout<< "Permeability Coupling:: Current time (s) = " << time << std::endl;
        this->SimulationData()->SetTime(time);

        
    }
    
}

/** @brief Compute the strain and the stress at x euclidean point for each time */
void TPZPoroPermAnalysis::AppendStrain_Stress(TPZVec<REAL> & x)
{
    
    
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    long n_elemenst = geometry->NElements();
    
    int sx_var = 2;
    int sy_var = 3;
//    int eey_var = 12;
//    int epy_var = 15;
    int eey_var = 11;
    int epy_var = 14;
    TPZVec<STATE> sx;
    TPZVec<STATE> sy;
    TPZVec<STATE> eey;
    TPZVec<STATE> epy;
    
    std::pair<REAL,REAL> duplet;
    
    TPZVec<REAL> parametric_space(dim,0.0);
    for (long iel = 0; iel < n_elemenst; iel++)
    {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel)
        {
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
            
            duplet.first = -duplet.first;
            duplet.second = fabs(duplet.second);
        }
        
    }

    m_strain_stress_duplets.Push(duplet);
}

/** @brief Compute the strain and the Pososity at x euclidean point for each time */
void TPZPoroPermAnalysis::AppendStrain_Pososity(TPZVec<REAL> & x)
{
    
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    long n_elemenst = geometry->NElements();
    
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
    for (long iel = 0; iel < n_elemenst; iel++)
    {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel)
        {
            DebugStop();
        }
#endif
        
        if (gel->Dimension() != dim)
        {
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
    
    m_strain_porosity_duplets.Push(duplet);
    
}

/** @brief Compute the strain and the Permeability at x euclidean point for each time */
void TPZPoroPermAnalysis::AppendStrain_Permeability(TPZVec<REAL> & x)
{
    
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    long n_elemenst = geometry->NElements();
    
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
    for (long iel = 0; iel < n_elemenst; iel++)
    {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel)
        {
            DebugStop();
        }
#endif
        
        if (gel->Dimension() != dim)
        {
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
    
   m_strain_permeability_duplets.Push(duplet);
    
}

/** @brief Compute the strain and the Permeability at x euclidean point for each time */
void TPZPoroPermAnalysis::AppendStrain_Pressure(TPZVec<REAL> & x)
{
    
    // Finding the geometic element that x bleongs to.
    REAL Tol = 1.0e-4;
    TPZGeoMesh * geometry = this->Mesh()->Reference();
    this->Mesh()->LoadReferences();
    
    int dim = geometry->Dimension();
    bool IsTargetElementQ = false;
    long n_elemenst = geometry->NElements();
    
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
    for (long iel = 0; iel < n_elemenst; iel++)
    {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if (!gel)
        {
            DebugStop();
        }
#endif
        
        if (gel->Dimension() != dim)
        {
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
    
    m_strain_pressure_duplets.Push(duplet);
    
}

/** @brief Compute the strain and the stress at x euclidean point for each time */
void TPZPoroPermAnalysis::PlotStrainStress(std::string file_name)
{
    
#ifdef PZDEBUG
    if (m_strain_stress_duplets.size() == 0)
    {
        DebugStop();
    }
#endif
    
    int n_data = m_strain_stress_duplets.size();
    TPZFMatrix<REAL> points(n_data,2,0.0);
    for(int i = 0; i < n_data; i++){
        points(i,0) = m_strain_stress_duplets[i].first;
        points(i,1) = m_strain_stress_duplets[i].second;
    }
    
    {
        std::ofstream out(file_name.c_str());
        points.Print("data = ",out,EMathematicaInput);
    }
}

/** @brief Compute the strain and the Porosity at x euclidean point for each time */
void TPZPoroPermAnalysis::PlotStrainPorosity(std::string file_name)
{
    
#ifdef PZDEBUG
    if (m_strain_porosity_duplets.size() == 0)
    {
        DebugStop();
    }
#endif
    
    int n_data = m_strain_porosity_duplets.size();
    TPZFMatrix<REAL> points(n_data,2,0.0);
    for(int i = 0; i < n_data; i++)
    {
        points(i,0) = m_strain_porosity_duplets[i].first;
        points(i,1) = m_strain_porosity_duplets[i].second;
    }
    
    {
        std::ofstream out(file_name.c_str());
        points.Print("phi = ",out,EMathematicaInput);
    }
    
}

/** @brief Compute the strain and the Porosity at x euclidean point for each time */
void TPZPoroPermAnalysis::PlotStrainPermeability(std::string file_name)
{
    
#ifdef PZDEBUG
    if (m_strain_permeability_duplets.size() == 0)
    {
        DebugStop();
    }
#endif
    
    int n_data = m_strain_permeability_duplets.size();
    TPZFMatrix<REAL> points(n_data,2,0.0);
    for(int i = 0; i < n_data; i++){
        points(i,0) = m_strain_permeability_duplets[i].first;
        points(i,1) = m_strain_permeability_duplets[i].second;
    }
    
    {
        std::ofstream out(file_name.c_str());
        points.Print("k = ",out,EMathematicaInput);
    }
    
}

/** @brief Compute the strain and the Porosity at x euclidean point for each time */
void TPZPoroPermAnalysis::PlotStrainPressure(std::string file_name)
{
    
#ifdef PZDEBUG
    if (m_strain_pressure_duplets.size() == 0)
    {
        DebugStop();
    }
#endif
    
    int n_data = m_strain_pressure_duplets.size();
    TPZFMatrix<REAL> points(n_data,2,0.0);
    for(int i = 0; i < n_data; i++)
    {
        points(i,0) = m_strain_pressure_duplets[i].first;
        points(i,1) = m_strain_pressure_duplets[i].second;
    }
    
    {
        std::ofstream out(file_name.c_str());
        points.Print("p = ",out,EMathematicaInput);
    }
    
}
