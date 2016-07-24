//
//  TRMMixedDarcy.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMMixedDarcy.h"


TRMMixedDarcy::TRMMixedDarcy() : TPZMatWithMem<TRMMemory, TPZDiscontinuousGalerkin>()
{
    
}

TRMMixedDarcy::TRMMixedDarcy(int matid) : TPZMatWithMem<TRMMemory, TPZDiscontinuousGalerkin>(matid)
{
    
}


TRMMixedDarcy::TRMMixedDarcy(const TRMMixedDarcy &mat) : TPZMatWithMem<TRMMemory, TPZDiscontinuousGalerkin>(mat)
{
    
}

TRMMixedDarcy::~TRMMixedDarcy()
{
    
}

void TRMMixedDarcy::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
    
}

void TRMMixedDarcy::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TRMMixedDarcy::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

int TRMMixedDarcy::VariableIndex(const std::string &name) {
    if (!strcmp("p", name.c_str())) return 0;
    if (!strcmp("u", name.c_str())) return 1;
    if (!strcmp("div_u", name.c_str())) return 2;
    if (!strcmp("s_a", name.c_str())) return 2;    
    return TPZMatWithMem::VariableIndex(name);
}

int TRMMixedDarcy::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return 3; // Vector
        case 2:
            return 1; // Scalar
    }
    return TPZMatWithMem::NSolutionVariables(var);
}

void TRMMixedDarcy::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            Solution_a(datavec, var, Solout);
        }
            break;
        case 2:
        {
            Solution_ab(datavec, var, Solout);
        }
            break;
        case 3:
        {
            DebugStop();
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
}

// Jacobian contribution
void TRMMixedDarcy::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            Contribute_a(datavec, weight, ek, ef);
        }
            break;
        case 2:
        {
            Contribute_ab(datavec, weight, ek, ef);
        }
            break;
        case 3:
        {
            DebugStop();
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
}


void TRMMixedDarcy::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            Contribute_a(datavec, weight, ef);
        }
            break;
        case 2:
        {
            Contribute_ab(datavec, weight, ef);
        }
            break;
        case 3:
        {
            DebugStop();
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
}

void TRMMixedDarcy::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            DebugStop();
        }
            break;
        case 3:
        {
            DebugStop();
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
}

void TRMMixedDarcy::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef)
{
    std::cout << " This method should be called only for Capillary pressure terms " << std::endl;
    DebugStop();
}


void TRMMixedDarcy::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            DebugStop();
        }
            break;
        case 3:
        {
            DebugStop();
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
}

void TRMMixedDarcy::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    std::cout << " This method should be called only for Capillary pressure terms " << std::endl;
    DebugStop();
}


void TRMMixedDarcy::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            ContributeBC_a(datavec, weight, ek, ef, bc);
        }
            break;
        case 2:
        {
            ContributeBC_ab(datavec, weight, ek, ef, bc);
        }
            break;
        case 3:
        {
            DebugStop();
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
}


// Divergence on master element

void TRMMixedDarcy::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU)
{
    int ublock = 0;
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;   // For H1  test functions Q
    TPZFMatrix<STATE> dphiuH1       = datavec[ublock].dphi; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiuH1axes   = datavec[ublock].dphix; // Derivative For H1  test functions
    TPZFNMatrix<9,STATE> gradu = datavec[ublock].dsol[0];
    TPZFNMatrix<9,STATE> graduMaster;
    gradu.Transpose();
    
    TPZFNMatrix<660> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[ublock].axes);
    
    int nphiuHdiv = datavec[ublock].fVecShapeIndex.NElements();
    
    DivergenceofPhi.Resize(nphiuHdiv,1);
    
    REAL JacobianDet = datavec[ublock].detjac;
    
    TPZFMatrix<STATE> Qaxes = datavec[ublock].axes;
    TPZFMatrix<STATE> QaxesT;
    TPZFMatrix<STATE> Jacobian = datavec[ublock].jacobian;
    TPZFMatrix<STATE> JacobianInverse = datavec[ublock].jacinv;
    
    TPZFMatrix<STATE> GradOfX;
    TPZFMatrix<STATE> GradOfXInverse;
    TPZFMatrix<STATE> VectorOnMaster;
    TPZFMatrix<STATE> VectorOnXYZ(3,1,0.0);
    Qaxes.Transpose(&QaxesT);
    QaxesT.Multiply(Jacobian, GradOfX);
    JacobianInverse.Multiply(Qaxes, GradOfXInverse);
    
    int ivectorindex = 0;
    int ishapeindex = 0;
    
    if (HDivPiola == 1)
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            VectorOnXYZ(0,0) = datavec[ublock].fNormalVec(0,ivectorindex);
            VectorOnXYZ(1,0) = datavec[ublock].fNormalVec(1,ivectorindex);
            VectorOnXYZ(2,0) = datavec[ublock].fNormalVec(2,ivectorindex);
            
            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
            VectorOnMaster *= JacobianDet;
            
            /* Contravariant Piola mapping preserves the divergence */
            // the division by the jacobianDet is to make the integral on the master element???
            DivergenceofPhi(iq,0) = ( dphiuH1(0,ishapeindex)*VectorOnMaster(0,0) +
                                     dphiuH1(1,ishapeindex)*VectorOnMaster(1,0) +
                                     dphiuH1(2,ishapeindex)*VectorOnMaster(2,0) );
        }
        
        GradOfXInverse.Multiply(gradu, graduMaster);
        graduMaster *= JacobianDet;
        DivergenceofU = (graduMaster(0,0)+graduMaster(1,1)+graduMaster(2,2));
    }
    else
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            /* Computing the divergence for constant jacobian elements */
            DivergenceofPhi(iq,0) =  datavec[ublock].fNormalVec(0,ivectorindex)*GradphiuH1(0,ishapeindex) +
            datavec[ublock].fNormalVec(1,ivectorindex)*GradphiuH1(1,ishapeindex) +
            datavec[ublock].fNormalVec(2,ivectorindex)*GradphiuH1(2,ishapeindex) ;
        }
    }
    
    return;
    
}

// ------------------------------------------------------------------- //
// one phase flow case
// ------------------------------------------------------------------- //

void TRMMixedDarcy::Contribute_a(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    
    int nvars = 4; // {p,sa,sb,t}
    
    int ub = 0;
    int pb = 1;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,STATE> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,STATE> dphi_us      = datavec[ub].dphix;
    TPZFNMatrix<100,STATE> dphi_ps      = datavec[pb].dphix;
    
    TPZFNMatrix<40,STATE> div_on_master;
    STATE divflux;
    this->ComputeDivergenceOnMaster(datavec, div_on_master,divflux);
    REAL jac_det = datavec[ub].detjac;
    
    int nphiu       = datavec[ub].fVecShapeIndex.NElements();
    int nphip       = phi_ps.Rows();
    int firstu      = 0;
    int firstp      = nphiu + firstu;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    REAL p                  = datavec[pb].sol[0][0];
    
    TPZFNMatrix<10,STATE> Graduaxes = datavec[ub].dsol[0];
    
    // Time
    STATE dt = fSimulationData->dt();
    
    //  Average values p_a
    
    // Get the pressure at the integrations points
    long global_point_index = datavec[0].intGlobPtIndex;
    TRMMemory &point_memory = GetMemory()[global_point_index];
    REAL p_avg_n    = point_memory.p_avg_n();
    
    //  Computing closure relationship at given average values
    
    TPZManVector<STATE, 10> v(nvars),v_avg(nvars);
    v[0]        = p;
    v_avg[0]    = p_avg_n;
    
    // Fluid parameters
    TPZManVector<STATE, 10> rho,l;
    fSimulationData->AlphaProp()->Density(rho, v);
    fSimulationData->PetroPhysics()->l(l, v_avg);
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    TPZManVector<STATE, 10> phi;
    fSimulationData->Map()->Kappa(datavec[ub].x, K, Kinv, v);
    fSimulationData->Map()->phi(datavec[ub].x, phi, v);
    
    // Defining local variables
    TPZFNMatrix<3,STATE> lambda_K_inv_u(3,1),lambda_dp_K_inv_u(3,1), lambda_K_inv_phi_u_j(3,1);
    TPZManVector<STATE,3> Gravity = fSimulationData->Gravity();
    
    for (int i = 0; i < u.size(); i++) {
        STATE dot = 0.0;
        for (int j =0; j < u.size(); j++) {
            dot += Kinv(i,j)*u[j];
        }
        lambda_K_inv_u(i,0)     = (1.0/l[0]) * dot;
        lambda_dp_K_inv_u(i,0)  = (-l[1]/(l[0]*l[0])) * dot;
    }
    
    // Integration point contribution
    STATE divu = 0.0;
    TPZFNMatrix<3,STATE> phi_u_i(3,1), phi_u_j(3,1);
    
    int s_i, s_j;
    int v_i, v_j;
    
    if(! fSimulationData->IsCurrentStateQ()){
        
        REAL p_a    = point_memory.p_avg();
        
        TPZManVector<STATE, 10> v(nvars);
        v[0] = p;//p_a;
        
        // Fluid parameters
        TPZManVector<STATE, 10> rho,l;
        fSimulationData->AlphaProp()->Density(rho, v);
        
        // Rock parameters
        TPZManVector<STATE, 10> phi;
        fSimulationData->Map()->phi(datavec[ub].x, phi, v);
        
        for (int ip = 0; ip < nphip; ip++)
        {
            
            ef(ip + firstp) += -1.0 * weight * (-1.0/dt) * rho[0] * phi[0] * phi_ps(ip,0);
            
        }
        
        return;
    }
    
    for (int iu = 0; iu < nphiu; iu++)
    {
        
        v_i = datavec[ub].fVecShapeIndex[iu].first;
        s_i = datavec[ub].fVecShapeIndex[iu].second;
        
        STATE Kl_inv_dot_u = 0.0, Kl_dp_inv_dot_u = 0.0, rho_g_dot_phi_u = 0.0, rho_dp_g_dot_phi_u = 0.0;
        for (int i = 0; i < u.size(); i++) {
            phi_u_i(i,0) = phi_us(s_i,0) * datavec[ub].fNormalVec(i,v_i);
            Kl_inv_dot_u        += lambda_K_inv_u(i,0)*phi_u_i(i,0);
            Kl_dp_inv_dot_u     += lambda_dp_K_inv_u(i,0)*phi_u_i(i,0);
            rho_g_dot_phi_u     += rho[0]*Gravity[i]*phi_u_i(i,0);
            rho_dp_g_dot_phi_u  += rho[1]*Gravity[i]*phi_u_i(i,0);
        }
        
        ef(iu + firstu) += weight * ( Kl_inv_dot_u - (1.0/jac_det) * (p) * div_on_master(iu,0) - rho_g_dot_phi_u);
        
        for (int ju = 0; ju < nphiu; ju++)
        {
            
            v_j = datavec[ub].fVecShapeIndex[ju].first;
            s_j = datavec[ub].fVecShapeIndex[ju].second;
            
            STATE Kl_inv_phi_u_j_dot_phi_u_j = 0.0;
            for (int j = 0; j < u.size(); j++) {
                phi_u_j(j,0) = phi_us(s_j,0) * datavec[ub].fNormalVec(j,v_j);
                STATE dot = 0.0;
                for (int k = 0; k < u.size(); k++) {
                    dot += (1.0/l[0]) * Kinv(j,k)*phi_u_j(k,0);
                }
                lambda_K_inv_phi_u_j(j,0) = dot;
                Kl_inv_phi_u_j_dot_phi_u_j += lambda_K_inv_phi_u_j(j,0)*phi_u_i(j,0);
            }
            
            
            ek(iu + firstu,ju + firstu) += weight * Kl_inv_phi_u_j_dot_phi_u_j;
        }
        
        for (int jp = 0; jp < nphip; jp++)
        {
            ek(iu + firstu, jp + firstp) += weight * ( Kl_dp_inv_dot_u - (1.0/jac_det) * div_on_master(iu,0) + rho_dp_g_dot_phi_u) * phi_ps(jp,0);
        }
        
    }
    
    
    TPZManVector<STATE,1> f(1,0.0);
    if(fForcingFunction)
    {
        fForcingFunction->Execute(datavec[pb].x,f);
    }
    
    
    divu = (Graduaxes(0,0) + Graduaxes(1,1) + Graduaxes(2,2));
    
    for (int ip = 0; ip < nphip; ip++)
    {
        
        ef(ip + firstp) += -1.0 * weight * (divu + (1.0/dt) * rho[0] * phi[0] - f[0]) * phi_ps(ip,0);
        
        for (int ju = 0; ju < nphiu; ju++)
        {
            ek(ip + firstp, ju + firstu) += -1.0 * weight * (1.0/jac_det) * div_on_master(ju,0) * phi_ps(ip,0);
        }
        
        for (int jp = 0; jp < nphip; jp++)
        {
            ek(ip + firstp, jp + firstp) += -1.0 * weight * ( (1.0/dt) * (rho[0] * phi[1] + rho[1] * phi[0]) * phi_ps(ip,0) ) * phi_ps(jp,0);
        }
        
    }
    
}

void TRMMixedDarcy::Contribute_a(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    int nvars = 4; // {p,sa,sb,t}
    
    int ub = 0;
    int pb = 1;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,STATE> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,STATE> dphi_us      = datavec[ub].dphix;
    TPZFNMatrix<100,STATE> dphi_ps      = datavec[pb].dphix;
    
    TPZFNMatrix<40,STATE> div_on_master;
    STATE divflux;
    this->ComputeDivergenceOnMaster(datavec, div_on_master,divflux);
    REAL jac_det = datavec[ub].detjac;
    
    int nphiu       = datavec[ub].fVecShapeIndex.NElements();
    int nphip       = phi_ps.Rows();
    int firstu      = 0;
    int firstp      = nphiu + firstu;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    REAL p                  = datavec[pb].sol[0][0];
    
    TPZFNMatrix<10,STATE> Graduaxes = datavec[ub].dsol[0];
    
    // Time
    STATE dt = fSimulationData->dt();
    
    //  Average values p_a
    // Get the pressure at the integrations points
    long global_point_index = datavec[0].intGlobPtIndex;
    TRMMemory &point_memory = GetMemory()[global_point_index];
    
    REAL p_avg_n  = point_memory.p_avg_n();
    
    //  Computing closure relationship at given average values
    
    TPZManVector<STATE, 10> v(nvars),v_avg(nvars);
    v[0]        = p;
    v_avg[0]    = p_avg_n;
    
    // Fluid parameters
    TPZManVector<STATE, 10> rho,l;
    fSimulationData->AlphaProp()->Density(rho, v);
    fSimulationData->PetroPhysics()->l(l, v_avg);
    
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    TPZManVector<STATE, 10> phi;
    fSimulationData->Map()->Kappa(datavec[ub].x, K, Kinv, v);
    fSimulationData->Map()->phi(datavec[ub].x, phi, v);
    
    // Defining local variables
    TPZFNMatrix<3,STATE> lambda_K_inv_u(3,1);
    TPZManVector<STATE,3> Gravity = fSimulationData->Gravity();
    
    for (int i = 0; i < u.size(); i++) {
        STATE dot = 0.0;
        for (int j =0; j < u.size(); j++) {
            dot += Kinv(i,j)*u[j];
        }
        lambda_K_inv_u(i,0) = (1.0/l[0]) * dot;
    }
    
    
    // Integration point contribution
    STATE divu = 0.0;
    TPZFNMatrix<3,STATE> phi_u_i(3,1);
    
    int s_i;
    int v_i;
    
    if(! fSimulationData->IsCurrentStateQ()){
        
        REAL p_avg    = point_memory.p_avg();
        
        TPZManVector<STATE, 10> v(nvars);
        v[0] = p;//p_avg;
        
        // Fluid parameters
        TPZManVector<STATE, 10> rho,l;
        fSimulationData->AlphaProp()->Density(rho, v);
        
        // Rock parameters
        TPZManVector<STATE, 10> phi;
        fSimulationData->Map()->phi(datavec[ub].x, phi, v);
        
        for (int ip = 0; ip < nphip; ip++)
        {
            
            ef(ip + firstp) += -1.0 * weight * (-1.0/dt) * rho[0] * phi[0] * phi_ps(ip,0);
            
        }
        
        return;
    }
    
    
    
    for (int iu = 0; iu < nphiu; iu++)
    {
        
        v_i = datavec[ub].fVecShapeIndex[iu].first;
        s_i = datavec[ub].fVecShapeIndex[iu].second;
        
        STATE Kl_inv_dot_u = 0.0, rho_g_dot_phi_u = 0.0;
        for (int i = 0; i < u.size(); i++) {
            phi_u_i(i,0) = phi_us(s_i,0) * datavec[ub].fNormalVec(i,v_i);
            Kl_inv_dot_u += lambda_K_inv_u(i,0)*phi_u_i(i,0);
            rho_g_dot_phi_u += rho[0]*Gravity[i]*phi_u_i(i,0);
        }
        
        
        ef(iu + firstu) += weight * ( Kl_inv_dot_u - (1.0/jac_det) * (p) * div_on_master(iu,0) - rho_g_dot_phi_u);
        
    }
    
    
    TPZManVector<STATE,1> f(1,0.0);
    if(fForcingFunction)
    {
        fForcingFunction->Execute(datavec[pb].x,f);
    }
    
    
    divu = (Graduaxes(0,0) + Graduaxes(1,1) + Graduaxes(2,2));
    
    for (int ip = 0; ip < nphip; ip++)
    {
        
        ef(ip + firstp) += -1.0 * weight * (divu + (1.0/dt) * rho[0] * phi[0] - f[0]) * phi_ps(ip,0);
        
    }
    
    return;
    
}

void TRMMixedDarcy::ContributeBC_a(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int ub = 0;
    int pb = 1;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    
    int nphiu       = phi_us.Rows();
    int firstu      = 0;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    
    REAL Value = bc.Val2()(0,0);
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<STATE,2> f(1);
        TPZFMatrix<double> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavec[pb].x, time, f, gradf);
        Value = f[0];
    }
    else{
        Value = bc.Val2()(0,0);
    }
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD
        {
            STATE p_D = Value;
            for (int iu = 0; iu < nphiu; iu++)
            {
                ef(iu + firstu) += weight * p_D * phi_us(iu,0);
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN
        {
            
            for (int iu = 0; iu < nphiu; iu++)
            {
                STATE un_N = Value, un = u[0];
                ef(iu + firstu) += weight * gBigNumber * (un - un_N) * phi_us(iu,0);
                
                for (int ju = 0; ju < nphiu; ju++)
                {
                    
                    ek(iu + firstu,ju + firstu) += weight * gBigNumber * phi_us(ju,0) * phi_us(iu,0);
                }
                
            }
            
        }
            break;
            
        case 2 :    // Neumann BC  Impervious bc
        {
            
            for (int iu = 0; iu < nphiu; iu++)
            {
                STATE un = u[0];
                ef(iu + firstu) += weight * 100000.0 * gBigNumber * (un - 0.0) * phi_us(iu,0);
                
                for (int ju = 0; ju < nphiu; ju++)
                {
                    
                    ek(iu + firstu,ju + firstu) += weight * 100000.0 * gBigNumber * phi_us(ju,0) * phi_us(iu,0);
                }
                
            }
            
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            
            DebugStop();
        }
            break;
    }
    
    return;
    
}

void TRMMixedDarcy::Solution_a(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int ub = 0;
    int pb = 1;
    
    TPZManVector<REAL,3> u = datavec[ub].sol[0];
    REAL p = datavec[pb].sol[0][0];
    
    TPZFMatrix<STATE> dudx = datavec[ub].dsol[0];
    TPZFMatrix<STATE> dpdx = datavec[pb].dsol[0];
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0:
        {
            Solout[0] = p;
        }
            break;
        case 1:
        {
            Solout[0] = u[0]; // Bulk mass velocity
            Solout[1] = u[1]; // Bulk mass velocity
            Solout[2] = u[2]; // Bulk mass velocity
        }
            break;
        case 2:
        {
            Solout[0] = dudx(0,0) + dudx(1,1) + dudx(2,2);
        }
            break;
        default:
        {
            DebugStop();
        }
    }
}

// ------------------------------------------------------------------- //
// two phase flow case
// ------------------------------------------------------------------- //


void TRMMixedDarcy::Contribute_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int nvars = 4; // {p,sa,sb,t}
    
    int ub      = 0;
    int pb      = 1;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,STATE> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,STATE> dphi_us      = datavec[ub].dphix;
    TPZFNMatrix<100,STATE> dphi_ps      = datavec[pb].dphix;
    
    TPZFNMatrix<40,STATE> div_on_master;
    STATE divflux;
    this->ComputeDivergenceOnMaster(datavec, div_on_master,divflux);
    REAL jac_det = datavec[ub].detjac;
    
    int nphiu       = datavec[ub].fVecShapeIndex.NElements();
    int nphip       = phi_ps.Rows();
    int firstu      = 0;
    int firstp      = nphiu + firstu;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    REAL p                  = datavec[pb].sol[0][0];
    
    TPZFNMatrix<10,STATE> Graduaxes = datavec[ub].dsol[0];
    
    // Time
    STATE dt = fSimulationData->dt();
    
    //  Average values p_a
    
    REAL s = 0.0;
    
    REAL p_a    = p;
    REAL s_a    = s;
    
    //  Computing closure relationship at given average values
    
    TPZManVector<STATE, 10> v(nvars);
    v[0] = p_a;
    v[1] = s_a;
    
    // Fluid parameters
    TPZManVector<STATE, 10> rho_a,rho_b,l;
    fSimulationData->AlphaProp()->Density(rho_a, v);
    fSimulationData->BetaProp()->Density(rho_b, v);
    fSimulationData->PetroPhysics()->l(l, v);
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    TPZManVector<STATE, 10> phi;
    fSimulationData->Map()->Kappa(datavec[ub].x, K, Kinv, v);
    fSimulationData->Map()->phi(datavec[ub].x, phi, v);
    
    // Defining local variables
    TPZFNMatrix<3,STATE> lambda_K_inv_u(3,1),lambda_dp_K_inv_u(3,1), lambda_ds_K_inv_u(3,1), lambda_K_inv_phi_u_j(3,1);
    TPZManVector<STATE,3> Gravity = fSimulationData->Gravity();
    
    for (int i = 0; i < u.size(); i++) {
        STATE dot = 0.0;
        for (int j =0; j < u.size(); j++) {
            dot += Kinv(i,j)*u[j];
        }
        lambda_K_inv_u(i,0)     = (1.0/l[0]) * dot;
        lambda_dp_K_inv_u(i,0)  = (-l[1]/(l[0]*l[0])) * dot;
        lambda_ds_K_inv_u(i,0)  = (-l[2]/(l[0]*l[0])) * dot;
    }
    
    // Integration point contribution
    STATE divu = 0.0;
    TPZFNMatrix<3,STATE> phi_u_i(3,1), phi_u_j(3,1);
    
    int s_i, s_j;
    int v_i, v_j;
    
    if(! fSimulationData->IsCurrentStateQ()){
        for (int ip = 0; ip < nphip; ip++)
        {
            
            ef(ip + firstp) += -1.0 * weight * (-1.0/dt) * (s*rho_a[0]+(1.0-s)*rho_b[0]) * phi[0] * phi_ps(ip,0);
            
        }
        
        return;
    }
    
    for (int iu = 0; iu < nphiu; iu++)
    {
        
        v_i = datavec[ub].fVecShapeIndex[iu].first;
        s_i = datavec[ub].fVecShapeIndex[iu].second;
        
        STATE Kl_inv_dot_u = 0.0, Kl_dp_inv_dot_u = 0.0, Kl_ds_inv_dot_u = 0.0, rho_g_dot_phi_u = 0.0, rho_dp_g_dot_phi_u = 0.0, rho_ds_g_dot_phi_u = 0.0;
        for (int i = 0; i < u.size(); i++) {
            phi_u_i(i,0) = phi_us(s_i,0) * datavec[ub].fNormalVec(i,v_i);
            Kl_inv_dot_u        += lambda_K_inv_u(i,0)*phi_u_i(i,0);
            Kl_dp_inv_dot_u     += lambda_dp_K_inv_u(i,0)*phi_u_i(i,0);
            Kl_ds_inv_dot_u     += lambda_ds_K_inv_u(i,0)*phi_u_i(i,0);
            rho_g_dot_phi_u     += (s*rho_a[0]+(1.0-s)*rho_b[0])*Gravity[i]*phi_u_i(i,0);
            rho_dp_g_dot_phi_u  += (s*rho_a[1]+(1.0-s)*rho_b[1])*Gravity[i]*phi_u_i(i,0);
            rho_ds_g_dot_phi_u  += (rho_a[0]-rho_b[0])*Gravity[i]*phi_u_i(i,0);
        }
        
        ef(iu + firstu) += weight * ( Kl_inv_dot_u - (1.0/jac_det) * (p) * div_on_master(iu,0) - rho_g_dot_phi_u);
        
        for (int ju = 0; ju < nphiu; ju++)
        {
            
            v_j = datavec[ub].fVecShapeIndex[ju].first;
            s_j = datavec[ub].fVecShapeIndex[ju].second;
            
            STATE Kl_inv_phi_u_j_dot_phi_u_i = 0.0;
            for (int j = 0; j < u.size(); j++) {
                phi_u_j(j,0) = phi_us(s_j,0) * datavec[ub].fNormalVec(j,v_j);
                STATE dot = 0.0;
                for (int k = 0; k < u.size(); k++) {
                    dot += (1.0/l[0]) * Kinv(j,k)*phi_u_j(k,0);
                }
                lambda_K_inv_phi_u_j(j,0) = dot;
                Kl_inv_phi_u_j_dot_phi_u_i += lambda_K_inv_phi_u_j(j,0)*phi_u_i(j,0);
            }
            
            
            ek(iu + firstu,ju + firstu) += weight * Kl_inv_phi_u_j_dot_phi_u_i;
        }
        
        for (int jp = 0; jp < nphip; jp++)
        {
            ek(iu + firstu, jp + firstp) += weight * ( Kl_dp_inv_dot_u - (1.0/jac_det) * div_on_master(iu,0) + rho_dp_g_dot_phi_u) * phi_ps(jp,0);
        }
        
    }
    
    
    TPZManVector<STATE,1> f(1,0.0);
    if(fForcingFunction)
    {
        fForcingFunction->Execute(datavec[pb].x,f);
    }
    
    
    divu = (Graduaxes(0,0) + Graduaxes(1,1) + Graduaxes(2,2));
    
    for (int ip = 0; ip < nphip; ip++)
    {
        
        ef(ip + firstp) += -1.0 * weight * (divu + (1.0/dt) *(s*rho_a[0]+(1.0-s)*rho_b[0]) * phi[0] - f[0]) * phi_ps(ip,0);
        
        for (int ju = 0; ju < nphiu; ju++)
        {
            ek(ip + firstp, ju + firstu) += -1.0 * weight * (1.0/jac_det) * div_on_master(ju,0) * phi_ps(ip,0);
        }
        
        for (int jp = 0; jp < nphip; jp++)
        {
            ek(ip + firstp, jp + firstp) += -1.0 * weight * ( (1.0/dt) * ((s*rho_a[0]+(1.0-s)*rho_b[0]) * phi[1] + (s*rho_a[1]+(1.0-s)*rho_b[1]) * phi[0]) * phi_ps(ip,0) ) * phi_ps(jp,0);
        }
        
    }

    
}

void TRMMixedDarcy::Contribute_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    int nvars = 4; // {p,sa,sb,t}
    
    int ub      = 0;
    int pb      = 1;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,STATE> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,STATE> dphi_us      = datavec[ub].dphix;
    TPZFNMatrix<100,STATE> dphi_ps      = datavec[pb].dphix;
    
    TPZFNMatrix<40,STATE> div_on_master;
    STATE divflux;
    this->ComputeDivergenceOnMaster(datavec, div_on_master,divflux);
    REAL jac_det = datavec[ub].detjac;
    
    int nphiu       = datavec[ub].fVecShapeIndex.NElements();
    int nphip       = phi_ps.Rows();
    int firstu      = 0;
    int firstp      = nphiu + firstu;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    REAL p                  = datavec[pb].sol[0][0];
    
    TPZFNMatrix<10,STATE> Graduaxes = datavec[ub].dsol[0];
    
    // Time
    STATE dt = fSimulationData->dt();
    
    REAL s = 0.0;
    
    //  Average values p_a
    
    REAL p_a    = p;
    REAL s_a    = s;
    
    //  Computing closure relationship at given average values
    
    TPZManVector<STATE, 10> v(nvars);
    v[0] = p_a;
    v[1] = s_a;
    
    // Fluid parameters
    TPZManVector<STATE, 10> rho_a,rho_b,l;
    fSimulationData->AlphaProp()->Density(rho_a, v);
    fSimulationData->BetaProp()->Density(rho_b, v);
    fSimulationData->PetroPhysics()->l(l, v);
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    TPZManVector<STATE, 10> phi;
    fSimulationData->Map()->Kappa(datavec[ub].x, K, Kinv, v);
    fSimulationData->Map()->phi(datavec[ub].x, phi, v);
    
    // Defining local variables
    TPZFNMatrix<3,STATE> lambda_K_inv_u(3,1);
    TPZManVector<STATE,3> Gravity = fSimulationData->Gravity();
    
    for (int i = 0; i < u.size(); i++) {
        STATE dot = 0.0;
        for (int j =0; j < u.size(); j++) {
            dot += Kinv(i,j)*u[j];
        }
        lambda_K_inv_u(i,0) = (1.0/l[0]) * dot;
    }
    
    
    // Integration point contribution
    STATE divu = 0.0;
    TPZFNMatrix<3,STATE> phi_u_i(3,1);
    
    int s_i;
    int v_i;
    
    if(! fSimulationData->IsCurrentStateQ()){
        for (int ip = 0; ip < nphip; ip++)
        {
            
            ef(ip + firstp) += -1.0 * weight * (-1.0/dt) * (s*rho_a[0]+(1.0-s)*rho_b[0]) * phi[0] * phi_ps(ip,0);
            
        }
        
        
        return;
    }
    
    
    
    for (int iu = 0; iu < nphiu; iu++)
    {
        
        v_i = datavec[ub].fVecShapeIndex[iu].first;
        s_i = datavec[ub].fVecShapeIndex[iu].second;
        
        STATE Kl_inv_dot_u = 0.0, rho_g_dot_phi_u = 0.0;
        for (int i = 0; i < u.size(); i++) {
            phi_u_i(i,0) = phi_us(s_i,0) * datavec[ub].fNormalVec(i,v_i);
            Kl_inv_dot_u += lambda_K_inv_u(i,0)*phi_u_i(i,0);
            rho_g_dot_phi_u += (s*rho_a[0]+(1.0-s)*rho_b[0])*Gravity[i]*phi_u_i(i,0);
        }
        
        
        ef(iu + firstu) += weight * ( Kl_inv_dot_u - (1.0/jac_det) * (p) * div_on_master(iu,0) - rho_g_dot_phi_u);
        
    }
    
    
    TPZManVector<STATE,1> f(1,0.0);
    if(fForcingFunction)
    {
        fForcingFunction->Execute(datavec[pb].x,f);
    }
    
    
    divu = (Graduaxes(0,0) + Graduaxes(1,1) + Graduaxes(2,2));
    
    for (int ip = 0; ip < nphip; ip++)
    {
        
        ef(ip + firstp) += -1.0 * weight * (divu + (1.0/dt) * (s*rho_a[0]+(1.0-s)*rho_b[0]) * phi[0] - f[0]) * phi_ps(ip,0);
        
    }
    
    return;
    
}

void TRMMixedDarcy::ContributeBC_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int ub = 0;
    int pb = 1;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    
    int nphiu       = phi_us.Rows();
    int firstu      = 0;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    
    REAL Value_m    = 0.0;
    REAL Value_s    = 0.0;
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<STATE,2> f(2);
        TPZFMatrix<double> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavec[pb].x, time, f, gradf);
        Value_m = f[0];
        Value_s = f[1];
    }
    else{
        Value_m = bc.Val2()(0,0);
        Value_s = bc.Val2()(1,0);
    }
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD outlet
        {
            STATE p_D = Value_m;
            for (int iu = 0; iu < nphiu; iu++)
            {
                ef(iu + firstu) += weight * p_D * phi_us(iu,0);
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN outlet
        {
            
            for (int iu = 0; iu < nphiu; iu++)
            {
                STATE un_N = Value_m, un = u[0];
                ef(iu + firstu) += weight * gBigNumber * (un - un_N) * phi_us(iu,0);
                
                for (int ju = 0; ju < nphiu; ju++)
                {
                    
                    ek(iu + firstu,ju + firstu) += weight * gBigNumber * phi_us(ju,0) * phi_us(iu,0);
                }
                
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD inlet
        {
            STATE p_D = Value_m;
            for (int iu = 0; iu < nphiu; iu++)
            {
                ef(iu + firstu) += weight * p_D * phi_us(iu,0);
            }
        }
            break;
            
        case 3 :    // Neumann BC  QN inlet
        {
            
            for (int iu = 0; iu < nphiu; iu++)
            {
                STATE un_N = Value_m, un = u[0];
                ef(iu + firstu) += weight * gBigNumber * (un - un_N) * phi_us(iu,0);
                
                for (int ju = 0; ju < nphiu; ju++)
                {
                    
                    ek(iu + firstu,ju + firstu) += weight * gBigNumber * phi_us(ju,0) * phi_us(iu,0);
                }
                
            }
            
        }
            break;
            
        case 4 :    // Neumann BC  Impervious bc
        {
            
            for (int iu = 0; iu < nphiu; iu++)
            {
                STATE un = u[0];
                ef(iu + firstu) += weight * 100000.0 * gBigNumber * (un - 0.0) * phi_us(iu,0);
                
                for (int ju = 0; ju < nphiu; ju++)
                {
                    
                    ek(iu + firstu,ju + firstu) += weight * 100000.0 * gBigNumber * phi_us(ju,0) * phi_us(iu,0);
                }
                
            }
            
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            
            DebugStop();
        }
            break;
    }
    
    return;
    
}

void TRMMixedDarcy::ContributeBCInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    DebugStop();
}

void TRMMixedDarcy::ContributeBCInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    DebugStop();
}

void TRMMixedDarcy::ContributeInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    DebugStop();
    
}

void TRMMixedDarcy::ContributeInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
    
    DebugStop();
    
    
}

void TRMMixedDarcy::Solution_ab(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int ub = 0;
    int pb = 1;
    int sb = 2;
    
    TPZManVector<REAL,3> u = datavec[ub].sol[0];
    REAL p = datavec[pb].sol[0][0];
    REAL s = 0.0;//datavec[sb].sol[0][0];
    
    TPZFMatrix<STATE> dudx = datavec[ub].dsol[0];
    TPZFMatrix<STATE> dpdx = datavec[pb].dsol[0];
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0:
        {
            Solout[0] = p;
        }
            break;
        case 1:
        {
            Solout[0] = u[0]; // Bulk mass velocity
            Solout[1] = u[1]; // Bulk mass velocity
            Solout[2] = u[2]; // Bulk mass velocity
        }
            break;
        case 2:
        {
            Solout[0] = dudx(0,0) + dudx(1,1) + dudx(2,2);
        }
            break;
        case 3:
        {
            Solout[0] = s;
        }
            break;
        default:
        {
            DebugStop();
        }
    }
}


// ------------------------------------------------------------------- //
// three phase flow case
// ------------------------------------------------------------------- //


void TRMMixedDarcy::Contribute_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
}

void TRMMixedDarcy::Contribute_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
}

void TRMMixedDarcy::ContributeBC_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
}

void TRMMixedDarcy::ContributeBCInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
}

void TRMMixedDarcy::ContributeBCInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
}

void TRMMixedDarcy::ContributeInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
}

void TRMMixedDarcy::ContributeInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
    
}


int TRMMixedDarcy::ClassId() const {
    return -63786378;
}

// -------------------------------------------------------------------------------------------

void TRMMixedDarcy::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
}

// -------------------------------------------------------------------------------------------

void TRMMixedDarcy::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}

// Update element memory by copying the n+1 data to the n data
void TRMMixedDarcy::UpdateMemory()
{
    DebugStop();
}

