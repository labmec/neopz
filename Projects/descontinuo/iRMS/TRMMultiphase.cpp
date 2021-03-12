//
//  TRMMultiphase.cpp
//  PZ
//
//  Created by Omar on 5/15/16.
//
//

#include "TRMMultiphase.h"



TRMMultiphase::TRMMultiphase() : TPZMatWithMem<TRMMemory, TPZMaterial>()
{

    fdimension = 0;
    
}

TRMMultiphase::TRMMultiphase(int matid, int dimension) : TPZMatWithMem<TRMMemory, TPZMaterial>(matid)
{
    fdimension = dimension;
}


TRMMultiphase::TRMMultiphase(const TRMMultiphase &mat) : TPZMatWithMem<TRMMemory, TPZMaterial>(mat)
{
    this->fdimension = mat.fdimension;
}

TRMMultiphase::~TRMMultiphase()
{
    
}

void TRMMultiphase::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

void TRMMultiphase::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TRMMultiphase::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

int TRMMultiphase::VariableIndex(const std::string &name) {
    if (!strcmp("u", name.c_str())) return 0;
    if (!strcmp("q", name.c_str())) return 1;
    if (!strcmp("p", name.c_str())) return 2;
    if (!strcmp("div_u", name.c_str())) return 3;
    if (!strcmp("div_q", name.c_str())) return 4;
    if (!strcmp("s_xx", name.c_str())) return 5;
    if (!strcmp("s_yy", name.c_str())) return 6;
    if (!strcmp("s_xy", name.c_str())) return 7;
    if (!strcmp("s_a", name.c_str())) return 8;
    if (!strcmp("s_b", name.c_str())) return 9;
    if (!strcmp("s_c", name.c_str())) return 10;
    return TPZMatWithMem::VariableIndex(name);
}

int TRMMultiphase::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return fdimension; // Vector
        case 1:
            return fdimension; // Vector
        case 2:
            return 1; // Scalar
        case 3:
            return 1; // Scalar
        case 4:
            return 1; // Scalar
        case 5:
            return 1; // Scalar
        case 6:
            return 1; // Scalar
        case 7:
            return 1; // Scalar
        case 8:
            return 1; // Scalar
        case 9:
            return 1; // Scalar
        case 10:
            return 1; // Scalar
    }
    return TPZMatWithMem::NSolutionVariables(var);
}

void TRMMultiphase::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {

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
            Solution_abc(datavec, var, Solout);
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
void TRMMultiphase::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
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
            Contribute_abc(datavec, weight, ek, ef);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
}


void TRMMultiphase::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
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
            Contribute_abc(datavec, weight, ef);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
}

void TRMMultiphase::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            ContributeInterface_ab(data, datavecleft, datavecright, weight, ek, ef);
        }
            break;
        case 3:
        {
            ContributeInterface_abc(data, datavecleft, datavecright, weight, ek, ef);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
}

void TRMMultiphase::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef)
{
    std::cout << " This method should be called only for Capillary pressure terms " << std::endl;
    DebugStop();
}


void TRMMultiphase::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            ContributeBCInterface_ab(data, datavecleft, weight, ek, ef, bc);
        }
            break;
        case 3:
        {
            ContributeBCInterface_abc(data, datavecleft, weight, ek, ef, bc);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
}

void TRMMultiphase::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    std::cout << " This method should be called only for Capillary pressure terms " << std::endl;
    DebugStop();
}


void TRMMultiphase::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
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
            ContributeBC_abc(datavec, weight, ek, ef, bc);
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

void TRMMultiphase::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU)
{
    int qb = 1;
    int dim = this->Dimension();
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[qb].phi;   // For H1  test functions Q
    TPZFMatrix<STATE> dphiuH1       = datavec[qb].dphi; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiuH1axes   = datavec[qb].dphix; // Derivative For H1  test functions
    TPZFNMatrix<9,STATE> gradu = datavec[qb].dsol[0];
    TPZFNMatrix<9,STATE> graduMaster;
    gradu.Transpose();
    
    TPZFNMatrix<660> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[qb].axes);
    
    int nphiuHdiv = datavec[qb].fVecShapeIndex.NElements();
    
    DivergenceofPhi.Resize(nphiuHdiv,1);
    
    REAL JacobianDet = datavec[qb].detjac;
    
    TPZFMatrix<STATE> Qaxes = datavec[qb].axes;
    TPZFMatrix<STATE> QaxesT;
    TPZFMatrix<STATE> Jacobian = datavec[qb].jacobian;
    TPZFMatrix<STATE> JacobianInverse = datavec[qb].jacinv;
    
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
            ivectorindex = datavec[qb].fVecShapeIndex[iq].first;
            ishapeindex = datavec[qb].fVecShapeIndex[iq].second;
            
            for (int k = 0; k < dim; k++) {
                VectorOnXYZ(k,0) = datavec[qb].fDeformedDirections(k,ivectorindex);
            }
            
            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
            VectorOnMaster *= JacobianDet;
            
            /* Contravariant Piola mapping preserves the divergence */
            for (int k = 0; k < dim; k++) {
                DivergenceofPhi(iq,0) +=  dphiuH1(k,ishapeindex)*VectorOnMaster(k,0);
            }
        }
        
        GradOfXInverse.Multiply(gradu, graduMaster);
        graduMaster *= JacobianDet;
        for (int k = 0; k < dim; k++) {
            DivergenceofU += graduMaster(k,k);
        }

    }
    else
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[qb].fVecShapeIndex[iq].first;
            ishapeindex = datavec[qb].fVecShapeIndex[iq].second;
            
            /* Computing the divergence for constant jacobian elements */
            for (int k = 0; k < dim; k++) {
                DivergenceofPhi(iq,0) +=  datavec[qb].fDeformedDirections(k,ivectorindex)*GradphiuH1(k,ishapeindex);
            }
        }
    }
    
    return;
    
}

//** @brief Compute elastic stress */
void TRMMultiphase::Sigma(TPZManVector<STATE, 10> & l, TPZManVector<STATE, 10> & mu, TPZFMatrix<REAL> & Grad_u, TPZFMatrix<REAL> & S){
    
#ifdef PZDEBUG
    if(Grad_u.Rows() != fdimension && Grad_u.Cols() != fdimension){
        DebugStop();
    }
#endif

    TPZFMatrix<REAL> Grad_u_t = Grad_u, e;
    Grad_u.Transpose(&Grad_u_t);

    e = Grad_u + Grad_u_t;
    e *= 0.5;

    TPZFMatrix<REAL> I(fdimension,fdimension);
    I.Identity();
    
    REAL tr_e = 0.0;
    for (int i = 0; i < fdimension; i++) {
        tr_e += e(i,i);
    }
    
    S = 2.0 * mu[0] * e + l[0] * tr_e * I;

}

// ------------------------------------------------------------------- //
// one phase flow case
// ------------------------------------------------------------------- //

////** @brief Compute elastic stress */
//void TRMMultiphase::Porosity_correction(TPZManVector<STATE, 10> & l, TPZManVector<STATE, 10> & mu, TPZFMatrix<REAL> & Grad_u, TPZFMatrix<REAL> & S){
//    
//#ifdef PZDEBUG
//    if(Grad_u.Rows() != fdimension && Grad_u.Cols() != fdimension){
//        DebugStop();
//    }
//#endif
//    
//    TPZFMatrix<REAL> Grad_u_t = Grad_u, e;
//    Grad_u.Transpose(&Grad_u_t);
//    
//    e = Grad_u + Grad_u_t;
//    e *= 0.5;
//    
//    TPZFMatrix<REAL> I(fdimension,fdimension);
//    I.Identity();
//    
//    REAL tr_e = 0.0;
//    for (int i = 0; i < fdimension; i++) {
//        tr_e += e(i,i);
//    }
//    
//    S = 2.0 * mu[0] * e + l[0] * tr_e * I;
//    
//}

void TRMMultiphase::Contribute_a(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    
    int nvars = 4; // {p,sa,sb,t}
    
    int ub = 0;
    int qb = 1;
    int pb = 2;

    TPZFNMatrix<100,REAL> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].phi;
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,REAL> dphi_us      = datavec[ub].dphix;
    TPZFNMatrix<300,REAL> dphi_qs      = datavec[qb].dphix;
    TPZFNMatrix<100,REAL> dphi_ps      = datavec[pb].dphix;
    
    TPZFNMatrix<40,REAL> div_on_master;
    STATE divflux;
    this->ComputeDivergenceOnMaster(datavec, div_on_master,divflux);
    REAL jac_det = datavec[qb].detjac;

    int n_u         = fdimension;
    int nphiu       = phi_us.Rows();
    int nphiq       = datavec[qb].fVecShapeIndex.NElements();
    int nphip       = phi_ps.Rows();
    int firstu      = 0;
    int firstq      = nphiu*n_u + firstu;
    int firstp      = nphiq + firstq;

    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    TPZManVector<REAL,3> q  = datavec[qb].sol[0];
    REAL p                  = datavec[pb].sol[0][0];
    
    TPZFNMatrix <9,REAL> du         = datavec[ub].dsol[0];
    TPZFNMatrix <9,REAL> dp         = datavec[pb].dsol[0];
    TPZFNMatrix<10,STATE> Gradqaxes = datavec[qb].dsol[0];
    
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[ub].axes;
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[pb].axes;
    
    // Time
    STATE dt = fSimulationData->dt();
    
    //  Average values p_a
    REAL p_a    = p;
    
    //  Computing closure relationship at given average values
    
    TPZManVector<STATE, 10> v(nvars);
    v[0] = p_a;
    
    // Fluid parameters
    TPZManVector<STATE, 10> rho,l;
    fSimulationData->AlphaProp()->Density(rho, v);
    fSimulationData->PetroPhysics()->l(l, v);
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    TPZManVector<STATE, 10> phi, lambda, lambda_u, mu, alpha;
    fSimulationData->Map()->Kappa(datavec[qb].x, K, Kinv, v);
    fSimulationData->Map()->phi(datavec[qb].x, phi, v);
    
    // Defining local variables
    TPZFNMatrix<3,STATE> lambda_K_inv_q(3,1),lambda_dp_K_inv_q(3,1), lambda_K_inv_phi_q_j(3,1);
    TPZManVector<STATE,3> Gravity = fSimulationData->Gravity();
    TPZFMatrix<STATE> S_geo;
    
    fSimulationData->Map()->S_0(datavec[ub].x, S_geo);
    fSimulationData->Map()->lambda(datavec[ub].x, lambda, v);
    fSimulationData->Map()->lambda_u(datavec[ub].x, lambda_u, v);
    fSimulationData->Map()->mu(datavec[ub].x, mu, v);
    fSimulationData->Map()->alpha(datavec[ub].x, alpha, v);
    
    TPZFNMatrix<300,REAL> Grad_us_xy;
    TPZFNMatrix<9,REAL> axes_u_T, Gradu_xy;
    
    axes_u.Transpose(&axes_u_T);
    axes_u.Multiply(dphi_us,Grad_us_xy,1/* Transpose axes_u */);
    axes_u.Multiply(du,Gradu_xy,1/* Transpose axes_u */);
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_p(n_u,n_u,0.0);
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0); // dp/dx
    Grad_p(0,1) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1); // dp/dy
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_u(n_u,n_u,0.0),S;
    
    for (int i = 0 ; i < fdimension; i++) {
        for (int j = 0 ; j < fdimension; j++) {
            Grad_u(i,j) = Gradu_xy(i,j);
        }
    }
    Sigma(lambda, mu, Grad_u, S);
    
    TPZManVector<STATE, 10> phi_star(nvars+1);// @omar:: problems with the definition of state var now is reuired introduce u
    REAL rho_s = 2500.0;
    REAL div_u = Grad_u(0,0) + Grad_u(1,1);
    REAL Se = 1.0e-8; // average compressibility
    phi_star[0] = phi[0] + alpha[0] * div_u + Se * p;
    phi_star[1] = Se;
    
    REAL rho_solid_avg =  rho_s*(1.0-phi_star[0])+rho[0]*phi_star[0];

    TPZManVector<REAL , 3 > b(3,0.0);
    b[0] = rho_solid_avg * Gravity[0];
    b[1] = rho_solid_avg * Gravity[1];
    
    for (int i = 0; i < q.size(); i++) {
        REAL dot = 0.0;
        for (int j =0; j < q.size(); j++) {
            dot += Kinv(i,j)*q[j];
        }
        lambda_K_inv_q(i,0)     = (1.0/l[0]) * dot;
        lambda_dp_K_inv_q(i,0)  = (-l[1]/(l[0]*l[0])) * dot;
    }
    
    // Integration point contribution
    REAL divq = 0.0;
    TPZFNMatrix<3,STATE> phi_q_i(3,1), phi_q_j(3,1);
    
    int s_i, s_j;
    int v_i, v_j;
    
    if(! fSimulationData->IsCurrentStateQ()){
        for (int ip = 0; ip < nphip; ip++)
        {
            
//            ef(ip + firstp) += -1.0 * weight * (-1.0/dt) * rho[0] * phi[0] * phi_ps(ip,0);
            ef(ip + firstp) += -1.0 * weight * (-1.0/dt) * rho[0] * phi_star[0] * phi_ps(ip,0);
            
        }
        
        return;
    }
    
    
    for (int iu = 0; iu < nphiu; iu++)
    {
        
        ef(n_u*iu + 0 + firstu, 0)  += weight * (S(0,0) * Grad_us_xy(0,iu) + S(0,1) * Grad_us_xy(1,iu) - alpha[0] * p * Grad_us_xy(0,iu) - b[0] * phi_us(iu, 0));
        ef(n_u*iu + 1 + firstu, 0)	+= weight * (S(1,0) * Grad_us_xy(0,iu) + S(1,1) * Grad_us_xy(1,iu) - alpha[0] * p * Grad_us_xy(1,iu) - b[1] * phi_us(iu, 0));
        
        for (int ju = 0; ju < nphiu; ju++) {
            
            ek(n_u*iu + firstu, n_u*ju + firstu)        += weight * ( (2.0*mu[0] + lambda[0]) * Grad_us_xy(0,iu) * Grad_us_xy(0,ju) + mu[0] * Grad_us_xy(1,iu) * Grad_us_xy(1,ju) );
            ek(n_u*iu + firstu, n_u*ju+1 + firstu)      += weight * (  lambda[0] * Grad_us_xy(0,iu) * Grad_us_xy(1,ju)  + mu[0] * Grad_us_xy(1,iu) * Grad_us_xy(0,ju) );
            ek(n_u*iu+1 + firstu, n_u*ju + firstu)      += weight * (  mu[0] * Grad_us_xy(0,iu) * Grad_us_xy(1,ju) + lambda[0] * Grad_us_xy(1,iu) * Grad_us_xy(0,ju) );
            ek(n_u*iu+1 + firstu, n_u*ju+1 + firstu)	+= weight * ( (2.0*mu[0] + lambda[0]) * Grad_us_xy(1,iu) * Grad_us_xy(1,ju) + mu[0] *  Grad_us_xy(0,iu) * Grad_us_xy(0,ju) );
            
        }
        
        for (int jp = 0; jp < nphip; jp++) {
            
            ek(n_u*iu + 0 + firstu, jp + firstp)      += weight * ( -1.0 * alpha[0] * phi_ps(jp,0) ) * Grad_us_xy(0,iu);
            ek(n_u*iu + 1 + firstu, jp + firstp)      += weight * ( -1.0 * alpha[0] * phi_ps(jp,0) ) * Grad_us_xy(1,iu);
            
        }
        
    }

    for (int iq = 0; iq < nphiq; iq++)
    {
        
        v_i = datavec[qb].fVecShapeIndex[iq].first;
        s_i = datavec[qb].fVecShapeIndex[iq].second;
        
        REAL Kl_inv_dot_q = 0.0, Kl_dp_inv_dot_q = 0.0, rho_g_dot_phi_q = 0.0, rho_dp_g_dot_phi_q = 0.0;
        for (int i = 0; i < q.size(); i++) {
            phi_q_i(i,0) = phi_qs(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
            Kl_inv_dot_q        += lambda_K_inv_q(i,0)*phi_q_i(i,0);
            Kl_dp_inv_dot_q     += lambda_dp_K_inv_q(i,0)*phi_q_i(i,0);
            rho_g_dot_phi_q     += rho[0]*Gravity[i]*phi_q_i(i,0);
            rho_dp_g_dot_phi_q  += rho[1]*Gravity[i]*phi_q_i(i,0);
        }
        
        ef(iq + firstq) += weight * ( Kl_inv_dot_q - (1.0/jac_det) * (p) * div_on_master(iq,0) - rho_g_dot_phi_q);
        
        for (int jq = 0; jq < nphiq; jq++)
        {
            
            v_j = datavec[qb].fVecShapeIndex[jq].first;
            s_j = datavec[qb].fVecShapeIndex[jq].second;
            
            REAL Kl_inv_phi_q_j_dot_phi_q_j = 0.0;
            for (int j = 0; j < q.size(); j++) {
                phi_q_j(j,0) = phi_qs(s_j,0) * datavec[qb].fDeformedDirections(j,v_j);
                STATE dot = 0.0;
                for (int k = 0; k < q.size(); k++) {
                    dot += (1.0/l[0]) * Kinv(j,k)*phi_q_j(k,0);
                }
                lambda_K_inv_phi_q_j(j,0) = dot;
                Kl_inv_phi_q_j_dot_phi_q_j += lambda_K_inv_phi_q_j(j,0)*phi_q_i(j,0);
            }
        
            
            ek(iq + firstq,jq + firstq) += weight * Kl_inv_phi_q_j_dot_phi_q_j;
        }
        
        for (int jp = 0; jp < nphip; jp++)
        {
            ek(iq + firstq, jp + firstp) += weight * ( Kl_dp_inv_dot_q - (1.0/jac_det) * div_on_master(iq,0) + rho_dp_g_dot_phi_q) * phi_ps(jp,0);
        }
        
    }
    
    
    TPZManVector<STATE,1> f(1,0.0);
    if(fForcingFunction)
    {
        fForcingFunction->Execute(datavec[pb].x,f);
    }
    
    divq = (Gradqaxes(0,0) + Gradqaxes(1,1) + Gradqaxes(2,2));
    
    for (int ip = 0; ip < nphip; ip++)
    {
        
//        ef(ip + firstp) += -1.0 * weight * (divq + (1.0/dt) * rho[0] * phi[0] - f[0]) * phi_ps(ip,0);
        ef(ip + firstp) += -1.0 * weight * (divq + (1.0/dt) * rho[0] * phi_star[0] - f[0]) * phi_ps(ip,0);

        
        for (int ju = 0; ju < nphiu; ju++)
        {
            
            ek(ip + firstp, n_u*ju + 0 + firstu) += -1.0 * weight * ( (1.0/dt) * alpha[0] * Grad_us_xy(0,ju) ) * phi_ps(ip,0);
            ek(ip + firstp, n_u*ju + 1 + firstu) += -1.0 * weight * ( (1.0/dt) * alpha[0] * Grad_us_xy(1,ju) ) * phi_ps(ip,0);
        }
        
        for (int jq = 0; jq < nphiq; jq++)
        {
            ek(ip + firstp, jq + firstq) += -1.0 * weight * (1.0/jac_det) * div_on_master(jq,0) * phi_ps(ip,0);
        }
        
        for (int jp = 0; jp < nphip; jp++)
        {
//            ek(ip + firstp, jp + firstp) += -1.0 * weight * ( (1.0/dt) * (rho[0] * phi[1] + rho[1] * phi[0]) * phi_ps(ip,0) ) * phi_ps(jp,0);
            ek(ip + firstp, jp + firstp) += -1.0 * weight * ( (1.0/dt) * (rho[0] * phi_star[1] + rho[1] * phi_star[0]) * phi_ps(ip,0) ) * phi_ps(jp,0);

        }
        
    }
    
}

void TRMMultiphase::Contribute_a(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    
    TPZFMatrix<STATE> ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute_a(datavec, weight, ek_fake, ef);
    return;
    
    int nvars = 4; // {p,sa,sb,t}

    int ub = 0;
    int qb = 1;
    int pb = 2;

    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,STATE> phi_qs       = datavec[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,STATE> dphi_us      = datavec[qb].dphix;
    TPZFNMatrix<100,STATE> dphi_ps      = datavec[pb].dphix;
    
    TPZFNMatrix<40,STATE> div_on_master;
    STATE divflux;
    this->ComputeDivergenceOnMaster(datavec, div_on_master,divflux);
    REAL jac_det = datavec[qb].detjac;
    
    int n_u         = fdimension;
    int nphiu       = phi_us.Rows();
    int nphiq       = datavec[qb].fVecShapeIndex.NElements();
    int nphip       = phi_ps.Rows();
    int firstu      = 0;
    int firstq      = nphiu*n_u + firstu;
    int firstp      = nphiq + firstq;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    TPZManVector<REAL,3> q  = datavec[qb].sol[0];
    REAL p                  = datavec[pb].sol[0][0];

    TPZFNMatrix <9,REAL> du         = datavec[ub].dsol[0];
    TPZFNMatrix<10,STATE> Graduaxes = datavec[qb].dsol[0];
    
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[ub].axes;

    // Time
    STATE dt = fSimulationData->dt();
    
    //  Average values p_a
    
    REAL p_a    = p;
    
    //  Computing closure relationship at given average values

    TPZManVector<STATE, 10> v(nvars);
    v[0] = p_a;

    // Fluid parameters
    TPZManVector<STATE, 10> rho,l;
    fSimulationData->AlphaProp()->Density(rho, v);
    fSimulationData->PetroPhysics()->l(l, v);
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    TPZManVector<STATE, 10> phi, lambda, lambda_u, mu, alpha;
    fSimulationData->Map()->Kappa(datavec[qb].x, K, Kinv, v);
    fSimulationData->Map()->phi(datavec[qb].x, phi, v);
    
    fSimulationData->Map()->lambda(datavec[ub].x, lambda, v);
    fSimulationData->Map()->lambda_u(datavec[ub].x, lambda_u, v);
    fSimulationData->Map()->mu(datavec[ub].x, mu, v);
    fSimulationData->Map()->alpha(datavec[ub].x, alpha, v);
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_u(n_u,n_u,0.0),S;
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    
    Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    Sigma(lambda, mu, Grad_u, S);
    
    // Defining local variables
    TPZFNMatrix<3,STATE> lambda_K_inv_q(3,1);
    TPZManVector<STATE,3> Gravity = fSimulationData->Gravity();
    
    for (int i = 0; i < q.size(); i++) {
        STATE dot = 0.0;
        for (int j =0; j < q.size(); j++) {
             dot += Kinv(i,j)*q[j];
        }
        lambda_K_inv_q(i,0) = (1.0/l[0]) * dot;
    }
    
    
    // Integration point contribution
    STATE divq = 0.0;
    TPZFNMatrix<3,STATE> phi_q_i(3,1);
    
    int s_i;
    int v_i;
    
    if(! fSimulationData->IsCurrentStateQ()){
        for (int ip = 0; ip < nphip; ip++)
        {
            
            ef(ip + firstp) += -1.0 * weight * (-1.0/dt) * rho[0] * phi[0] * phi_ps(ip,0);
            
        }
        
        return;
    }
    
    TPZFNMatrix<6,REAL> Grad_vx_i(2,1,0.0);
    TPZFNMatrix<6,REAL> Grad_vy_i(2,1,0.0);
    
    TPZFNMatrix<6,REAL> Grad_v(2,2,0.0);
    TPZFNMatrix<6,REAL> Grad_vx_j(2,1,0.0);
    TPZFNMatrix<6,REAL> Grad_vy_j(2,1,0.0);
    
    for (int iu = 0; iu < nphiu; iu++)
    {
        
        // Computing Gradient of the test function for each component
        Grad_vx_i(0,0) = dphi_us(0,iu)*axes_u(0,0)+dphi_us(1,iu)*axes_u(1,0); // dvx/dx
        Grad_vx_i(1,0) = dphi_us(0,iu)*axes_u(0,1)+dphi_us(1,iu)*axes_u(1,1); // dvx/dy
        
        Grad_vy_i(0,0) = dphi_us(0,iu)*axes_u(0,0)+dphi_us(1,iu)*axes_u(1,0); // dvy/dx
        Grad_vy_i(1,0) = dphi_us(0,iu)*axes_u(0,1)+dphi_us(1,iu)*axes_u(1,1); // dvy/dy
        
        ef(n_u*iu + firstu, 0)      += weight * (S(0,0) * Grad_vx_i(0,0) + S(0,1) * Grad_vx_i(1,0) + 0.0 * phi_us(iu, 0));
        ef(n_u*iu+1 + firstu, 0)	+= weight * (S(1,0) * Grad_vy_i(0,0) + S(1,1) * Grad_vy_i(1,0) + 0.0 * phi_us(iu, 0));
        
    }
    
    
    for (int iq = 0; iq < nphiq; iq++)
    {
        
        v_i = datavec[qb].fVecShapeIndex[iq].first;
        s_i = datavec[qb].fVecShapeIndex[iq].second;
        
        STATE Kl_inv_dot_q = 0.0, rho_g_dot_phi_q = 0.0;
        for (int i = 0; i < q.size(); i++) {
            phi_q_i(i,0) = phi_qs(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
            Kl_inv_dot_q += lambda_K_inv_q(i,0)*phi_q_i(i,0);
            rho_g_dot_phi_q += rho[0]*Gravity[i]*phi_q_i(i,0);
        }
        
        
        ef(iq + firstq) += weight * ( Kl_inv_dot_q - (1.0/jac_det) * (p) * div_on_master(iq,0) - rho_g_dot_phi_q);
        
    }
    
    
    TPZManVector<STATE,1> f(1,0.0);
    if(fForcingFunction)
    {
        fForcingFunction->Execute(datavec[pb].x,f);
    }
    
    
    divq = (Graduaxes(0,0) + Graduaxes(1,1) + Graduaxes(2,2));
    
    for (int ip = 0; ip < nphip; ip++)
    {
        
        ef(ip + firstp) += -1.0 * weight * (divq + (1.0/dt) * rho[0] * phi[0] - f[0]) * phi_ps(ip,0);

    }
    
    return;
    
}

void TRMMultiphase::ContributeBC_a(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    

    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    switch (bc.Type()) {
        case 0 :
        {
            this->apply_ux(datavec, weight, ek, ef, bc);
            this->apply_uy(datavec, weight, ek, ef, bc);
//            this->apply_uz(datavec, weight, ek, ef, bc);
            this->apply_p(datavec, weight, ek, ef, bc);
        }
            break;
        case 1 :
        {
            this->apply_ux(datavec, weight, ek, ef, bc);
            this->apply_uy(datavec, weight, ek, ef, bc);
            //            this->apply_uz(datavec, weight, ek, ef, bc);
            this->apply_q(datavec, weight, ek, ef, bc);
        }
            break;
        case 2 :
        {
            this->apply_ux(datavec, weight, ek, ef, bc);
            this->apply_p(datavec, weight, ek, ef, bc);
        }
            break;
        case 3 :
        {
            this->apply_ux(datavec, weight, ek, ef, bc);
            this->apply_q(datavec, weight, ek, ef, bc);
        }
            break;
        case 4 :
        {
            this->apply_uy(datavec, weight, ek, ef, bc);
            this->apply_p(datavec, weight, ek, ef, bc);
        }
            break;
        case 5 :
        {
            this->apply_uy(datavec, weight, ek, ef, bc);
            this->apply_q(datavec, weight, ek, ef, bc);
        }
            break;
        case 6 :
        {
            this->apply_uz(datavec, weight, ek, ef, bc);
            this->apply_p(datavec, weight, ek, ef, bc);
        }
            break;
        case 7 :
        {
            this->apply_uz(datavec, weight, ek, ef, bc);
            this->apply_q(datavec, weight, ek, ef, bc);
        }
            break;
        case 8 :
        {
            this->apply_tn(datavec, weight, ek, ef, bc);
            this->apply_p(datavec, weight, ek, ef, bc);
        }
            break;
        case 9 :
        {
            this->apply_tn(datavec, weight, ek, ef, bc);
            this->apply_q(datavec, weight, ek, ef, bc);
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

void TRMMultiphase::Solution_a(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {

    int ub = 0;
    int qb = 1;
    int pb = 2;
    
    TPZManVector<REAL,3> u = datavec[ub].sol[0];
    TPZManVector<REAL,3> q = datavec[qb].sol[0];
    REAL p = datavec[pb].sol[0][0];
    
    TPZFMatrix<STATE> dqdx = datavec[qb].dsol[0];
    TPZFMatrix<STATE> dpdx = datavec[pb].dsol[0];
    
    TPZFNMatrix <9,REAL> du         = datavec[ub].dsol[0];
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[ub].axes;
    
    Solout.Resize(this->NSolutionVariables(var));

    int nvars = 4;
    int n_u = fdimension;
    
    //  Computing closure relationship at given average values

    TPZManVector<STATE, 10> v(nvars);
    v[0] = p;
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    TPZManVector<STATE, 10> phi, lambda, lambda_u, mu, alpha;
    fSimulationData->Map()->Kappa(datavec[qb].x, K, Kinv, v);
    fSimulationData->Map()->phi(datavec[qb].x, phi, v);
    
    // Defining local variables
    TPZFNMatrix<3,STATE> lambda_K_inv_q(3,1),lambda_dp_K_inv_q(3,1), lambda_K_inv_phi_q_j(3,1);
    TPZManVector<STATE,3> Gravity = fSimulationData->Gravity();
    
    fSimulationData->Map()->lambda(datavec[ub].x, lambda, v);
    fSimulationData->Map()->lambda_u(datavec[ub].x, lambda_u, v);
    fSimulationData->Map()->mu(datavec[ub].x, mu, v);
    fSimulationData->Map()->alpha(datavec[ub].x, alpha, v);
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_u(n_u,n_u,0.0),S;
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    
    Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    Sigma(lambda, mu, Grad_u, S);
    
    // convertions
    REAL to_MPa = 1.0e-6;
    
    switch(var) {
        case 0:
        {
            // Displacement
            for (int i = 0; i < fdimension; i++) {
                Solout[i] = u[i];
            }

        }
            break;
        case 1:
        {
            // Bulk mass velocity
            for (int i = 0; i < fdimension; i++) {
                Solout[i] = q[i];
            }
        }
            break;
        case 2:
        {
            Solout[0] = p * to_MPa;
        }
            break;
        case 3:
        {
            // div_u
            for (int i = 0; i < fdimension; i++) {
                Solout[0] += Grad_u(i,i);
            }
        }
            break;
        case 4:
        {
            // div_q
            for (int i = 0; i < fdimension; i++) {
                Solout[0] += dqdx(i,i);
            }
        }
            break;
        case 5:
        {
            Solout[0] = S(0,0) * to_MPa;
        }
            break;
        case 6:
        {
            Solout[0] = S(1,1) * to_MPa;
        }
            break;
        case 7:
        {
            Solout[0] = S(1,0) * to_MPa;
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


void TRMMultiphase::Contribute_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int nvars = 4; // {p,sa,sb,t}
    
    int qb      = 0;
    int pb      = 1;
    int sb_a    = 2;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<100,STATE> phi_ss       = datavec[sb_a].phi;
    TPZFNMatrix<300,STATE> dphi_us      = datavec[qb].dphix;
    TPZFNMatrix<100,STATE> dphi_ps      = datavec[pb].dphix;
    
    TPZFNMatrix<40,STATE> div_on_master;
    STATE divflux;
    this->ComputeDivergenceOnMaster(datavec, div_on_master,divflux);
    REAL jac_det = datavec[qb].detjac;
    
    int nphiu       = datavec[qb].fVecShapeIndex.NElements();
    int nphip       = phi_ps.Rows();
    int nphis_a     = phi_ss.Rows();
    int firstu      = 0;
    int firstp      = nphiu + firstu;
    int firsts_a    = nphip + firstp;
    
    TPZManVector<REAL,3> u  = datavec[qb].sol[0];
    REAL p                  = datavec[pb].sol[0][0];
    REAL s                  = datavec[sb_a].sol[0][0];
    
    TPZFNMatrix<10,STATE> Graduaxes = datavec[qb].dsol[0];
    
    // Time
    STATE dt = fSimulationData->dt();
    
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
    fSimulationData->Map()->Kappa(datavec[qb].x, K, Kinv, v);
    fSimulationData->Map()->phi(datavec[qb].x, phi, v);
    
    // Defining local variables
    TPZFNMatrix<3,REAL> lambda_K_inv_u(3,1),lambda_dp_K_inv_u(3,1), lambda_ds_K_inv_u(3,1), lambda_K_inv_phi_u_j(3,1);
    TPZManVector<REAL,3> Gravity = fSimulationData->Gravity();
    
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
        
        for (int is = 0; is < nphis_a; is++)
        {
            
            ef(is + firsts_a) += weight * (-1.0/dt) * s * rho_a[0] * phi[0] * phi_ss(is,0);
            
        }
        
        return;
    }
    
    for (int iu = 0; iu < nphiu; iu++)
    {
        
        v_i = datavec[qb].fVecShapeIndex[iu].first;
        s_i = datavec[qb].fVecShapeIndex[iu].second;
        
        REAL Kl_inv_dot_u = 0.0, Kl_dp_inv_dot_u = 0.0, Kl_ds_inv_dot_u = 0.0, rho_g_dot_phi_u = 0.0, rho_dp_g_dot_phi_u = 0.0, rho_ds_g_dot_phi_u = 0.0;
        for (int i = 0; i < u.size(); i++) {
            phi_u_i(i,0) = phi_us(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
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
            
            v_j = datavec[qb].fVecShapeIndex[ju].first;
            s_j = datavec[qb].fVecShapeIndex[ju].second;
            
            REAL Kl_inv_phi_u_j_dot_phi_u_i = 0.0;
            for (int j = 0; j < u.size(); j++) {
                phi_u_j(j,0) = phi_us(s_j,0) * datavec[qb].fDeformedDirections(j,v_j);
                REAL dot = 0.0;
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
        
        for (int js = 0; js < nphis_a; js++)
        {
            ek(iu + firstu, js + firsts_a) += weight * ( Kl_ds_inv_dot_u + rho_ds_g_dot_phi_u) * phi_ss(js,0);
        }
        
    }
    
    
    TPZManVector<REAL,1> f(1,0.0);
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
        
        for (int js = 0; js < nphis_a; js++)
        {
            ek(ip + firstp, js + firsts_a) += -1.0 * weight * ( (1.0/dt) * ((rho_a[0]-rho_b[0]) * phi[0]) * phi_ps(ip,0) ) * phi_ss(js,0);
        }
        
    }
    
    
    for (int is = 0; is < nphis_a; is++)
    {
        
        ef(is + firsts_a) += weight * (1.0/dt) * s * rho_a[0] * phi[0] * phi_ss(is,0);
        
        for (int jp = 0; jp < nphip; jp++)
        {
            ek(is + firsts_a, jp + firstp) += weight * ( (1.0/dt) * ( s*rho_a[0] * phi[1] + s*rho_a[1] * phi[0]) * phi_ss(is,0) ) * phi_ps(jp,0);
        }
        
        for (int js = 0; js < nphis_a; js++)
        {
            ek(is + firsts_a, js + firsts_a) += weight * (1.0/dt) * rho_a[0] * phi[0] * phi_ss(js,0) * phi_ss(is,0);
        }
        
    }
    
}

void TRMMultiphase::Contribute_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    int nvars = 4; // {p,sa,sb,t}
    
    int qb      = 0;
    int pb      = 1;
    int sb_a    = 2;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<100,STATE> phi_ss       = datavec[sb_a].phi;
    TPZFNMatrix<300,STATE> dphi_us      = datavec[qb].dphix;
    TPZFNMatrix<100,STATE> dphi_ps      = datavec[pb].dphix;
    
    TPZFNMatrix<40,STATE> div_on_master;
    STATE divflux;
    this->ComputeDivergenceOnMaster(datavec, div_on_master,divflux);
    REAL jac_det = datavec[qb].detjac;
    
    int nphiu       = datavec[qb].fVecShapeIndex.NElements();
    int nphip       = phi_ps.Rows();
    int nphis_a     = phi_ss.Rows();
    int firstu      = 0;
    int firstp      = nphiu + firstu;
    int firsts_a    = nphip + firstp;
    
    TPZManVector<REAL,3> u  = datavec[qb].sol[0];
    REAL p                  = datavec[pb].sol[0][0];
    REAL s                  = datavec[sb_a].sol[0][0];
    
    TPZFNMatrix<10,STATE> Graduaxes = datavec[qb].dsol[0];
    
    // Time
    STATE dt = fSimulationData->dt();
    
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
    fSimulationData->Map()->Kappa(datavec[qb].x, K, Kinv, v);
    fSimulationData->Map()->phi(datavec[qb].x, phi, v);
    
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
        
        for (int is = 0; is < nphis_a; is++)
        {
            
            ef(is + firsts_a) += weight * (-1.0/dt) * s * rho_a[0] * phi[0] * phi_ss(is,0);
            
        }
        
        return;
    }
    
    
    
    for (int iu = 0; iu < nphiu; iu++)
    {
        
        v_i = datavec[qb].fVecShapeIndex[iu].first;
        s_i = datavec[qb].fVecShapeIndex[iu].second;
        
        STATE Kl_inv_dot_u = 0.0, rho_g_dot_phi_u = 0.0;
        for (int i = 0; i < u.size(); i++) {
            phi_u_i(i,0) = phi_us(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
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
    
    for (int is = 0; is < nphis_a; is++)
    {
        
        ef(is + firsts_a) += weight * (1.0/dt) * s * rho_a[0] * phi[0] * phi_ss(is,0);
        
    }
    
    return;
    
}

void TRMMultiphase::ContributeBC_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int qb = 0;
    int pb = 1;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[qb].phi;
    
    int nphiu       = phi_us.Rows();
    int firstu      = 0;
    
    TPZManVector<REAL,3> u  = datavec[qb].sol[0];
    
    REAL Value_m    = 0.0;
    REAL Value_s    = 0.0;
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<REAL,2> f(2);
        TPZFMatrix<REAL> gradf;
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

void TRMMultiphase::ContributeBCInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    
    int qb      = 0;
    int pb      = 1;
    int sb_a    = 2;
    
    TPZFNMatrix<100,STATE> phi_us_l       = datavecleft[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps_l       = datavecleft[pb].phi;
    TPZFNMatrix<100,STATE> phi_ss_l       = datavecleft[sb_a].phi;
    TPZFNMatrix<300,STATE> dphi_us_l      = datavecleft[qb].dphix;
    TPZFNMatrix<100,STATE> dphi_ps_l      = datavecleft[pb].dphix;
    
    int nphiu_l       = datavecleft[qb].fVecShapeIndex.NElements();
    int nphip_l       = phi_ps_l.Rows();
    int nphis_a_l     = phi_ss_l.Rows();
    int firstu_l      = 0;
    int firstp_l      = nphiu_l + firstu_l;
    int firsts_a_l    = nphip_l + firstp_l;
    
    TPZManVector<STATE,3> n = data.normal;
    TPZManVector<REAL,3> u_l  = datavecleft[qb].sol[0];
    REAL p_l                  = datavecleft[pb].sol[0][0];
    REAL s_l                  = datavecleft[sb_a].sol[0][0];
    
    STATE un_l = 0.0;
    
    for (int i = 0; i < u_l.size(); i++) {
        un_l += u_l[i]*n[i];
    }
    
    //  Average values p_a
    
    STATE p_a_l    = p_l;
    STATE s_a_l    = s_l;
    
    STATE beta = 0.0;

    
    TPZManVector<STATE, 10> fa_l,v_l(nvars+1);
    
    TPZFNMatrix<3,STATE> phi_u_i_l(3,1);
    STATE phi_un_l = 0.0;
    int s_j;
    int v_j;
    
    REAL Value_m    = 0.0;
    REAL Value_s    = 0.0;
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<REAL,2> f(2);
        TPZFMatrix<REAL> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavecleft[qb].x, time, f, gradf);
        Value_m = f[0];
        Value_s = f[1];
    }
    else{
        Value_m = bc.Val2()(0,0);
    }
    
    switch (bc.Type()) {
            
        case 0 :    // Dirichlet BC  PD outlet
        {
            
            // upwinding
            if (un_l > 0) {
                beta = 1.0;
            }

            STATE p_D = Value_m;
            
            v_l[0] = p_D;
            v_l[1] = s_a_l;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0])*phi_ss_l(is,0)*un_l;
                
//                for (int ju = 0; ju < nphiu_l; ju++) {
//                    
//                    v_j = datavecleft[qb].fVecShapeIndex[ju].first;
//                    s_j = datavecleft[qb].fVecShapeIndex[ju].second;
//                    
//                    for (int j = 0; j < u_l.size(); j++) {
//                        phi_u_i_l(j,0) = phi_us_l(s_j,0) * datavecleft[qb].fDeformedDirections(j,v_j);
//                        phi_un_l += phi_u_i_l(j,0)*n[j];
//                    }
//                    
//                    ek(is + firsts_a_l, ju + firstu_l) += +1.0*weight * beta*fa_l[0] * phi_ss_l(is,0)*phi_un_l;
//                }
                
                for (int js = 0; js < nphis_a_l; js++) {
                    ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ss_l(js,0) * phi_ss_l(is,0)*un_l;
                }
                
                
            }
            
        }
            break;
            
        case 1 :    // Neumann BC  QN outlet
        {
            
            // upwinding
            if (Value_m > 0) {
                beta = 1.0;
            }
            
            STATE un_N = Value_m;
            
            v_l[0] = p_a_l;
            v_l[1] = s_a_l;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ss_l(is,0)*un_N;
                
                
//                for (int jp = 0; jp < nphip_l; jp++) {
//                    ek(is + firsts_a_l, jp + firstp_l) += +1.0*weight * beta * fa_l[1] * phi_ps_l(jp,0) * phi_ss_l(is,0)*un_N;
//                }
                
                
                for (int js = 0; js < nphis_a_l; js++) {
                    ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ss_l(js,0) * phi_ss_l(is,0)*un_N;
                }
                
                
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD inlet
        {
            
            // upwinding
            if (un_l > 0) {
                beta = 1.0;
            }
            
            STATE p_D = Value_m;
            
            v_l[0] = p_D;
            v_l[1] = Value_s;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);

            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ss_l(is,0)*un_l;
                
//                for (int ju = 0; ju < nphiu_l; ju++) {
//                    
//                    v_j = datavecleft[qb].fVecShapeIndex[ju].first;
//                    s_j = datavecleft[qb].fVecShapeIndex[ju].second;
//                    
//                    for (int j = 0; j < u_l.size(); j++) {
//                        phi_u_i_l(j,0) = phi_us_l(s_j,0) * datavecleft[qb].fDeformedDirections(j,v_j);
//                        phi_un_l += phi_u_i_l(j,0)*n[j];
//                    }
//                    
//                    ek(is + firsts_a_l, ju + firstu_l) += +1.0*weight * beta*fa_l[0] * phi_ss_l(is,0)*phi_un_l;
//                }
    
                
            }
            
        }
            break;
            
        case 3 :    // Neumann BC  QN inlet
        {
            
            // upwinding
            if (Value_m < 0) {
                beta = 1.0;
            }
            
            STATE un_N = Value_m;
            
            v_l[0] = p_a_l;
            v_l[1] = Value_s;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);

            for (int is = 0; is < nphis_a_l; is++) {

                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ss_l(is,0)*un_N;
                
//                for (int jp = 0; jp < nphip_l; jp++) {
//                    ek(is + firsts_a_l, jp + firstp_l) += +1.0*weight * beta * fa_l[1] * phi_ps_l(jp,0) * phi_ss_l(is,0)*un_N;
//                }
                
            }
            
        }
            break;
            
        case 4 :    // Neumann BC  Impervious bc
        {
            
            // upwinding
            beta = 1.0;
            
            STATE un_N = 0.0;
            
            v_l[0] = p_a_l;
            v_l[1] = Value_s;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ss_l(is,0)*un_N;
                
//                for (int jp = 0; jp < nphip_l; jp++) {
//                    ek(is + firsts_a_l, jp + firstp_l) += +1.0*weight * beta * fa_l[1] * phi_ps_l(jp,0) * phi_ss_l(is,0)*un_N;
//                }
                
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

void TRMMultiphase::ContributeBCInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    DebugStop();    
}

void TRMMultiphase::ContributeInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    
    int qb      = 0;
    int pb      = 1;
    int sb_a    = 2;
    
    TPZFNMatrix<100,STATE> phi_us_l       = datavecleft[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps_l       = datavecleft[pb].phi;
    TPZFNMatrix<100,STATE> phi_ss_l       = datavecleft[sb_a].phi;
    TPZFNMatrix<300,STATE> dphi_us_l      = datavecleft[qb].dphix;
    TPZFNMatrix<100,STATE> dphi_ps_l      = datavecleft[pb].dphix;
    
    TPZFNMatrix<100,STATE> phi_us_r       = datavecright[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps_r       = datavecright[pb].phi;
    TPZFNMatrix<100,STATE> phi_ss_r       = datavecright[sb_a].phi;
    TPZFNMatrix<300,STATE> dphi_us_r      = datavecright[qb].dphix;
    TPZFNMatrix<100,STATE> dphi_ps_r      = datavecright[pb].dphix;
    
    int nphiu_l       = datavecleft[qb].fVecShapeIndex.NElements();
    int nphip_l       = phi_ps_l.Rows();
    int nphis_a_l     = phi_ss_l.Rows();
    int firstu_l      = 0;
    int firstp_l      = nphiu_l + firstu_l;
    int firsts_a_l    = nphip_l + firstp_l;
    
    int nphiu_r       = datavecright[qb].fVecShapeIndex.NElements();
    int nphip_r       = phi_ps_r.Rows();
    int nphis_a_r     = phi_ss_r.Rows();
    int firstu_r      = firsts_a_l + nphis_a_l;
    int firstp_r      = nphiu_r + firstu_r;
    int firsts_a_r    = nphip_r + firstp_r;
    
    TPZManVector<STATE,3> n = data.normal;
    TPZManVector<REAL,3> u_l  = datavecleft[qb].sol[0];
    REAL p_l                  = datavecleft[pb].sol[0][0];
    REAL s_l                  = datavecleft[sb_a].sol[0][0];
    
    TPZManVector<REAL,3> u_r  = datavecright[qb].sol[0];
    REAL p_r                  = datavecright[pb].sol[0][0];
    REAL s_r                  = datavecright[sb_a].sol[0][0];
    
    STATE un_l = 0.0, un_r = 0.0;
    
    for (int i = 0; i < u_l.size(); i++) {
        un_l += u_l[i]*n[i];
        un_r += u_r[i]*n[i];
    }
    
    
    //  Average values p_a
    
    STATE p_a_l    = p_l;
    STATE s_a_l    = s_l;
    STATE p_a_r    = p_r;
    STATE s_a_r    = s_r;
    
    STATE beta = 0.0;
    // upwinding
    if (un_l > 0.0) {
        beta = 1.0;
    }
    
    TPZManVector<STATE, 10> fa_l,v_l(nvars+1),fa_r,v_r(nvars+1);
    v_l[0] = p_a_l;
    v_l[1] = s_a_l;
    v_r[0] = p_a_r;
    v_r[1] = s_a_r;
    
    this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
    this->fSimulationData->PetroPhysics()->fa(fa_r, v_r);
    
    TPZFNMatrix<3,STATE> phi_u_i_l(3,1);
//    STATE phi_un_l = 0.0;
//    int s_j;
//    int v_j;
    
    for (int is = 0; is < nphis_a_l; is++) {
        
        ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ss_l(is,0)*un_l;
        
//        for (int ju = 0; ju < nphiu_l; ju++) {
//            
//            v_j = datavecleft[qb].fVecShapeIndex[ju].first;
//            s_j = datavecleft[qb].fVecShapeIndex[ju].second;
//            
//            for (int j = 0; j < u_l.size(); j++) {
//                phi_u_i_l(j,0) = phi_us_l(s_j,0) * datavecleft[qb].fDeformedDirections(j,v_j);
//                phi_un_l += phi_u_i_l(j,0)*n[j];
//            }
//            
//            ek(is + firsts_a_l, ju + firstu_l) += +1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0]) * phi_ss_l(is,0)*phi_un_l;
//        }
//        
//        for (int jp = 0; jp < nphip_l; jp++) {
//            ek(is + firsts_a_l, jp + firstp_l) += +1.0*weight * beta * fa_l[1] * phi_ps_l(jp,0) * phi_ss_l(is,0)*un_l;
//        }
//        
//        for (int jp = 0; jp < nphip_r; jp++) {
//            ek(is + firsts_a_l, jp + firstp_r) += +1.0*weight * (1.0-beta) * fa_r[1] * phi_ps_r(jp,0) * phi_ss_l(is,0)*un_l;
//        }
        
        for (int js = 0; js < nphis_a_l; js++) {
            ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ss_l(js,0) * phi_ss_l(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_a_l, js + firsts_a_r) += +1.0*weight * (1.0-beta) * fa_r[2] * phi_ss_r(js,0) * phi_ss_l(is,0)*un_l;
        }
        
    }
    
    for (int is = 0; is < nphis_a_r; is++) {
        
        ef(is + firsts_a_r) += -1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ss_r(is,0)*un_l;
        
//        for (int ju = 0; ju < nphiu_l; ju++) {
//            
//            v_j = datavecleft[qb].fVecShapeIndex[ju].first;
//            s_j = datavecleft[qb].fVecShapeIndex[ju].second;
//            
//            for (int j = 0; j < u_l.size(); j++) {
//                phi_u_i_l(j,0) = phi_us_l(s_j,0) * datavecleft[qb].fDeformedDirections(j,v_j);
//                phi_un_l += phi_u_i_l(j,0)*n[j];
//            }
//
//            ek(is + firsts_a_r, ju + firstu_l) += -1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ss_r(is,0)*phi_un_l;
//        }
//        
//        for (int jp = 0; jp < nphip_l; jp++) {
//            ek(is + firsts_a_r, jp + firstp_l) += -1.0*weight * beta * fa_l[1] * phi_ps_l(jp,0) * phi_ss_r(is,0)*un_l;
//        }
//        
//        for (int jp = 0; jp < nphip_r; jp++) {
//            ek(is + firsts_a_r, jp + firstp_r) += -1.0*weight * (1.0-beta) * fa_r[1] * phi_ps_r(jp,0) * phi_ss_r(is,0)*un_l;
//        }
        
        for (int js = 0; js < nphis_a_l; js++) {
            ek(is + firsts_a_r, js + firsts_a_l) += -1.0*weight * beta * fa_l[2] * phi_ss_l(js,0) * phi_ss_r(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_a_r, js + firsts_a_r) += -1.0*weight * (1.0-beta) * fa_r[2] * phi_ss_r(js,0) * phi_ss_r(is,0)*un_l;
        }
        
    }
    
}

void TRMMultiphase::ContributeInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    
    int qb      = 0;
    int pb      = 1;
    int sb_a    = 2;
    
    TPZFNMatrix<100,STATE> phi_us_l       = datavecleft[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps_l       = datavecleft[pb].phi;
    TPZFNMatrix<100,STATE> phi_ss_l       = datavecleft[sb_a].phi;
    
    TPZFNMatrix<100,STATE> phi_us_r       = datavecright[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps_r       = datavecright[pb].phi;
    TPZFNMatrix<100,STATE> phi_ss_r       = datavecright[sb_a].phi;
    
    int nphiu_l       = datavecleft[qb].fVecShapeIndex.NElements();
    int nphip_l       = phi_ps_l.Rows();
    int nphis_a_l     = phi_ss_l.Rows();
    int firstu_l      = 0;
    int firstp_l      = nphiu_l + firstu_l;
    int firsts_a_l    = nphip_l + firstp_l;
    
    int nphiu_r       = datavecright[qb].fVecShapeIndex.NElements();
    int nphip_r       = phi_ps_r.Rows();
    int nphis_a_r     = phi_ss_r.Rows();
    int firstu_r      = firsts_a_l + nphis_a_l;
    int firstp_r      = nphiu_r + firstu_r;
    int firsts_a_r    = nphip_r + firstp_r;
    
    TPZManVector<STATE,3> n = data.normal;
    TPZManVector<REAL,3> u_l  = datavecleft[qb].sol[0];
    REAL p_l                  = datavecleft[pb].sol[0][0];
    REAL s_l                  = datavecleft[sb_a].sol[0][0];

    TPZManVector<REAL,3> u_r  = datavecright[qb].sol[0];
    REAL p_r                  = datavecright[pb].sol[0][0];
    REAL s_r                  = datavecright[sb_a].sol[0][0];
    
    STATE un_l = 0.0, un_r = 0.0;
    
    for (int i = 0; i < u_l.size(); i++) {
        un_l += u_l[i]*n[i];
        un_r += u_r[i]*n[i];
    }
    
    
    //  Average values p_a
    
    STATE p_a_l    = p_l;
    STATE s_a_l    = s_l;
    STATE p_a_r    = p_r;
    STATE s_a_r    = s_r;
    
    STATE beta = 0.0;
    // upwinding
    if (un_l > 0) {
        beta = 1.0;
    }
    
    TPZManVector<STATE, 10> fa_l,v_l(nvars+1),fa_r,v_r(nvars+1);
    v_l[0] = p_a_l;
    v_l[1] = s_a_l;
    v_r[0] = p_a_r;
    v_r[1] = s_a_r;
    
    this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
    this->fSimulationData->PetroPhysics()->fa(fa_r, v_r);
    
    for (int is = 0; is < nphis_a_l; is++) {
        
        ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ss_l(is,0)*un_l;
    }
    
    for (int is = 0; is < nphis_a_r; is++) {
        
        ef(is + firsts_a_r) += -1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ss_r(is,0)*un_l;
    }
    
    
}

void TRMMultiphase::Solution_ab(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int qb = 0;
    int pb = 1;
    int sb = 2;
    
    TPZManVector<REAL,3> u = datavec[qb].sol[0];
    REAL p = datavec[pb].sol[0][0];
    REAL s = datavec[sb].sol[0][0];
    
    TPZFMatrix<STATE> dudx = datavec[qb].dsol[0];
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
        case 4:
        {
            Solout[0] = 1.0-s;
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


void TRMMultiphase::Contribute_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int nvars = 4; // {p,sa,sb,t}
    
    int qb      = 0;
    int pb      = 1;
    int sb_a    = 2;
    int sb_b    = 3;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<100,STATE> phi_ssa       = datavec[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb       = datavec[sb_b].phi;
    TPZFNMatrix<300,STATE> dphi_us      = datavec[qb].dphix;
    TPZFNMatrix<100,STATE> dphi_ps      = datavec[pb].dphix;
    
    TPZFNMatrix<40,STATE> div_on_master;
    STATE divflux;
    this->ComputeDivergenceOnMaster(datavec, div_on_master,divflux);
    REAL jac_det = datavec[qb].detjac;
    
    int nphiu       = datavec[qb].fVecShapeIndex.NElements();
    int nphip       = phi_ps.Rows();
    int nphis_a     = phi_ssa.Rows();
    int nphis_b     = phi_ssb.Rows();
    int firstu      = 0;
    int firstp      = nphiu + firstu;
    int firsts_a    = nphip + firstp;
    int firsts_b    = nphis_a + firsts_a;
    
    TPZManVector<REAL,3> u  = datavec[qb].sol[0];
    REAL p                  = datavec[pb].sol[0][0];
    REAL sa                 = datavec[sb_a].sol[0][0];
    REAL sb                 = datavec[sb_b].sol[0][0];
    
    TPZFNMatrix<10,STATE> Graduaxes = datavec[qb].dsol[0];
    
    // Time
    STATE dt = fSimulationData->dt();
    
    //  Average values p_a
    
    REAL p_a    = p;
    REAL s_a    = sa;
    REAL s_b    = sb;
    
    //  Computing closure relationship at given average values
    
    TPZManVector<STATE, 10> v(nvars);
    v[0] = p_a;
    v[1] = s_a;
    v[2] = s_b;
    
    // Fluid parameters
    TPZManVector<STATE, 10> rho_a,rho_b,rho_c,l;
    fSimulationData->AlphaProp()->Density(rho_a, v);
    fSimulationData->BetaProp()->Density(rho_b, v);
    fSimulationData->GammaProp()->Density(rho_c, v);
    fSimulationData->PetroPhysics()->l(l, v);
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    TPZManVector<STATE, 10> phi;
    fSimulationData->Map()->Kappa(datavec[qb].x, K, Kinv, v);
    fSimulationData->Map()->phi(datavec[qb].x, phi, v);
    
    // Defining local variables
    TPZFNMatrix<3,STATE> lambda_K_inv_u(3,1),lambda_dp_K_inv_u(3,1), lambda_dsa_K_inv_u(3,1), lambda_dsb_K_inv_u(3,1), lambda_K_inv_phi_u_j(3,1);
    TPZManVector<STATE,3> Gravity = fSimulationData->Gravity();
    
    for (int i = 0; i < u.size(); i++) {
        STATE dot = 0.0;
        for (int j =0; j < u.size(); j++) {
            dot += Kinv(i,j)*u[j];
        }
        lambda_K_inv_u(i,0)     = (1.0/l[0]) * dot;
        lambda_dp_K_inv_u(i,0)  = (-l[1]/(l[0]*l[0])) * dot;
        lambda_dsa_K_inv_u(i,0)  = (-l[2]/(l[0]*l[0])) * dot;
        lambda_dsb_K_inv_u(i,0)  = (-l[3]/(l[0]*l[0])) * dot;
    }
    
    // Integration point contribution
    STATE divu = 0.0;
    TPZFNMatrix<3,STATE> phi_u_i(3,1), phi_u_j(3,1);
    
    int s_i, s_j;
    int v_i, v_j;
    
    if(! fSimulationData->IsCurrentStateQ()){
        for (int ip = 0; ip < nphip; ip++)
        {
            
            ef(ip + firstp) += -1.0 * weight * (-1.0/dt) * (sa*rho_a[0]+sb*rho_b[0]+(1.0-sa-sb)*rho_c[0]) * phi[0] * phi_ps(ip,0);
            
        }
        
        for (int is = 0; is < nphis_a; is++)
        {
            
            ef(is + firsts_a) += weight * (-1.0/dt) * sa * rho_a[0] * phi[0] * phi_ssa(is,0);
            
        }
        
        for (int is = 0; is < nphis_b; is++)
        {
            
            ef(is + firsts_b) += weight * (-1.0/dt) * sb * rho_b[0] * phi[0] * phi_ssb(is,0);
            
        }
        
        return;
    }
    
    for (int iu = 0; iu < nphiu; iu++)
    {
        
        v_i = datavec[qb].fVecShapeIndex[iu].first;
        s_i = datavec[qb].fVecShapeIndex[iu].second;
        
        STATE Kl_inv_dot_u = 0.0, Kl_dp_inv_dot_u = 0.0, Kl_dsa_inv_dot_u = 0.0, Kl_dsb_inv_dot_u = 0.0, rho_g_dot_phi_u = 0.0, rho_dp_g_dot_phi_u = 0.0, rho_dsa_g_dot_phi_u = 0.0, rho_dsb_g_dot_phi_u = 0.0;
        for (int i = 0; i < u.size(); i++) {
            phi_u_i(i,0) = phi_us(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
            Kl_inv_dot_u        += lambda_K_inv_u(i,0)*phi_u_i(i,0);
            Kl_dp_inv_dot_u     += lambda_dp_K_inv_u(i,0)*phi_u_i(i,0);
            Kl_dsa_inv_dot_u    += lambda_dsa_K_inv_u(i,0)*phi_u_i(i,0);
            Kl_dsb_inv_dot_u    += lambda_dsb_K_inv_u(i,0)*phi_u_i(i,0);
            rho_g_dot_phi_u     += (sa*rho_a[0]+sb*rho_b[0]+(1.0-sa-sb)*rho_c[0])*Gravity[i]*phi_u_i(i,0);
            rho_dp_g_dot_phi_u  += (sa*rho_a[1]+sb*rho_b[1]+(1.0-sa-sb)*rho_c[1])*Gravity[i]*phi_u_i(i,0);
            rho_dsa_g_dot_phi_u += (rho_a[0]-rho_c[0])*Gravity[i]*phi_u_i(i,0);
            rho_dsb_g_dot_phi_u += (rho_b[0]-rho_c[0])*Gravity[i]*phi_u_i(i,0);
        }
        
        ef(iu + firstu) += weight * ( Kl_inv_dot_u - (1.0/jac_det) * (p) * div_on_master(iu,0) - rho_g_dot_phi_u);
        
        for (int ju = 0; ju < nphiu; ju++)
        {
            
            v_j = datavec[qb].fVecShapeIndex[ju].first;
            s_j = datavec[qb].fVecShapeIndex[ju].second;
            
            STATE Kl_inv_phi_u_j_dot_phi_u_i = 0.0;
            for (int j = 0; j < u.size(); j++) {
                phi_u_j(j,0) = phi_us(s_j,0) * datavec[qb].fDeformedDirections(j,v_j);
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
        
        for (int js = 0; js < nphis_a; js++)
        {
            ek(iu + firstu, js + firsts_a) += weight * ( Kl_dsa_inv_dot_u + rho_dsa_g_dot_phi_u) * phi_ssa(js,0);
        }
        
        for (int js = 0; js < nphis_b; js++)
        {
            ek(iu + firstu, js + firsts_b) += weight * ( Kl_dsb_inv_dot_u + rho_dsb_g_dot_phi_u) * phi_ssb(js,0);
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
        
        ef(ip + firstp) += -1.0 * weight * (divu + (1.0/dt) * (sa*rho_a[0]+sb*rho_b[0]+(1.0-sa-sb)*rho_c[0]) * phi[0] - f[0]) * phi_ps(ip,0);
        
        for (int ju = 0; ju < nphiu; ju++)
        {
            ek(ip + firstp, ju + firstu) += -1.0 * weight * (1.0/jac_det) * div_on_master(ju,0) * phi_ps(ip,0);
        }
        
        for (int jp = 0; jp < nphip; jp++)
        {
            ek(ip + firstp, jp + firstp) += -1.0 * weight * ( (1.0/dt) * ((sa*rho_a[0]+sb*rho_b[0]+(1.0-sa-sb)*rho_c[0]) * phi[1] + (sa*rho_a[1]+sb*rho_b[1]+(1.0-sa-sb)*rho_c[1]) * phi[0]) * phi_ps(ip,0) ) * phi_ps(jp,0);
        }
        
        for (int js = 0; js < nphis_a; js++)
        {
            ek(ip + firstp, js + firsts_a) += -1.0 * weight * ( (1.0/dt) * ((rho_a[0]-rho_c[0]) * phi[0]) * phi_ps(ip,0) ) * phi_ssa(js,0);
        }
        
        for (int js = 0; js < nphis_b; js++)
        {
            ek(ip + firstp, js + firsts_b) += -1.0 * weight * ( (1.0/dt) * ((rho_b[0]-rho_c[0]) * phi[0]) * phi_ps(ip,0) ) * phi_ssb(js,0);
        }
        
    }
    
    
    for (int is = 0; is < nphis_a; is++)
    {
        
        ef(is + firsts_a) += weight * (1.0/dt) * sa * rho_a[0] * phi[0] * phi_ssa(is,0);
        
        for (int jp = 0; jp < nphip; jp++)
        {
            ek(is + firsts_a, jp + firstp) += weight * ( (1.0/dt) * ( sa*rho_a[0] * phi[1] + sa*rho_a[1] * phi[0]) * phi_ssa(is,0) ) * phi_ps(jp,0);
        }
        
        for (int js = 0; js < nphis_a; js++)
        {
            ek(is + firsts_a, js + firsts_a) += weight * (1.0/dt) * rho_a[0] * phi[0] * phi_ssa(js,0) * phi_ssa(is,0);
        }
        
    }
    
    for (int is = 0; is < nphis_b; is++)
    {
        
        ef(is + firsts_b) += weight * (1.0/dt) * sb * rho_b[0] * phi[0] * phi_ssb(is,0);
        
        for (int jp = 0; jp < nphip; jp++)
        {
            ek(is + firsts_b, jp + firstp) += weight * ( (1.0/dt) * ( sb*rho_b[0] * phi[1] + sb*rho_b[1] * phi[0]) * phi_ssb(is,0) ) * phi_ps(jp,0);
        }
        
        for (int js = 0; js < nphis_b; js++)
        {
            ek(is + firsts_b, js + firsts_b) += weight * (1.0/dt) * rho_b[0] * phi[0] * phi_ssb(js,0) * phi_ssb(is,0);
        }
        
    }
    
}

void TRMMultiphase::Contribute_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    // @omar:: fake ek
    TPZFMatrix<STATE> ek;
    this->Contribute_abc(datavec, weight, ek, ef);
    
}

void TRMMultiphase::ContributeBC_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int qb = 0;
    int pb = 1;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[qb].phi;
    
    int nphiu       = phi_us.Rows();
    int firstu      = 0;
    
    TPZManVector<REAL,3> u  = datavec[qb].sol[0];
    
    REAL Value_m    = 0.0;
    REAL Value_sa    = 0.0;
    REAL Value_sb    = 0.0;
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<REAL,2> f(3);
        TPZFMatrix<REAL> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavec[pb].x, time, f, gradf);
        Value_m = f[0];
        Value_sa = f[1];
        Value_sb = f[2];
    }
    else{
        Value_m  = bc.Val2()(0,0);
        Value_sa = bc.Val2()(1,0);
        Value_sb = bc.Val2()(2,0);
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

void TRMMultiphase::ContributeBCInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    
    int qb      = 0;
    int pb      = 1;
    int sb_a    = 2;
    int sb_b    = 3;
    
    TPZFNMatrix<100,STATE> phi_us_l       = datavecleft[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps_l       = datavecleft[pb].phi;
    TPZFNMatrix<100,STATE> phi_ssa_l      = datavecleft[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb_l      = datavecleft[sb_b].phi;
    TPZFNMatrix<300,STATE> dphi_us_l      = datavecleft[qb].dphix;
    TPZFNMatrix<100,STATE> dphi_ps_l      = datavecleft[pb].dphix;
    
    int nphiu_l       = datavecleft[qb].fVecShapeIndex.NElements();
    int nphip_l       = phi_ps_l.Rows();
    int nphis_a_l     = phi_ssa_l.Rows();
    int nphis_b_l     = phi_ssb_l.Rows();
    int firstu_l      = 0;
    int firstp_l      = nphiu_l + firstu_l;
    int firsts_a_l    = nphip_l + firstp_l;
    int firsts_b_l    = nphis_a_l + firsts_a_l;
    
    TPZManVector<STATE,3> n = data.normal;
    TPZManVector<REAL,3> u_l  = datavecleft[qb].sol[0];
    REAL p_l                  = datavecleft[pb].sol[0][0];
    REAL sa_l                  = datavecleft[sb_a].sol[0][0];
    REAL sb_l                  = datavecleft[sb_b].sol[0][0];
    
    STATE un_l = 0.0;
    
    for (int i = 0; i < u_l.size(); i++) {
        un_l += u_l[i]*n[i];
    }
    
    //  Average values p_a
    
    REAL p_a_l    = p_l;
    REAL s_a_l    = sa_l;
    REAL s_b_l    = sb_l;
    
    REAL beta = 0.0;
    
    
    TPZManVector<STATE, 10> fa_l,fb_l,v_l(nvars+1);
    
    TPZFNMatrix<3,STATE> phi_u_i_l(3,1);
//    STATE phi_un_l = 0.0;
//    int s_j;
//    int v_j;
    
    REAL Value_m    = 0.0;
    REAL Value_sa    = 0.0;
    REAL Value_sb    = 0.0;
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<REAL,2> f(3);
        TPZFMatrix<REAL> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavecleft[qb].x, time, f, gradf);
        Value_m = f[0];
        Value_sa = f[1];
        Value_sb = f[2];
    }
    else{
        Value_m = bc.Val2()(0,0);
        Value_sa = bc.Val2()(1,0);
        Value_sb = bc.Val2()(2,0);
    }
    
    switch (bc.Type()) {
            
        case 0 :    // Dirichlet BC  PD outlet
        {
            
            // upwinding
            if (un_l > 0) {
                beta = 1.0;
            }
            
            STATE p_D = Value_m;
            
            v_l[0] = p_D;
            v_l[1] = s_a_l;
            v_l[2] = s_b_l;
            
            this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
            this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
            
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0])*phi_ssa_l(is,0)*un_l;
                
//                for (int ju = 0; ju < nphiu_l; ju++) {
//                    
//                    v_j = datavecleft[qb].fVecShapeIndex[ju].first;
//                    s_j = datavecleft[qb].fVecShapeIndex[ju].second;
//                    
//                    for (int j = 0; j < u_l.size(); j++) {
//                        phi_u_i_l(j,0) = phi_us_l(s_j,0) * datavecleft[qb].fDeformedDirections(j,v_j);
//                        phi_un_l += phi_u_i_l(j,0)*n[j];
//                    }
//                    
//                    ek(is + firsts_a_l, ju + firstu_l) += +1.0*weight * beta * fa_l[0] * phi_ssa_l(is,0)*phi_un_l;
//                }
                
                for (int js = 0; js < nphis_a_l; js++) {
                    ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ssa_l(js,0) * phi_ssa_l(is,0)*un_l;
                }
                
                
            }
            
            for (int is = 0; is < nphis_b_l; is++) {
                
                ef(is + firsts_b_l) += +1.0*weight * (beta*fb_l[0])*phi_ssb_l(is,0)*un_l;
                
//                for (int ju = 0; ju < nphiu_l; ju++) {
//                    
//                    v_j = datavecleft[qb].fVecShapeIndex[ju].first;
//                    s_j = datavecleft[qb].fVecShapeIndex[ju].second;
//                    
//                    for (int j = 0; j < u_l.size(); j++) {
//                        phi_u_i_l(j,0) = phi_us_l(s_j,0) * datavecleft[qb].fDeformedDirections(j,v_j);
//                        phi_un_l += phi_u_i_l(j,0)*n[j];
//                    }
//                    
//                    ek(is + firsts_b_l, ju + firstu_l) += +1.0*weight * beta*fb_l[0] * phi_ssb_l(is,0)*phi_un_l;
//                }
                
                for (int js = 0; js < nphis_b_l; js++) {
                    ek(is + firsts_b_l, js + firsts_b_l) += +1.0*weight * beta * fb_l[3] * phi_ssb_l(js,0) * phi_ssb_l(is,0)*un_l;
                }
                
                
            }
            
            
        }
            break;
            
        case 1 :    // Neumann BC  QN outlet
        {
            
            // upwinding
            if (Value_m > 0) {
                beta = 1.0;
            }
            
            STATE un_N = Value_m;
            
            v_l[0] = p_a_l;
            v_l[1] = s_a_l;
            v_l[2] = s_b_l;
            
            this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
            this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ssa_l(is,0)*un_N;
                
                
//                for (int jp = 0; jp < nphip_l; jp++) {
//                    ek(is + firsts_a_l, jp + firstp_l) += +1.0*weight * beta * fa_l[1] * phi_ps_l(jp,0) * phi_ssa_l(is,0)*un_N;
//                }
                
                
                for (int js = 0; js < nphis_a_l; js++) {
                    ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ssa_l(js,0) * phi_ssa_l(is,0)*un_N;
                }
                
                
            }
            
            for (int is = 0; is < nphis_b_l; is++) {
                
                ef(is + firsts_b_l) += +1.0*weight * beta*fb_l[0]*phi_ssb_l(is,0)*un_N;
                
                
//                for (int jp = 0; jp < nphip_l; jp++) {
//                    ek(is + firsts_b_l, jp + firstp_l) += +1.0*weight * beta * fb_l[1] * phi_ps_l(jp,0) * phi_ssb_l(is,0)*un_N;
//                }
                
                
                for (int js = 0; js < nphis_b_l; js++) {
                    ek(is + firsts_b_l, js + firsts_b_l) += +1.0*weight * beta * fb_l[3] * phi_ssb_l(js,0) * phi_ssb_l(is,0)*un_N;
                }
                
                
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD inlet
        {
            
            // upwinding
            if (un_l > 0) {
                beta = 1.0;
            }
            
            STATE p_D = Value_m;
            
            v_l[0] = p_D;
            v_l[1] = Value_sa;
            v_l[2] = Value_sb;
            
            this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
            this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ssa_l(is,0)*un_l;
                
//                for (int ju = 0; ju < nphiu_l; ju++) {
//                    
//                    v_j = datavecleft[qb].fVecShapeIndex[ju].first;
//                    s_j = datavecleft[qb].fVecShapeIndex[ju].second;
//                    
//                    for (int j = 0; j < u_l.size(); j++) {
//                        phi_u_i_l(j,0) = phi_us_l(s_j,0) * datavecleft[qb].fDeformedDirections(j,v_j);
//                        phi_un_l += phi_u_i_l(j,0)*n[j];
//                    }
//                    
//                    ek(is + firsts_a_l, ju + firstu_l) += +1.0*weight * beta*fa_l[0] * phi_ssa_l(is,0)*phi_un_l;
//                }
                
                
            }
            
            for (int is = 0; is < nphis_b_l; is++) {
                
                ef(is + firsts_b_l) += +1.0*weight * beta*fb_l[0]*phi_ssb_l(is,0)*un_l;
                
//                for (int ju = 0; ju < nphiu_l; ju++) {
//                    
//                    v_j = datavecleft[qb].fVecShapeIndex[ju].first;
//                    s_j = datavecleft[qb].fVecShapeIndex[ju].second;
//                    
//                    for (int j = 0; j < u_l.size(); j++) {
//                        phi_u_i_l(j,0) = phi_us_l(s_j,0) * datavecleft[qb].fDeformedDirections(j,v_j);
//                        phi_un_l += phi_u_i_l(j,0)*n[j];
//                    }
//                    
//                    ek(is + firsts_b_l, ju + firstu_l) += +1.0*weight * beta*fb_l[0] * phi_ssb_l(is,0)*phi_un_l;
//                }
                
                
            }
            
        }
            break;
            
        case 3 :    // Neumann BC  QN inlet
        {
            
            // upwinding
            if (Value_m < 0) {
                beta = 1.0;
            }
            
            STATE un_N = Value_m;
            
            v_l[0] = p_a_l;
            v_l[1] = Value_sa;
            v_l[2] = Value_sb;
            
            this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
            this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ssa_l(is,0)*un_N;
                
//                for (int jp = 0; jp < nphip_l; jp++) {
//                    ek(is + firsts_a_l, jp + firstp_l) += +1.0*weight * beta * fa_l[1] * phi_ps_l(jp,0) * phi_ssa_l(is,0)*un_N;
//                }
                
            }
            
            for (int is = 0; is < nphis_b_l; is++) {
                
                ef(is + firsts_b_l) += +1.0*weight * beta*fb_l[0]*phi_ssb_l(is,0)*un_N;
                
//                for (int jp = 0; jp < nphip_l; jp++) {
//                    ek(is + firsts_b_l, jp + firstp_l) += +1.0*weight * beta * fb_l[1] * phi_ps_l(jp,0) * phi_ssb_l(is,0)*un_N;
//                }
                
            }
            
        }
            break;
            
        case 4 :    // Neumann BC  Impervious bc
        {
            
            // upwinding
            beta = 1.0;
            
            STATE un_N = 0.0;
            
            v_l[0] = p_a_l;
            v_l[1] = Value_sa;
            v_l[2] = Value_sb;
            
            this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
            this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ssa_l(is,0)*un_N;
                
//                for (int jp = 0; jp < nphip_l; jp++) {
//                    ek(is + firsts_a_l, jp + firstp_l) += +1.0*weight * beta * fa_l[1] * phi_ps_l(jp,0) * phi_ssa_l(is,0)*un_N;
//                }
                
            }
            
            for (int is = 0; is < nphis_b_l; is++) {
                
                ef(is + firsts_b_l) += +1.0*weight * beta*fb_l[0]*phi_ssb_l(is,0)*un_N;
                
//                for (int jp = 0; jp < nphip_l; jp++) {
//                    ek(is + firsts_b_l, jp + firstp_l) += +1.0*weight * beta * fb_l[1] * phi_ps_l(jp,0) * phi_ssb_l(is,0)*un_N;
//                }
                
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

void TRMMultiphase::ContributeBCInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    DebugStop();
    
}

void TRMMultiphase::ContributeInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    
    int qb      = 0;
    int pb      = 1;
    int sb_a    = 2;
    int sb_b    = 3;
    
    TPZFNMatrix<100,STATE> phi_us_l       = datavecleft[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps_l       = datavecleft[pb].phi;
    TPZFNMatrix<100,STATE> phi_ssa_l      = datavecleft[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb_l      = datavecleft[sb_b].phi;
    TPZFNMatrix<300,STATE> dphi_us_l      = datavecleft[qb].dphix;
    TPZFNMatrix<100,STATE> dphi_ps_l      = datavecleft[pb].dphix;
    
    TPZFNMatrix<100,STATE> phi_us_r       = datavecright[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps_r       = datavecright[pb].phi;
    TPZFNMatrix<100,STATE> phi_ssa_r      = datavecright[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb_r      = datavecright[sb_b].phi;
    TPZFNMatrix<300,STATE> dphi_us_r      = datavecright[qb].dphix;
    TPZFNMatrix<100,STATE> dphi_ps_r      = datavecright[pb].dphix;
    
    int nphiu_l       = datavecleft[qb].fVecShapeIndex.NElements();
    int nphip_l       = phi_ps_l.Rows();
    int nphis_a_l     = phi_ssa_l.Rows();
    int nphis_b_l     = phi_ssb_l.Rows();
    int firstu_l      = 0;
    int firstp_l      = nphiu_l + firstu_l;
    int firsts_a_l    = nphip_l + firstp_l;
    int firsts_b_l    = nphis_a_l + firsts_a_l;
    
    int nphiu_r       = datavecright[qb].fVecShapeIndex.NElements();
    int nphip_r       = phi_ps_r.Rows();
    int nphis_a_r     = phi_ssa_r.Rows();
    int nphis_b_r     = phi_ssb_r.Rows();
    int firstu_r      = firsts_b_l + nphis_b_l;
    int firstp_r      = nphiu_r + firstu_r;
    int firsts_a_r    = nphip_r + firstp_r;
    int firsts_b_r    = nphis_a_r + firsts_a_r;
    
    TPZManVector<STATE,3> n = data.normal;
    TPZManVector<REAL,3> u_l  = datavecleft[qb].sol[0];
    REAL p_l                  = datavecleft[pb].sol[0][0];
    REAL sa_l                 = datavecleft[sb_a].sol[0][0];
    REAL sb_l                 = datavecleft[sb_b].sol[0][0];
    
    TPZManVector<REAL,3> u_r  = datavecright[qb].sol[0];
    REAL p_r                  = datavecright[pb].sol[0][0];
    REAL sa_r                 = datavecright[sb_a].sol[0][0];
    REAL sb_r                 = datavecright[sb_b].sol[0][0];
    
    STATE un_l = 0.0, un_r = 0.0;
    
    for (int i = 0; i < u_l.size(); i++) {
        un_l += u_l[i]*n[i];
        un_r += u_r[i]*n[i];
    }
    
    
    //  Average values p_a
    
    STATE p_a_l    = p_l;
    STATE s_a_l    = sa_l;
    STATE s_b_l    = sb_l;
    STATE p_a_r    = p_r;
    STATE s_a_r    = sa_r;
    STATE s_b_r    = sb_r;
    
    STATE beta = 0.0;
    // upwinding
    if (un_l > 0.0) {
        beta = 1.0;
    }
    
    TPZManVector<STATE, 10> fa_l,fb_l,v_l(nvars+1),fa_r,fb_r,v_r(nvars+1);
    v_l[0] = p_a_l;
    v_l[1] = s_a_l;
    v_l[2] = s_b_l;
    v_r[0] = p_a_r;
    v_r[1] = s_a_r;
    v_r[2] = s_b_r;
    
    this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
    this->fSimulationData->PetroPhysics()->fa_3p(fa_r, v_r);
    
    this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
    this->fSimulationData->PetroPhysics()->fb_3p(fb_r, v_r);
    
    TPZFNMatrix<3,STATE> phi_u_i_l(3,1);
//    STATE phi_un_l = 0.0;
//    int s_j;
//    int v_j;
    
    for (int is = 0; is < nphis_a_l; is++) {
        
        ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ssa_l(is,0)*un_l;
        
//        for (int ju = 0; ju < nphiu_l; ju++) {
//            
//            v_j = datavecleft[qb].fVecShapeIndex[ju].first;
//            s_j = datavecleft[qb].fVecShapeIndex[ju].second;
//            
//            for (int j = 0; j < u_l.size(); j++) {
//                phi_u_i_l(j,0) = phi_us_l(s_j,0) * datavecleft[qb].fDeformedDirections(j,v_j);
//                phi_un_l += phi_u_i_l(j,0)*n[j];
//            }
//            
//            ek(is + firsts_a_l, ju + firstu_l) += +1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0]) * phi_ssa_l(is,0)*phi_un_l;
//        }
//        
//        for (int jp = 0; jp < nphip_l; jp++) {
//            ek(is + firsts_a_l, jp + firstp_l) += +1.0*weight * beta * fa_l[1] * phi_ps_l(jp,0) * phi_ssa_l(is,0)*un_l;
//        }
//        
//        for (int jp = 0; jp < nphip_r; jp++) {
//            ek(is + firsts_a_l, jp + firstp_r) += +1.0*weight * (1.0-beta) * fa_r[1] * phi_ps_r(jp,0) * phi_ssa_l(is,0)*un_l;
//        }
        
        for (int js = 0; js < nphis_a_l; js++) {
            ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ssa_l(js,0) * phi_ssa_l(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_a_l, js + firsts_a_r) += +1.0*weight * (1.0-beta) * fa_r[2] * phi_ssa_r(js,0) * phi_ssa_l(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_b_l; js++) {
            ek(is + firsts_a_l, js + firsts_b_l) += +1.0*weight * beta * fa_l[3] * phi_ssb_l(js,0) * phi_ssa_l(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_a_l, js + firsts_b_r) += +1.0*weight * (1.0-beta) * fa_r[3] * phi_ssb_r(js,0) * phi_ssa_l(is,0)*un_l;
        }
        
    }
    
    for (int is = 0; is < nphis_a_r; is++) {
        
        ef(is + firsts_a_r) += -1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ssa_r(is,0)*un_l;
        
//        for (int ju = 0; ju < nphiu_l; ju++) {
//            
//            v_j = datavecleft[qb].fVecShapeIndex[ju].first;
//            s_j = datavecleft[qb].fVecShapeIndex[ju].second;
//            
//            for (int j = 0; j < u_l.size(); j++) {
//                phi_u_i_l(j,0) = phi_us_l(s_j,0) * datavecleft[qb].fDeformedDirections(j,v_j);
//                phi_un_l += phi_u_i_l(j,0)*n[j];
//            }
//            
//            ek(is + firsts_a_r, ju + firstu_l) += -1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ssa_r(is,0)*phi_un_l;
//        }
//        
//        for (int jp = 0; jp < nphip_l; jp++) {
//            ek(is + firsts_a_r, jp + firstp_l) += -1.0*weight * beta * fa_l[1] * phi_ps_l(jp,0) * phi_ssa_r(is,0)*un_l;
//        }
//        
//        for (int jp = 0; jp < nphip_r; jp++) {
//            ek(is + firsts_a_r, jp + firstp_r) += -1.0*weight * (1.0-beta) * fa_r[1] * phi_ps_r(jp,0) * phi_ssa_r(is,0)*un_l;
//        }
        
        for (int js = 0; js < nphis_a_l; js++) {
            ek(is + firsts_a_r, js + firsts_a_l) += -1.0*weight * beta * fa_l[2] * phi_ssa_l(js,0) * phi_ssa_r(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_a_r, js + firsts_a_r) += -1.0*weight * (1.0-beta) * fa_r[2] * phi_ssa_r(js,0) * phi_ssa_r(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_b_l; js++) {
            ek(is + firsts_a_r, js + firsts_b_l) += -1.0*weight * beta * fa_l[3] * phi_ssb_l(js,0) * phi_ssa_r(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_b_r; js++) {
            ek(is + firsts_a_r, js + firsts_b_r) += -1.0*weight * (1.0-beta) * fa_r[3] * phi_ssb_r(js,0) * phi_ssa_r(is,0)*un_l;
        }
        
    }
    
    for (int is = 0; is < nphis_b_l; is++) {
        
        ef(is + firsts_b_l) += +1.0*weight * (beta*fb_l[0] + (1.0-beta)*fb_r[0])*phi_ssb_l(is,0)*un_l;
        
//        for (int ju = 0; ju < nphiu_l; ju++) {
//            
//            v_j = datavecleft[qb].fVecShapeIndex[ju].first;
//            s_j = datavecleft[qb].fVecShapeIndex[ju].second;
//            
//            for (int j = 0; j < u_l.size(); j++) {
//                phi_u_i_l(j,0) = phi_us_l(s_j,0) * datavecleft[qb].fDeformedDirections(j,v_j);
//                phi_un_l += phi_u_i_l(j,0)*n[j];
//            }
//            
//            ek(is + firsts_b_l, ju + firstu_l) += +1.0*weight * (beta*fb_l[0] + (1.0-beta)*fb_r[0]) * phi_ssb_l(is,0)*phi_un_l;
//        }
//        
//        for (int jp = 0; jp < nphip_l; jp++) {
//            ek(is + firsts_b_l, jp + firstp_l) += +1.0*weight * beta * fb_l[1] * phi_ps_l(jp,0) * phi_ssb_l(is,0)*un_l;
//        }
//        
//        for (int jp = 0; jp < nphip_r; jp++) {
//            ek(is + firsts_b_l, jp + firstp_r) += +1.0*weight * (1.0-beta) * fb_r[1] * phi_ps_r(jp,0) * phi_ssb_l(is,0)*un_l;
//        }
        
        for (int js = 0; js < nphis_a_l; js++) {
            ek(is + firsts_b_l, js + firsts_a_l) += +1.0*weight * beta * fb_l[2] * phi_ssa_l(js,0) * phi_ssb_l(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_b_l, js + firsts_a_r) += +1.0*weight * (1.0-beta) * fb_r[2] * phi_ssa_r(js,0) * phi_ssb_l(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_b_l; js++) {
            ek(is + firsts_b_l, js + firsts_b_l) += +1.0*weight * beta * fb_l[3] * phi_ssb_l(js,0) * phi_ssb_l(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_b_r; js++) {
            ek(is + firsts_b_l, js + firsts_b_r) += +1.0*weight * (1.0-beta) * fb_r[3] * phi_ssb_r(js,0) * phi_ssb_l(is,0)*un_l;
        }
        
    }
    
    for (int is = 0; is < nphis_b_r; is++) {
        
        ef(is + firsts_b_r) += -1.0*weight * (beta*fb_l[0] + (1.0-beta)*fb_r[0])*phi_ssb_r(is,0)*un_l;
        
//        for (int ju = 0; ju < nphiu_l; ju++) {
//            
//            v_j = datavecleft[qb].fVecShapeIndex[ju].first;
//            s_j = datavecleft[qb].fVecShapeIndex[ju].second;
//            
//            for (int j = 0; j < u_l.size(); j++) {
//                phi_u_i_l(j,0) = phi_us_l(s_j,0) * datavecleft[qb].fDeformedDirections(j,v_j);
//                phi_un_l += phi_u_i_l(j,0)*n[j];
//            }
//            
//            ek(is + firsts_b_r, ju + firstu_l) += -1.0*weight * (beta*fb_l[0] + (1.0-beta)*fb_r[0])*phi_ssb_r(is,0)*phi_un_l;
//        }
//        
//        for (int jp = 0; jp < nphip_l; jp++) {
//            ek(is + firsts_b_r, jp + firstp_l) += -1.0*weight * beta * fb_l[1] * phi_ps_l(jp,0) * phi_ssb_r(is,0)*un_l;
//        }
//        
//        for (int jp = 0; jp < nphip_r; jp++) {
//            ek(is + firsts_b_r, jp + firstp_r) += -1.0*weight * (1.0-beta) * fb_r[1] * phi_ps_r(jp,0) * phi_ssb_r(is,0)*un_l;
//        }
        
        for (int js = 0; js < nphis_a_l; js++) {
            ek(is + firsts_b_r, js + firsts_a_l) += -1.0*weight * beta * fb_l[2] * phi_ssa_l(js,0) * phi_ssb_r(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_b_r, js + firsts_a_r) += -1.0*weight * (1.0-beta) * fb_r[2] * phi_ssa_r(js,0) * phi_ssb_r(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_b_l; js++) {
            ek(is + firsts_b_r, js + firsts_b_l) += -1.0*weight * beta * fb_l[3] * phi_ssb_l(js,0) * phi_ssb_r(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_b_r; js++) {
            ek(is + firsts_b_r, js + firsts_b_r) += -1.0*weight * (1.0-beta) * fb_r[3] * phi_ssb_r(js,0) * phi_ssb_r(is,0)*un_l;
        }
        
    }
    
}

void TRMMultiphase::ContributeInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    
    int qb      = 0;
    int pb      = 1;
    int sb_a    = 2;
    int sb_b    = 3;
    
    TPZFNMatrix<100,STATE> phi_us_l       = datavecleft[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps_l       = datavecleft[pb].phi;
    TPZFNMatrix<100,STATE> phi_ssa_l      = datavecleft[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb_l      = datavecleft[sb_b].phi;
    TPZFNMatrix<300,STATE> dphi_us_l      = datavecleft[qb].dphix;
    TPZFNMatrix<100,STATE> dphi_ps_l      = datavecleft[pb].dphix;
    
    TPZFNMatrix<100,STATE> phi_us_r       = datavecright[qb].phi;
    TPZFNMatrix<100,STATE> phi_ps_r       = datavecright[pb].phi;
    TPZFNMatrix<100,STATE> phi_ssa_r      = datavecright[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb_r      = datavecright[sb_b].phi;
    TPZFNMatrix<300,STATE> dphi_us_r      = datavecright[qb].dphix;
    TPZFNMatrix<100,STATE> dphi_ps_r      = datavecright[pb].dphix;
    
    int nphiu_l       = datavecleft[qb].fVecShapeIndex.NElements();
    int nphip_l       = phi_ps_l.Rows();
    int nphis_a_l     = phi_ssa_l.Rows();
    int nphis_b_l     = phi_ssb_l.Rows();
    int firstu_l      = 0;
    int firstp_l      = nphiu_l + firstu_l;
    int firsts_a_l    = nphip_l + firstp_l;
    int firsts_b_l    = nphis_a_l + firsts_a_l;
    
    int nphiu_r       = datavecright[qb].fVecShapeIndex.NElements();
    int nphip_r       = phi_ps_r.Rows();
    int nphis_a_r     = phi_ssa_r.Rows();
    int nphis_b_r     = phi_ssb_r.Rows();
    int firstu_r      = firsts_b_l + nphis_b_l;
    int firstp_r      = nphiu_r + firstu_r;
    int firsts_a_r    = nphip_r + firstp_r;
    int firsts_b_r    = nphis_a_r + firsts_a_r;
    
    TPZManVector<STATE,3> n = data.normal;
    TPZManVector<REAL,3> u_l  = datavecleft[qb].sol[0];
    REAL p_l                  = datavecleft[pb].sol[0][0];
    REAL sa_l                 = datavecleft[sb_a].sol[0][0];
    REAL sb_l                 = datavecleft[sb_b].sol[0][0];
    
    TPZManVector<REAL,3> u_r  = datavecright[qb].sol[0];
    REAL p_r                  = datavecright[pb].sol[0][0];
    REAL sa_r                 = datavecright[sb_a].sol[0][0];
    REAL sb_r                 = datavecright[sb_b].sol[0][0];
    
    STATE un_l = 0.0, un_r = 0.0;
    
    for (int i = 0; i < u_l.size(); i++) {
        un_l += u_l[i]*n[i];
        un_r += u_r[i]*n[i];
    }
    
    
    //  Average values p_a
    
    STATE p_a_l    = p_l;
    STATE s_a_l    = sa_l;
    STATE s_b_l    = sb_l;
    STATE p_a_r    = p_r;
    STATE s_a_r    = sa_r;
    STATE s_b_r    = sb_r;
    
    STATE beta = 0.0;
    // upwinding
    if (un_l > 0.0) {
        beta = 1.0;
    }
    
    TPZManVector<STATE, 10> fa_l,fb_l,v_l(nvars+1),fa_r,fb_r,v_r(nvars+1);
    v_l[0] = p_a_l;
    v_l[1] = s_a_l;
    v_l[2] = s_b_l;
    v_r[0] = p_a_r;
    v_r[1] = s_a_r;
    v_r[2] = s_b_r;
    
    this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
    this->fSimulationData->PetroPhysics()->fa_3p(fa_r, v_r);
    
    this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
    this->fSimulationData->PetroPhysics()->fb_3p(fb_r, v_r);

    
    for (int is = 0; is < nphis_a_l; is++) {
        
        ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ssa_l(is,0)*un_l;
        
    }
    
    for (int is = 0; is < nphis_a_r; is++) {
        
        ef(is + firsts_a_r) += -1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ssa_r(is,0)*un_l;
        
    }
    
    for (int is = 0; is < nphis_b_l; is++) {
        
        ef(is + firsts_b_l) += +1.0*weight * (beta*fb_l[0] + (1.0-beta)*fb_r[0])*phi_ssb_l(is,0)*un_l;
        
    }
    
    for (int is = 0; is < nphis_b_r; is++) {
        
        ef(is + firsts_b_r) += -1.0*weight * (beta*fb_l[0] + (1.0-beta)*fb_r[0])*phi_ssb_r(is,0)*un_l;
        
    }
    
}

void TRMMultiphase::Solution_abc(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
    
    int qb = 0;
    int pb = 1;
    int sba = 2;
    int sbb = 3;
    
    TPZManVector<REAL,3> u = datavec[qb].sol[0];
    REAL p = datavec[pb].sol[0][0];
    REAL sa = datavec[sba].sol[0][0];
    REAL sb = datavec[sbb].sol[0][0];
    
    TPZFMatrix<STATE> dudx = datavec[qb].dsol[0];
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
            Solout[0] = sa;
        }
            break;
        case 4:
        {
            Solout[0] = sb;
        }
            break;
        case 5:
        {
            Solout[0] = 1.0 - sa - sb;
        }
            break;
        default:
        {
            DebugStop();
        }
    }
    
}


int TRMMultiphase::ClassId() const{
    return -6378;
}

// -------------------------------------------------------------------------------------------

void TRMMultiphase::Write(TPZStream &buf, int withclassid) const{
    
    TPZMaterial::Write(buf, withclassid);
    
}

// -------------------------------------------------------------------------------------------

void TRMMultiphase::Read(TPZStream &buf, void *context) {
    TPZMaterial::Read(buf, context);
    
}

// Update element memory by copying the n+1 data to the n data
void TRMMultiphase::UpdateMemory()
{
    DebugStop();
}

void TRMMultiphase::apply_ux(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    int ub = 0;
    int qb = 1;
    int pb = 2;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,STATE> phi_qs       = datavec[qb].phi;
    
    int n_u         = fdimension;
    int nphiu       = phi_us.Rows();
    int firstu      = 0;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    TPZManVector<REAL,3> q  = datavec[qb].sol[0];
    
    REAL Value = 0.0;
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<REAL,2> f(4);
        TPZFMatrix<REAL> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavec[pb].x, time, f, gradf);
        Value = f[0];
    }
    else{
        Value = bc.Val2()(0,0);
        std::cout<< "Define the boundary as x-t function " << std::endl;
        DebugStop();
    }
    
    REAL u_x = Value;
    for (int iu = 0; iu < nphiu; iu++)
    {
        ef(n_u*iu + 0 + firstu) += weight * gBigNumber * (u[0] - u_x) * phi_us(iu,0);
        
        for (int ju = 0; ju < nphiu; ju++)
        {
            
            ek(n_u*iu + 0 + firstu, n_u*ju + 0 + firstu) += weight * gBigNumber * phi_us(ju,0) * phi_us(iu,0);
        }
    }
    
}

void TRMMultiphase::apply_uy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    int ub = 0;
    int qb = 1;
    int pb = 2;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,STATE> phi_qs       = datavec[qb].phi;
    
    int n_u         = fdimension;
    int nphiu       = phi_us.Rows();
    int firstu      = 0;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    TPZManVector<REAL,3> q  = datavec[qb].sol[0];
    
    REAL Value = 0.0;
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<REAL,2> f(4);
        TPZFMatrix<REAL> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavec[pb].x, time, f, gradf);
        Value = f[1];
    }
    else{
        Value = bc.Val2()(0,0);
        std::cout<< "Define the boundary as x-t function " << std::endl;
        DebugStop();
    }
    
    REAL u_y = Value;
    for (int iu = 0; iu < nphiu; iu++)
    {
        ef(n_u*iu + 1 + firstu) += weight * gBigNumber * (u[1] - u_y) * phi_us(iu,0);
        
        for (int ju = 0; ju < nphiu; ju++)
        {
            
            ek(n_u*iu + 1 + firstu, n_u*ju + 1 + firstu) += weight * gBigNumber * phi_us(ju,0) * phi_us(iu,0);
        }
    }
    
}

void TRMMultiphase::apply_uz(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    int ub = 0;
    int qb = 1;
    int pb = 2;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,STATE> phi_qs       = datavec[qb].phi;
    
    int n_u         = fdimension;
    int nphiu       = phi_us.Rows();
    int firstu      = 0;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    
    REAL Value = 0.0;
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<REAL,2> f(4);
        TPZFMatrix<REAL> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavec[pb].x, time, f, gradf);
        Value = f[2];
    }
    else{
        Value = bc.Val2()(0,0);
        std::cout<< "Define the boundary as x-t function " << std::endl;
        DebugStop();
    }
    
    REAL u_z = Value;
    for (int iu = 0; iu < nphiu; iu++)
    {
        ef(n_u*iu + 2 + firstu) += weight * gBigNumber * (u[3] - u_z) * phi_us(iu,0);
        
        for (int ju = 0; ju < nphiu; ju++)
        {
            
            ek(n_u*iu + 2 + firstu, n_u*ju + 2 + firstu) += weight * gBigNumber * phi_us(ju,0) * phi_us(iu,0);
        }
    }
    
}

void TRMMultiphase::apply_tn(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    int ub = 0;
    int qb = 1;
    int pb = 2;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,STATE> phi_qs       = datavec[qb].phi;
    
    int n_u         = fdimension;
    int nphiu       = phi_us.Rows();
    int firstu      = 0;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    
    REAL Tx = 0.0;
    REAL Ty = 0.0;
    REAL Tz = 0.0;
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<REAL,2> f(4);
        TPZFMatrix<REAL> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavec[pb].x, time, f, gradf);
        Tx = f[0];
        Ty = f[1];
        Tz = f[2];
    }
    else{
        Tx = bc.Val2()(0,0);
        std::cout<< "Define the boundary as x-t function " << std::endl;
        DebugStop();
    }
    
    for (int iu = 0; iu < nphiu; iu++)
    {
        ef(n_u*iu + 0 + firstu) += -1.0 * weight * Tx * phi_us(iu,0);
        ef(n_u*iu + 1 + firstu) += -1.0 * weight * Ty * phi_us(iu,0);
//        ef(n_u*iu + 2 + firstu) += -1.0 * weight * Tz * phi_us(iu,0);
    }
    
}

void TRMMultiphase::apply_p(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    int ub = 0;
    int qb = 1;
    int pb = 2;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,STATE> phi_qs       = datavec[qb].phi;
    
    int n_u         = fdimension;
    int nphiu       = phi_us.Rows();
    int nphiq       = phi_qs.Rows();
    int firstu      = 0;
    int firstq      = nphiu*n_u + firstu;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    TPZManVector<REAL,3> q  = datavec[qb].sol[0];
    
    REAL Value = 0.0;
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<REAL,2> f(4);
        TPZFMatrix<REAL> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavec[pb].x, time, f, gradf);
        Value = f[3];
    }
    else{
        Value = bc.Val2()(0,0);
        std::cout<< "Define the boundary as x-t function " << std::endl;
        DebugStop();
    }
    
    REAL p_D = Value;
    for (int iq = 0; iq < nphiq; iq++)
    {
        ef(iq + firstq) += weight * p_D * phi_qs(iq,0);
    }
    
}

void TRMMultiphase::apply_q(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    int ub = 0;
    int qb = 1;
    int pb = 2;
    
    TPZFNMatrix<100,STATE> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,STATE> phi_qs       = datavec[qb].phi;
    
    int n_u         = fdimension;
    int nphiu       = phi_us.Rows();
    int nphiq       = phi_qs.Rows();
    int firstu      = 0;
    int firstq      = nphiu*n_u + firstu;
    
    TPZManVector<REAL,3> u  = datavec[ub].sol[0];
    TPZManVector<REAL,3> q  = datavec[qb].sol[0];
    
    REAL Value = 0.0;
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<REAL,2> f(4);
        TPZFMatrix<REAL> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavec[pb].x, time, f, gradf);
        Value = f[3];
    }
    else{
        Value = bc.Val2()(0,0);
    }
    
    for (int iq = 0; iq < nphiq; iq++)
    {
        REAL qn_N = Value, qn = q[0];
        
        if (qn_N == 0.0) {

            ef(iq + firstq) += weight * 100000.0 * gBigNumber * (qn - qn_N) * phi_qs(iq,0);
            
            for (int jq = 0; jq < nphiq; jq++)
            {
                
                ek(iq + firstq,jq + firstq) += weight * 100000.0 * gBigNumber * phi_qs(jq,0) * phi_qs(iq,0);
            }
            
        }
        else{
            
            ef(iq + firstq) += weight * gBigNumber * (qn - qn_N) * phi_qs(iq,0);
            
            for (int jq = 0; jq < nphiq; jq++)
            {
                
                ek(iq + firstq,jq + firstq) += weight * gBigNumber * phi_qs(jq,0) * phi_qs(iq,0);
            }
        }
        
    }
    
}
