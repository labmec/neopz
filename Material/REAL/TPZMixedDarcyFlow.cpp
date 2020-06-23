//
//  TPZMixedDarcyFlow.cpp
//  HDiv
//
//  Created by Omar Durán on 3/29/19.
//

#include "TPZMixedDarcyFlow.h"

TPZMixedDarcyFlow::TPZMixedDarcyFlow() : TPZMaterial(){
    
    m_kappa.Resize(0, 0);
    m_kappa_inv.Resize(0, 0);
    m_gravity.resize(3);
    m_d = 0.0;
    m_dim = 0;
    
}

TPZMixedDarcyFlow::~TPZMixedDarcyFlow(){
    
}

TPZMixedDarcyFlow::TPZMixedDarcyFlow(int mat_id, int dim) :  TPZMaterial(mat_id){
    m_kappa.Resize(3, 3);
    m_kappa_inv.Resize(3, 3);
    m_gravity.resize(3);
    m_d = 0.0;
    m_dim = dim;
}


TPZMixedDarcyFlow::TPZMixedDarcyFlow(const TPZMixedDarcyFlow &other) : TPZMaterial(other){

    m_kappa         = other.m_kappa;
    m_kappa_inv     = other.m_kappa_inv;
    m_d             = other.m_d;
    m_dim           = other.m_dim;
    m_gravity       =other.m_gravity;
}

TPZMaterial * TPZMixedDarcyFlow::NewMaterial(){
    return new TPZMixedDarcyFlow(*this);
}

TPZMixedDarcyFlow & TPZMixedDarcyFlow::operator=(const TPZMixedDarcyFlow &other){
    if (this != & other) // prevent self-assignment
    {
        TPZMaterial::operator=(other);
        m_kappa         = other.m_kappa;
        m_kappa_inv     = other.m_kappa_inv;
        m_d             = other.m_d;
        m_dim           = other.m_dim;
        m_gravity       =other.m_gravity;
    }
    return *this;
}

int TPZMixedDarcyFlow::Dimension() const {
    return m_dim;
}

int TPZMixedDarcyFlow::NStateVariables() const{
    return 1;
}

void TPZMixedDarcyFlow::Print(std::ostream & out){
    m_kappa.Print(out);
    m_kappa_inv.Print(out);
    out << m_d << std::endl;
    out << m_dim << std::endl;
}

std::string TPZMixedDarcyFlow::Name(){
    return "TPZMixedDarcyFlow";
}

void TPZMixedDarcyFlow::FillDataRequirements(TPZVec<TPZMaterialData> &datavec){
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

void TPZMixedDarcyFlow::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec){
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

int TPZMixedDarcyFlow::ClassId() const{
    return Hash("TPZMixedDarcyFlow") ^ TPZMaterial::ClassId() << 1;
}

void TPZMixedDarcyFlow::Write(TPZStream &buf, int withclassid) const{
    DebugStop();
}

void TPZMixedDarcyFlow::Read(TPZStream &buf, void *context){
    DebugStop();
}

void TPZMixedDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    int qb = 0;
    int pb = 1;
    
    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].phi;
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,REAL> dphi_qs      = datavec[qb].dphix;
    TPZFNMatrix<100,REAL> dphi_ps      = datavec[pb].dphix;
    
    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = datavec[0].x[0];
    STATE val = x_spatial(0,0);
    REAL r = Norm(x_spatial);
    
    
    auto &div_phi = datavec[qb].divphi;
    REAL div_q = datavec[qb].divsol[0][0];
    
    int nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    int first_q      = 0;
    int first_p      = nphi_q + first_q;
    
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    STATE p                  = datavec[pb].sol[0][0];
    
    //axisimetria
    REAL s = 1.0;
    if (0) {
        s *= 2.0*M_PI*r;
        q[0] *= (1.0/s);
        q[1] *= (1.0/s);
   //     q[2] *= (1.0/s);
        dphi_qs *= (1.0/s);
    }
     //axisimetria
    
    TPZFNMatrix<3,STATE> phi_q_i(3,1,0.0), kappa_inv_phi_q_j(3,1,0.0), kappa_inv_q(3,1,0.0);
    
    int s_i, s_j;
    int v_i, v_j;
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            kappa_inv_q(i,0) += m_kappa_inv(i,j)*q[j];
        }
    }
    
    for (int iq = 0; iq < nphi_q; iq++)
    {
        
        v_i = datavec[qb].fVecShapeIndex[iq].first;
        s_i = datavec[qb].fVecShapeIndex[iq].second;
        
        STATE kappa_inv_q_dot_phi_q_i = 0.0;
        STATE g_dot_phi_q_i = 0.0;
        for (int i = 0; i < 3; i++) {
            phi_q_i(i,0) = phi_qs(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
            kappa_inv_q_dot_phi_q_i        += kappa_inv_q(i,0)*phi_q_i(i,0);
            g_dot_phi_q_i                  += m_gravity[i]*phi_q_i(i,0);
        }
        //errado
//        ef(iq + first_q) += -1.0 * weight * ( kappa_inv_q_dot_phi_q_i - p * div_phi(iq,0));
     
        
        ef(iq + first_q) += -1.0 * weight * (  g_dot_phi_q_i );
        
        for (int jq = 0; jq < nphi_q; jq++)
        {
            
            v_j = datavec[qb].fVecShapeIndex[jq].first;
            s_j = datavec[qb].fVecShapeIndex[jq].second;
            
            kappa_inv_phi_q_j.Zero();
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    kappa_inv_phi_q_j(i,0) += m_kappa_inv(i,j) * phi_qs(s_j,0) * datavec[qb].fDeformedDirections(j,v_j);
                }
            }
            
            STATE kappa_inv_phi_q_j_dot_phi_q_i = 0.0;
            for (int j = 0; j < 3; j++) {
                kappa_inv_phi_q_j_dot_phi_q_i += kappa_inv_phi_q_j(j,0)*phi_q_i(j,0);
            }
            
            ek(iq + first_q,jq + first_q) += weight * kappa_inv_phi_q_j_dot_phi_q_i;
        }
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(iq + first_q, jp + first_p) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);
        }
        
    }
    
    for (int ip = 0; ip < nphi_p; ip++)
    {
        STATE force = 0.0;
        if(fForcingFunction) {
            TPZManVector<STATE> res(1);
            fForcingFunction->Execute(datavec[0].x,res);
            force = res[0];
        }
        
        
        ef(ip + first_p) += -1.0 * weight * (force) * phi_ps(ip,0);
       
        for (int jq = 0; jq < nphi_q; jq++)
        {
            ek(ip + first_p, jq + first_q) += -1.0 * weight * div_phi(jq,0) * phi_ps(ip,0);
        }
        
    }
    
}

void TPZMixedDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ekfake, ef);
}

void TPZMixedDarcyFlow::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->ContributeBC(datavec, weight, ekfake, ef, bc);
}

void TPZMixedDarcyFlow::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    gBigNumber = 1e18;
    int qb = 0;
    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].phi;
    
    int nphi_q       = phi_qs.Rows();
    int first_q      = 0;
    
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    
    
    //
    
    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = datavec[0].x[0];
    REAL r = Norm(x_spatial);
    
    int npos = q.size();
    
    //axisimetria
    REAL s = 1.0;
    if (0) {
        s *= 2.0*M_PI*r;
        q[0] *= (1.0/s);
  //      q[1] *= (1.0/s);
  //      q[2] *= (1.0/s);
   
    }
    //axisimetria
    
    //
    
    
    TPZManVector<STATE,1> bc_data(1,0.0);
    bc_data[0] = bc.Val2()(0,0);
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD
        {
            STATE p_D = bc_data[0];
            for (int iq = 0; iq < nphi_q; iq++)
            {
                ef(iq + first_q) += -1.0 *  weight * p_D * phi_qs(iq,0);
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN
        {
            
            for (int iq = 0; iq < nphi_q; iq++)
            {
                REAL qn_N = bc_data[0], qn = q[0];
                ef(iq + first_q) += -1.0 * weight * gBigNumber * (qn - qn_N) * phi_qs(iq,0);
                
                for (int jq = 0; jq < nphi_q; jq++)
                {
                    
                    ek(iq + first_q,jq + first_q) += -1.0 * weight * gBigNumber * phi_qs(jq,0) * phi_qs(iq,0);
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

int TPZMixedDarcyFlow::VariableIndex(const std::string &name){
    if(!strcmp("q",name.c_str()))               return  1;
    if(!strcmp("p",name.c_str()))               return  2;
    if(!strcmp("div_q",name.c_str()))           return  3;
    if(!strcmp("kappa",name.c_str()))           return  4;
    return TPZMaterial::VariableIndex(name);
}

int TPZMixedDarcyFlow::NSolutionVariables(int var){
    if(var == 1) return 3;
    if(var == 2) return 1;
    if(var == 3) return 1;
    if(var == 4) return this->Dimension();
    
    return TPZMaterial::NSolutionVariables(var);
}

void TPZMixedDarcyFlow::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    int qb = 0;
    int pb = 1;
    Solout.Resize( this->NSolutionVariables(var));
    TPZManVector<STATE,3> p, q;
    
    q = datavec[qb].sol[0];
    p = datavec[pb].sol[0];
    REAL div_q = datavec[qb].divsol[0][0];
    
    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = datavec[0].x[0];
    REAL r = Norm(x_spatial);
    
    
    //axisimetria
    REAL s = 1.0;
    if (0) {
        s *= 2.0*M_PI*r;
        q[0] *= (1.0/s);
        q[1] *= (1.0/s);
  //      q[2] *= (1.0/s);
     
    }
    //axisimetria
    
    if(var == 1){
        for (int i=0; i < 3; i++)
        {
            Solout[i] = q[i];
        }
        return;
    }
    
    if(var == 2){
        Solout[0] = p[0];
        return;
    }
    
    if(var == 3){
        Solout[0] = div_q;
        return;
    }
    
    if(var == 4){
        for (int i  = 0; i < this->Dimension(); i++) {
            Solout[i] = m_kappa(i,i);
        }
        return;
    }
    
    DebugStop();
}
