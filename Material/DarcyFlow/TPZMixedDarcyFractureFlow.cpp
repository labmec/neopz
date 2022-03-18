//
//  TPZMixedDarcyFractureFlow.cpp
//

#include "TPZMixedDarcyFractureFlow.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"
#ifdef USING_MKL
#include "mkl.h"
#endif
#define USEBLAS

TPZMixedDarcyFractureFlow::TPZMixedDarcyFractureFlow() : TBase() {
}


TPZMixedDarcyFractureFlow::TPZMixedDarcyFractureFlow(int mat_id, int dimension) : TBase(mat_id, dimension) {
    
    
}


TPZMixedDarcyFractureFlow::TPZMixedDarcyFractureFlow(const TPZMixedDarcyFractureFlow & other) : TBase(other){
}


TPZMixedDarcyFractureFlow & TPZMixedDarcyFractureFlow::operator=(const TPZMixedDarcyFractureFlow & other){
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    return *this;
}


TPZMixedDarcyFractureFlow::~TPZMixedDarcyFractureFlow(){
    
}


void TPZMixedDarcyFractureFlow::FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        //datavec[idata].fNormalVec = true;
    }
}


void TPZMixedDarcyFractureFlow::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fDeformedDirections = true;
    }
}


//void TPZMixedDarcyFractureFlow::SetDataTransfer(TMRSDataTransfer & SimData){
//    this->mSimData = SimData;
//}


void TPZMixedDarcyFractureFlow::Print(std::ostream &out) const{
    TPZMaterial::Print(out);
}


int TPZMixedDarcyFractureFlow::VariableIndex(const std::string &name) const{
    if(!strcmp("Flux",name.c_str()))            return  1;
    if(!strcmp("Pressure",name.c_str()))        return  2;
    if(!strcmp("div_q",name.c_str()))           return  3;
    if(!strcmp("kappa",name.c_str()))           return  4;
    if(!strcmp("g_average",name.c_str()))        return  5;
    if(!strcmp("p_average",name.c_str()))        return  6;
    return TPZMaterial::VariableIndex(name);
}


int TPZMixedDarcyFractureFlow::NSolutionVariables(int var) const{
    if(var == 1) return 3;
    if(var == 2) return 1;
    if(var == 3) return 1;
    if(var == 4) return 1;
    if(var == 5) return 1;
    if(var == 6) return 1;
    return TBase::NSolutionVariables(var);
}


void TPZMixedDarcyFractureFlow::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout) {
    TBase::Solution(datavec,var,Solout);
}


void TPZMixedDarcyFractureFlow::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    // OBS: This function was taken from iMRS with very few corrections just to make it
    // compile. It might need to be optimized and refactored!
    
    int qb = 0;
    int pb = 1;
  //  int sb = 2;
    
    
    //    datavec[qb].Print(std::cout);
    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].phi;
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,REAL> dphi_qs      = datavec[qb].dphix;
    TPZFNMatrix<100,REAL> dphi_ps      = datavec[pb].dphix;
    
    
    TPZFNMatrix<40, REAL> div_phi = datavec[qb].divphi;
  //  REAL div_q = datavec[qb].divsol[0][0];
    
    int nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    int first_q      = 0;
    int first_p      = nphi_q + first_q;
    int nvecs = datavec[qb].fDeformedDirections.Cols();
    
    // first index in fVecShapeIndex corresponding to the first transverse flux
    int first_transverse_q = 0;
    // first index in fVecShapeIndex corresponding to the second transverse flux
    int second_transverse_q = 0;
    int nconnects = datavec[qb].fHDivNumConnectShape.size();
    second_transverse_q = nvecs-datavec[qb].fHDivNumConnectShape[nconnects-1];
    first_transverse_q = second_transverse_q-datavec[qb].fHDivNumConnectShape[nconnects-2];
    /*
    for(int i=0; i< nphi_q; i++)
    {
        if(first_transverse_q == 0 && datavec[qb].fVecShapeIndex[i].first == nvecs-2) first_transverse_q = i;
        if(second_transverse_q == 0 && datavec[qb].fVecShapeIndex[i].first == nvecs-1) second_transverse_q = i;
    }
     */
    if(first_transverse_q == 0 || second_transverse_q == 0 || first_transverse_q == second_transverse_q)
    {
        DebugStop();
    }
    
    
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    STATE p                  = datavec[pb].sol[0][0];
  
    
    // Get the data at integrations points
    TPZFNMatrix<3,STATE> phi_q_i(3,1,0.0), kappa_inv_phi_q_j(3,1,0.0), kappa_inv_q(3,1,0.0),kappa_inv_qFrac(3,1,0.0) ;
    
    const STATE perm = GetPermeability(datavec[0].x);
    const STATE inv_perm = 1 / perm;
    REAL kappaNormal = perm;
    REAL ad =1.0;
    REAL eps = ad;
    REAL factor = 1.0;
    
    int s_i, s_j;
    int v_i, v_j;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            kappa_inv_q(i,0) += inv_perm*q[j];
        }
    }
   
#if defined(USEBLAS) && defined(USING_MKL)
    TPZFNMatrix<3, REAL> ivec(3, first_transverse_q, 0.);
    for (int iq = 0; iq < first_transverse_q; iq++){
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        for (int id = 0; id < 3; id++) {
            ivec(id, iq) = datavec[0].fDeformedDirections(id, ivecind);
        }
    }
    {
        double *A, *B, *C;
        double alpha, beta;
        int m,n,k;
        m = first_transverse_q;
        n = first_transverse_q;
        k = 3;
        alpha = weight*(1.0/eps)*inv_perm;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDA = 3;
        LDB = 3;
        LDC = ek.Rows();
        A = &ivec(0,0);
        B = &ivec(0,0);
        C = &ek(0,0);
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
    //contribucion matrix B
    {
        double *A, *B, *C;
        double alpha, beta;
        int m,n,k;
        m = first_transverse_q;
        n = nphi_p;
        k = 1;
        alpha = -weight;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDA = first_transverse_q;
        LDB = nphi_p;
        LDC = ek.Rows();
        A = &div_phi(0,0);
        B = &phi_ps(0,0);
        C = &ek(0,first_p);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
    //contribucion matrix B^t
    {
        double *A, *B, *C;
        double alpha, beta;
        int m,n,k;
        m = nphi_p;
        n = first_transverse_q;
        k = 1;
        alpha = -weight;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDA = nphi_p;
        LDB = first_transverse_q;
        LDC = ek.Rows();
        A = &phi_ps(0,0);
        B = &div_phi(0,0);
        C = &ek(first_p,0);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
#else
    for (int iq = 0; iq < first_transverse_q; iq++) {
        //        datavec[qb].Print(std::cou        t);
        v_i = datavec[qb].fVecShapeIndex[iq].first;
        s_i = datavec[qb].fVecShapeIndex[iq].second;
        STATE kappa_inv_q_dot_phi_q_i = 0.0;
        for (int i = 0; i < 3; i++) {
            phi_q_i(i,0) = phi_qs(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
            kappa_inv_q_dot_phi_q_i        += kappa_inv_q(i,0)*phi_q_i(i,0);
        }
        
        ef(iq + first_q) += weight * ( kappa_inv_q_dot_phi_q_i - p * div_phi(iq,0));
        
        for (int jq = 0; jq < first_transverse_q; jq++) {
            
            v_j = datavec[qb].fVecShapeIndex[jq].first;
            s_j = datavec[qb].fVecShapeIndex[jq].second;
            //            if(v_j < nvecs-2){
            kappa_inv_phi_q_j.Zero();
            
            for (int j = 0; j < 3; j++) {
                kappa_inv_phi_q_j(j,0) += (1.0/eps)* inv_perm * phi_qs(s_j,0) * datavec[qb].fDeformedDirections(j,v_j);
            }
            
            STATE kappa_inv_phi_q_j_dot_phi_q_i = 0.0;
            for (int j = 0; j < 3; j++) {
                kappa_inv_phi_q_j_dot_phi_q_i +=  kappa_inv_phi_q_j(j,0)*phi_q_i(j,0);
            }
            
            ek(iq + first_q,jq + first_q) += weight * kappa_inv_phi_q_j_dot_phi_q_i;
        }
        if(v_j >=  nvecs-2){
            DebugStop();
        }
        for (int jp = 0; jp < nphi_p; jp++) {
            ek(iq + first_q, jp + first_p) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);
            ek(jp + first_p, iq + first_q) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);
        }
    }
#endif
    
    
//    
    // compute the contribution of the hdivbound
    for (int iq = first_transverse_q; iq < second_transverse_q; iq++)
    {
        //        datavec[qb].Print(std::cou        t);
        v_i = datavec[qb].fVecShapeIndex[iq].first;
        s_i = datavec[qb].fVecShapeIndex[iq].second;
        STATE kappa_inv_q_dot_phi_q_i = 0.0;

        for (int i = 0; i < 3; i++) {
            kappa_inv_q(i,0) += (2.0/kappaNormal)*q[i];
        }
        
        for (int i = 0; i < 3; i++) {
            phi_q_i(i,0) = phi_qs(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
            kappa_inv_q_dot_phi_q_i        += kappa_inv_q(i,0)*phi_q_i(i,0);
        }

        ef(iq + first_q) += weight * ( kappa_inv_q_dot_phi_q_i - p * div_phi(iq,0));
        for (int jq = first_transverse_q; jq < second_transverse_q; jq++)
        {

            v_j = datavec[qb].fVecShapeIndex[jq].first;
            s_j = datavec[qb].fVecShapeIndex[jq].second;
//            if(v_j != v_i) DebugStop();
            kappa_inv_phi_q_j.Zero();

            // kappanormal is the orthogonal permeability
            // phi_qs is the scalar function corresponding to jq (scalar)
            // s_j is the index of the scalar function corresponding to jq
            // v_j is the index of the vector corresponding to jq
            // kappa_inv_phi_q_j = vector equal to 1/kappa phi_qs * phi_qs * vector(v_j)
            // ad = 1
            // factor = 1
            for (int j = 0; j < 3; j++) {
                REAL KappaInvVal = (2.0/kappaNormal);
                kappa_inv_phi_q_j(j,0) = ad*KappaInvVal * factor * phi_qs(s_j,0) * datavec[qb].fDeformedDirections(j,v_j);
            }
            

            // kappa_inv_phi_q_j_dot_phi_q_i is the dot product of both vectors
            STATE kappa_inv_phi_q_j_dot_phi_q_i = 0.0;
            for (int j = 0; j < 3; j++) {
                kappa_inv_phi_q_j_dot_phi_q_i += kappa_inv_phi_q_j(j,0)*phi_q_i(j,0);
            }
            ek(iq + first_q,jq + first_q) += weight * kappa_inv_phi_q_j_dot_phi_q_i;
        }
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(iq + first_q, jp + first_p) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);
            ek(jp + first_p, iq + first_q) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);

        }
    }
//    // compute the contribution of the hdivbound
    for (int iq = second_transverse_q; iq < nphi_q; iq++)
    {
        //        datavec[qb].Print(std::cou        t);
        v_i = datavec[qb].fVecShapeIndex[iq].first;
        s_i = datavec[qb].fVecShapeIndex[iq].second;
        STATE kappa_inv_q_dot_phi_q_i = 0.0;

        for (int i = 0; i < 3; i++) {
            kappa_inv_q(i,0) += (2.0/kappaNormal)*q[i];
        }
        
        for (int i = 0; i < 3; i++) {
            phi_q_i(i,0) = phi_qs(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
            kappa_inv_q_dot_phi_q_i        += kappa_inv_q(i,0)*phi_q_i(i,0);
        }
        
        ef(iq + first_q) += weight * ( kappa_inv_q_dot_phi_q_i - p * div_phi(iq,0));
        for (int jq = second_transverse_q; jq < nphi_q; jq++)
        {
            v_j = datavec[qb].fVecShapeIndex[jq].first;
            s_j = datavec[qb].fVecShapeIndex[jq].second;
//            if(v_j != v_i) DebugStop();
            kappa_inv_phi_q_j.Zero();

            // kappanormal is the orthogonal permeability
            // phi_qs is the scalar function corresponding to jq (scalar)
            // s_j is the index of the scalar function corresponding to jq
            // v_j is the index of the vector corresponding to jq
            // kappa_inv_phi_q_j = vector equal to 1/kappa phi_qs * phi_qs * vector(v_j)
            for (int j = 0; j < 3; j++) {
                REAL KappaInvVal = (2.0/kappaNormal);
                kappa_inv_phi_q_j(j,0) += ad*KappaInvVal* factor *phi_qs(s_j,0) * datavec[qb].fDeformedDirections(j,v_j);
            }
            
            // kappa_inv_phi_q_j_dot_phi_q_i is the dot product of both vectors
            STATE kappa_inv_phi_q_j_dot_phi_q_i = 0.0;
            for (int j = 0; j < 3; j++) {
                kappa_inv_phi_q_j_dot_phi_q_i += kappa_inv_phi_q_j(j,0)*phi_q_i(j,0);
            }
            ek(iq + first_q,jq + first_q) += weight * kappa_inv_phi_q_j_dot_phi_q_i;

        }
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(iq + first_q, jp + first_p) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);
            ek(jp + first_p, iq + first_q) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);

        }
    }
    
  
    // Will we bring four spaces to pz?
//    if(this->mSimData.mTNumerics.m_four_approx_spaces_Q){
//        ContributeFourSpaces(datavec,weight,ek,ef);
//    }
}


void TPZMixedDarcyFractureFlow::ContributeFourSpaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    int qb = 0;
    int pb = 1;
    int g_avgb = 2;
    int p_avgb = 3;
    
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    int nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    
    int nphi_gb = datavec[g_avgb].phi.Rows();
    int nphi_pb = datavec[p_avgb].phi.Rows();
    if(nphi_q+nphi_p+nphi_gb+nphi_pb != ek.Rows())
    {
        DebugStop();
    }
    
    STATE p     = datavec[pb].sol[0][0];
    STATE g_avg = datavec[g_avgb].sol[0][0];
    STATE p_avg = datavec[p_avgb].sol[0][0];
    
    for(int ip=0; ip<nphi_p; ip++)
    {
        ef(nphi_q+ip,0) += weight * g_avg * phi_ps(ip,0);
        ek(nphi_q+ip,nphi_q+nphi_p) += weight * phi_ps(ip,0);
        
        ek(nphi_q+nphi_p,nphi_q+ip) += weight * phi_ps(ip,0);
    }
    
    ef(nphi_q+nphi_p+1,0) += -weight * g_avg;
    ek(nphi_q+nphi_p+1,nphi_q+nphi_p) += -weight;
    
    ef(nphi_q+nphi_p,0) += weight * (p - p_avg);
    ek(nphi_q+nphi_p,nphi_q+nphi_p+1) += -weight;
    
}


void TPZMixedDarcyFractureFlow::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    DebugStop();
    this->Contribute(datavec, weight, ef);
    
}


void TPZMixedDarcyFractureFlow::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->ContributeBC(datavec, weight, ekfake, ef, bc);
}


void TPZMixedDarcyFractureFlow::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    
    
    REAL gBigNumber = 1.0e12; //TPZMaterial::gBigNumber;
    
    int qb = 0;
    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].phi;
    
    int nphi_q       = phi_qs.Rows();
    int first_q      = 0;
    
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    
    TPZManVector<STATE,1> bc_data(1,0.0);
    bc_data[0] = bc.Val2()[0];
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD
        {
            STATE p_D = bc_data[0];
            for (int iq = 0; iq < nphi_q; iq++)
            {
                ef(iq + first_q) += weight * p_D * phi_qs(iq,0);
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN
        {
            
            for (int iq = 0; iq < nphi_q; iq++)
            {
                REAL qn_N = bc_data[0];
                REAL qn = 0.0;
                qn = q[0];
                
                ef(iq + first_q) += weight * gBigNumber * (qn - qn_N) * phi_qs(iq,0);
                for (int jq = 0; jq < nphi_q; jq++)
                {
                    ek(iq + first_q,jq + first_q) += weight * gBigNumber * phi_qs(jq,0) * phi_qs(iq,0);
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

























