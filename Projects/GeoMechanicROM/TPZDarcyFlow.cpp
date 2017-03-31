//
//  TPZDarcyFlow.cpp
//  PZ
//
//  Created by Omar on 3/5/17.
//
//

#include "TPZDarcyFlow.h"

/** @brief Default constructor */
TPZDarcyFlow::TPZDarcyFlow() : TPZMatWithMem<TPZDarcyFlowMemory, TPZDiscontinuousGalerkin>(){
    
    fdimension = 0;
    fCsymetric = -1.0;
    fIsMixedQ = false;
    
}

/** @brief Constructor based on a material id */
TPZDarcyFlow::TPZDarcyFlow(int matid, int dimension) : TPZMatWithMem<TPZDarcyFlowMemory, TPZDiscontinuousGalerkin>(matid){
    
    fdimension = dimension;
    fCsymetric = -1.0;
    fIsMixedQ = false;
}

/** @brief Constructor based on a Biot Poroelasticity  object */
TPZDarcyFlow::TPZDarcyFlow(const TPZDarcyFlow &mat) : TPZMatWithMem<TPZDarcyFlowMemory, TPZDiscontinuousGalerkin>(mat){
    
    this->fdimension    = mat.fdimension;
    this->fCsymetric    = mat.fCsymetric;
    this->fIsMixedQ     = mat.fIsMixedQ;
}


/** @brief Default destructor */
TPZDarcyFlow::~TPZDarcyFlow(){
    
}

/** @brief Set the required data at each integration point */
void TPZDarcyFlow::FillDataRequirements(TPZVec<TPZMaterialData> &datavec){
    
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = false;
    }
    
    if (fSimulationData->IsMixedApproxQ()) {
        datavec[0].fNeedsSol = true;
    }
    
}

/** @brief Set the required data at each integration point */
void TPZDarcyFlow::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec){
    
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsBasis = false;
        datavec[idata].fNeedsSol = false;
        datavec[idata].fNeedsNormal = false;
    }
    
    if (fSimulationData->IsMixedApproxQ()) {
        datavec[0].fNeedsSol = true;
    }
    
}


/** print out the data associated with the material */
void TPZDarcyFlow::Print(std::ostream &out){
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

/** returns the variable index associated with the name */
int TPZDarcyFlow::VariableIndex(const std::string &name){
    
    //	Diffusion Variables
    if(!strcmp("p_ex",name.c_str()))			return	0;
    if(!strcmp("v",name.c_str()))				return	1;
    
    //	Diffusion Variables
    if(!strcmp("pe_ex",name.c_str()))			return	2;
    if(!strcmp("ve",name.c_str()))				return	3;
    
    return -1;
    
}

/** returns the number of variables associated with the variable
 indexed by var.  var is obtained by calling VariableIndex */
int TPZDarcyFlow::NSolutionVariables(int var){
    
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return fdimension; // Vector
        case 2:
            return 1; // Scalar
        case 3:
            return fdimension; // Vector
    }
    return TPZMatWithMem::NSolutionVariables(var);
    
}

void TPZDarcyFlow::Compute_Sigma_qin(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_u){
    
    REAL trace;
    for (int i = 0; i < 3; i++) {
        trace = 0.0;
        for (int j = 0; j < 3; j++) {
            S(i,j) = fmu_quase_in * (Grad_u(i,j) + Grad_u(j,i));
            trace +=  Grad_u(j,j);
        }
        S(i,i) += flambda_quase_in * trace;
    }
    
    
}

void TPZDarcyFlow::Compute_Sigma(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_u){
    
    REAL trace;
    for (int i = 0; i < 3; i++) {
        trace = 0.0;
        for (int j = 0; j < 3; j++) {
            S(i,j) = fmu * (Grad_u(i,j) + Grad_u(j,i));
            trace +=  Grad_u(j,j);
        }
        S(i,i) += flambda * trace;
    }
    
    
}

void TPZDarcyFlow::ContributeUndrained(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    if (fIsMixedQ) {
        ContributeUndrainedMF(datavec, weight, ek, ef);
        return;
    }
    
    int p_b = 0;
    
    // Getting the space functions from memory
    long global_point_index = datavec[p_b].intGlobPtIndex;
    TPZDarcyFlowMemory &point_memory = GetMemory()[global_point_index];
    
    TPZFMatrix<REAL>    &phip   =   point_memory.phi_p();
    REAL   & p_n = point_memory.p_n();
    TPZFMatrix<REAL>    &grad_u_n   =   point_memory.grad_u_n();
    
    int nphi_p = phip.Rows();
    int first_p = 0;
    
    TPZFNMatrix<9,REAL> S(3,3),S_n(3,3);
    Compute_Sigma_qin(S_n, grad_u_n);
    REAL S_n_v = (S_n(0,0) + S_n(1,1) + S_n(2,2))/3.0;

    point_memory.Set_sigma_0(S_n);
    
    if (!fSimulationData->IsCurrentStateQ()) {
        
        return;
    }
    

    for (int ip = 0; ip < nphi_p; ip++) {
        
        ef(ip + first_p, 0)		+= weight * ( (p_n + S_n_v) * phip(ip,0) );
        
        for (int jp = 0; jp < nphi_p; jp++) {
            
            ek(ip + first_p, jp + first_p)  += weight * phip(jp,0) * phip(ip,0) ;
        }
        
    }
    
}

void TPZDarcyFlow::ContributeUndrainedMF(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int q_b = 0;
    int p_b = 1;
    
    // Getting the space functions from memory
    long global_point_index = datavec[q_b].intGlobPtIndex;
    TPZDarcyFlowMemory &point_memory = GetMemory()[global_point_index];
    
    TPZFMatrix<REAL>    &phip   =   point_memory.phi_p();
    REAL   & p_n = point_memory.p_n();
    TPZFMatrix<REAL>    &grad_u_n   =   point_memory.grad_u_n();
    
    // Getting the space functions
    TPZFMatrix<REAL>    &phiq   =   datavec[q_b].phi;
//    TPZFMatrix<REAL>    &phip   =   datavec[p_b].phi;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,3> q = datavec[q_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix<10,STATE> dq = datavec[q_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    // Size of matrix blocks
    int nphi_q = datavec[q_b].fVecShapeIndex.NElements();
    int nphi_p = phip.Rows();
    
    int first_q = 0;
    int first_p = first_q + nphi_q;

    
    TPZFNMatrix<9,REAL> S(3,3),S_n(3,3);
    Compute_Sigma(S_n, grad_u_n);
    REAL S_n_v = (S_n(0,0) + S_n(1,1) + S_n(2,2))/3.0;
    
    if (!fSimulationData->IsCurrentStateQ()) {
        
        return;
    }
    
    REAL c = fk/feta;
    TPZFNMatrix<3,STATE> phi_q_i(3,1), phi_q_j(3,1);
    int s_i, s_j;
    int v_i, v_j;
    REAL jac_det = datavec[q_b].detjac;
    TPZFNMatrix<40,STATE> div_on_master;
    REAL div_qq;
    this->ComputeDivergenceOnMaster(datavec, div_on_master,div_qq);
    
    for (int iq = 0; iq < nphi_q; iq++)
    {
        
        v_i = datavec[q_b].fVecShapeIndex[iq].first;
        s_i = datavec[q_b].fVecShapeIndex[iq].second;
        
        REAL Kl_inv_dot_q = 0.0;
        for (int k = 0; k < q.size(); k++) {
            phi_q_i(k,0) = phiq(s_i,0) * datavec[q_b].fNormalVec(k,v_i);
            Kl_inv_dot_q        += (1.0/c)*q[k]*phi_q_i(k,0);
        }
        
//        ef(iq + first_q) +=  weight * ( Kl_inv_dot_q - (1.0/jac_det) * p_n * div_on_master(iq,0));
        ef(iq + first_q) +=  weight * ( 0.0 );
        
        for (int jq = 0; jq < nphi_q; jq++)
        {
            
            v_j = datavec[q_b].fVecShapeIndex[jq].first;
            s_j = datavec[q_b].fVecShapeIndex[jq].second;
            
            STATE Kl_inv_phi_u_j_dot_phi_u_j = 0.0;
            for (int k = 0; k < q.size(); k++) {
                phi_q_j(k,0) = phiq(s_j,0) * datavec[q_b].fNormalVec(k,v_j);
                Kl_inv_phi_u_j_dot_phi_u_j += (1.0/c)*phi_q_j(k,0)*phi_q_i(k,0);
            }
            
            ek(iq + first_q,jq + first_q) +=  weight * Kl_inv_phi_u_j_dot_phi_u_j;
        }

        
    }

    
    for (int ip = 0; ip < nphi_p; ip++) {
        
        ef(ip + first_p, 0)		+= weight * ( (p_n + S_n_v) * phip(ip,0) );
        
        for (int jp = 0; jp < nphi_p; jp++) {
            
            ek(ip + first_p, jp + first_p)  += weight * phip(jp,0) * phip(ip,0) ;
        }
        
    }
    
}

// Contribute Methods being used
void TPZDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    if (fSimulationData->IsInitialStateQ()) {
        ContributeUndrained(datavec, weight, ek, ef);
        return;
    }
    
    if (fIsMixedQ) {
        ContributeMF(datavec, weight, ek, ef);
        return;
    }
    
    int p_b = 0;
    
    // Getting the space functions from memory
    long global_point_index = datavec[p_b].intGlobPtIndex;
    TPZDarcyFlowMemory &point_memory = GetMemory()[global_point_index];
    
    TPZFMatrix<REAL>    &phip   =   point_memory.phi_p();
    TPZFMatrix<REAL>    &grad_phi_p   =   point_memory.grad_phi_p();
    
    REAL   & p_n = point_memory.p_n();
    REAL   & p   = point_memory.p();
    
    // initial state
    REAL   & p_0   = point_memory.p_0();
    TPZFMatrix<REAL>    &S_0 = point_memory.sigma_0();
    
    TPZFMatrix<REAL>    &grad_u   =   point_memory.grad_u();
    TPZFMatrix<REAL>    &grad_u_n   =   point_memory.grad_u_n();
    
    TPZFMatrix<REAL>    &grad_p = point_memory.grad_p_n();
    
    int nphi_p = phip.Rows();
    int first_p = 0;
    
    REAL dt = fSimulationData->dt();
    REAL alpha = 1.0;
    REAL Se = 0.0;
    
    TPZFNMatrix<9,REAL> S(3,3),S_n(3,3);
    Compute_Sigma(S, grad_u);
    Compute_Sigma(S_n, grad_u_n);


    REAL div_u = grad_u(0,0) + grad_u(1,1);
    REAL div_u_n = grad_u_n(0,0) + grad_u_n(1,1);
    REAL lambda     = 8.333e3;
    REAL mu         = 12.50e3;
    REAL Kdr = lambda + (2.0/3.0)*mu;
    
    REAL S_v_0 = (S_0(0,0) + S_0(1,1) + S_0(2,2))/3.0;
    REAL S_v = (S(0,0) + S(1,1) + S(2,2))/3.0;
    REAL S_n_v = (S_n(0,0) + S_n(1,1) + S_n(2,2))/3.0;
    REAL Ss = (Se + alpha*alpha/Kdr);
    
    REAL phi = fphi_0 + alpha * (S_v - S_v_0) / Kdr + Ss * (p - p_0);
    REAL phi_n = fphi_0 + alpha * (S_n_v - S_v_0) / Kdr + Ss * (p_n - p_0);
    
//    REAL phi = fphi_0 + alpha * (S_n_v) / Kdr + Ss * (p);
//    REAL phi_n = fphi_0 + alpha * (S_n_v) / Kdr + Ss * (p_n);
    
    REAL ct = 1.0;
    
    if (!fSimulationData->IsCurrentStateQ()) {
        
        // Darcy mono-phascis flow
        for (int ip = 0; ip < nphi_p; ip++) {
            if(fSimulationData->CheckCouplingQ()){
                ef(ip + first_p, 0)		+=  weight *  (-1.0) * (1.0/dt) * (alpha * div_u + Se * p ) * phip(ip,0);
            }
            else{
                ef(ip + first_p, 0)		+=  weight *  (-1.0) * (1.0/dt) * (ct * phi) * phip(ip,0);
            }

            
        }
        return;
    }

    
    REAL c = fk/feta;
    // Darcy mono-phasic flow
    for (int ip = 0; ip < nphi_p; ip++) {
        
        REAL dot = 0.0;
        for (int i = 0;  i < fdimension; i++) {
            dot += grad_p(i,0) * grad_phi_p(i,ip);
        }
        
        if(fSimulationData->CheckCouplingQ()){
            
            ef(ip + first_p, 0)		+= weight * (c * dot + (1.0/dt) * (alpha * div_u_n + Se * p_n) * phip(ip,0) );
            
            for (int jp = 0; jp < nphi_p; jp++) {
                
                REAL dot = 0.0;
                for (int i = 0;  i < fdimension; i++) {
                    dot += grad_phi_p(i,jp) * grad_phi_p(i,ip);
                }
                
                ek(ip + first_p, jp + first_p)  += weight * ( c * dot + (1.0/dt) * (Se * phip(jp,0) ) * phip(ip,0) );
            }
            
        }
        else{
            
            ef(ip + first_p, 0)		+= weight * (c * dot + (1.0/dt) * (ct * phi_n) * phip(ip,0) );
            
            for (int jp = 0; jp < nphi_p; jp++) {
                
                REAL dot = 0.0;
                for (int i = 0;  i < fdimension; i++) {
                    dot += grad_phi_p(i,jp) * grad_phi_p(i,ip);
                }
                
                ek(ip + first_p, jp + first_p)  += weight * ( c * dot + (1.0/dt) * (ct * Ss * phip(jp,0) ) * phip(ip,0) );
            }
            
        }
        
    }
    
}


void TPZDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
}

void TPZDarcyFlow::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (fSimulationData->IsInitialStateQ()) {
        return;
    }
    
    if (fIsMixedQ) {
        this->ContributeMFBC(datavec, weight, ek, ef, bc);
        return;
    }
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int p_b = 0;
    
    // Getting the solutions and derivatives
    // Get the data at the integrations points
    TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond>  & material_mem = dynamic_cast<TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond > & >(bc);
    
    long global_point_index = datavec[p_b].intGlobPtIndex;
    TPZDarcyFlowMemory &point_memory = material_mem.GetMemory()[global_point_index];
    
    TPZFMatrix<REAL>  & phip = point_memory.phi_p();
    REAL p_n = point_memory.p_n();
    
    TPZManVector<REAL,1> p(1);
    p[0] = p_n;
    
    int phrp = phip.Rows();
    short in,jn;
    REAL v[1];
    v[0] = bc.Val2()(0,0);	//	Pressure or Flux
    
    REAL time = this->SimulationData()->t();
    if (bc.HasTimedependentBCForcingFunction()) {
        TPZManVector<REAL,3> f(1);
        TPZFMatrix<REAL> gradf;
        bc.TimedependentBCForcingFunction()->Execute(datavec[p_b].x, time, f, gradf);
        v[0] = f[0];	//	Pressure or flux
    }
    
    // Dirichlet in Pressure
    switch (bc.Type())
    {
        case 0 :
        {
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in,0)		+= gBigNumber*(p[0]-v[0])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in,jn)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
                }
            }
            break;
        }
            
        case 1 :
        {
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in,0)	+= -1.0 * weight * v[0] * phip(in,0);	// Qnormal
            }
            break;
        }
                        
        default:
        {
            DebugStop();
        }
            break;
    }
    
}


/** returns the solution associated with the var index based on the finite element approximation */
void TPZDarcyFlow::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
    
    if (fIsMixedQ) {
        this->SolutionMF(datavec, var, Solout);
        return;
    }
    
    Solout.Resize( this->NSolutionVariables(var));
    int p_b = 0;
    
    // Getting the space functions
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    
    REAL to_Mpa     = 1.0;//1.0e-6;
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_p(3,1,0.0),Grad_u(3,3,0.0),Grad_u_n(3,3,0.0),e_e(3,3,0.0),e_p(3,3,0.0),S;
    
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0); // dp/dx
    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1); // dp/dy
    
    //	p
    if(var == 0){
        Solout[0] = p[0]*to_Mpa;
        return;
    }
    
    //	darcy
    if(var == 1) {
        Solout[0] = -1.0*(fk/feta)*Grad_p(0,0);
        Solout[1] = -1.0*(fk/feta)*Grad_p(1,0);
        return;
    }
    
    //	p exact
    if(var == 2) {
        TPZManVector<STATE,5> f(1,0.0);
        TPZFNMatrix<4,STATE> df(4,1,0.0);
        if (this->HasTimedependentForcingFunction()) {
            REAL time = fSimulationData->t();
            this->fTimeDependentForcingFunction->Execute(datavec[p_b].x, time, f, df);
        }
        
        Solout[0] = f[2]*to_Mpa;
        
        return;
    }
    
    //	darcy exact
    if(var == 3) {
        TPZManVector<STATE,5> f(1,0.0);
        TPZFNMatrix<4,STATE> df(4,1,0.0);
        if (this->HasTimedependentForcingFunction()) {
            REAL time = fSimulationData->t();
            this->fTimeDependentForcingFunction->Execute(datavec[p_b].x, time, f, df);
        }
        
        Solout[0] = f[3];
        Solout[1] = f[4];
        return;
    }

    
    std::cout  << "not implemented. " << std::endl;
    DebugStop();
    
}


///////////////////////////////////// monophasic mixed methods ////////////////////////////////////////////////


/** Computes the divergence over the parametric space */
void TPZDarcyFlow::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &Divergence_of_q){
    
    int qblock = 0;
    int dim = this->Dimension();
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[qblock].phi;   // For H1  test functions Q
    TPZFMatrix<STATE> dphiuH1       = datavec[qblock].dphi; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiuH1axes   = datavec[qblock].dphix; // Derivative For H1  test functions
    TPZFNMatrix<9,STATE> gradu = datavec[qblock].dsol[0];
    TPZFNMatrix<9,STATE> graduMaster;
    gradu.Transpose();
    
    TPZFNMatrix<660> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[qblock].axes);
    
    int nphiuHdiv = datavec[qblock].fVecShapeIndex.NElements();
    
    DivergenceofPhi.Resize(nphiuHdiv,1);
    
    REAL JacobianDet = datavec[qblock].detjac;
    
    TPZFMatrix<STATE> Qaxes = datavec[qblock].axes;
    TPZFMatrix<STATE> QaxesT;
    TPZFMatrix<STATE> Jacobian = datavec[qblock].jacobian;
    TPZFMatrix<STATE> JacobianInverse = datavec[qblock].jacinv;
    
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
            ivectorindex = datavec[qblock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[qblock].fVecShapeIndex[iq].second;
            
            for (int k = 0; k < dim; k++) {
                VectorOnXYZ(k,0) = datavec[qblock].fNormalVec(k,ivectorindex);
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
            Divergence_of_q += graduMaster(k,k);
        }
        
    }
    else
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[qblock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[qblock].fVecShapeIndex[iq].second;
            
            /* Computing the divergence for constant jacobian elements */
            for (int k = 0; k < dim; k++) {
                DivergenceofPhi(iq,0) +=  datavec[qblock].fNormalVec(k,ivectorindex)*GradphiuH1(k,ishapeindex);
            }
        }
    }
    
    return;
}

void TPZDarcyFlow::ContributeMF(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int q_b = 0;
    int p_b = 1;
    
    long global_point_index = datavec[q_b].intGlobPtIndex;
    TPZDarcyFlowMemory &point_memory = GetMemory()[global_point_index];
    
    TPZFMatrix<REAL>    &phip   =   point_memory.phi_p();
    TPZFMatrix<REAL>    &grad_phi_p   =   point_memory.grad_phi_p();
    
    REAL   & p_n = point_memory.p_n();
    REAL   & p   = point_memory.p();
    
    TPZFMatrix<REAL>    &grad_u   =   point_memory.grad_u();
    TPZFMatrix<REAL>    &grad_u_n   =   point_memory.grad_u_n();
    
    TPZFMatrix<REAL>    &grad_p = point_memory.grad_p_n();
    
    // Getting the space functions
    TPZFMatrix<REAL>    &phiq   =   datavec[q_b].phi;
//    TPZFMatrix<REAL>    &phip   =   datavec[p_b].phi;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,3> q = datavec[q_b].sol[0];
//    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix<10,STATE> dq = datavec[q_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    // Size of matrix blocks
    int nphi_q = datavec[q_b].fVecShapeIndex.NElements();
    int nphi_p = phip.Rows();
    
    int first_q = 0;
    int first_p = first_q + nphi_q;
    
    REAL dt = fSimulationData->dt();
    REAL alpha = 1.0;
    REAL Se = 0.0;
    
    TPZFNMatrix<9,REAL> S(3,3),S_n(3,3);
    Compute_Sigma(S, grad_u_n);
    Compute_Sigma(S_n, grad_u_n);
    
    // Getting memory information
    
    REAL lambda     = 8.333e3;
    REAL mu         = 12.50e3;
    REAL Kdr = lambda + (2.0/3.0)*mu;
    REAL S_v = (S(0,0) + S(1,1) + S(2,2))/3.0;
    REAL S_n_v = (S_n(0,0) + S_n(1,1) + S_n(2,2))/3.0;
    REAL Ss = (Se + alpha*alpha/Kdr);
    S_v = S_n_v;
    
    if (!fSimulationData->IsCurrentStateQ()) {
        
        // Darcy mono-phascis flow
        for (int ip = 0; ip < nphi_p; ip++) {
            ef(ip + first_p, 0)		+= (fCsymetric) * weight *  (-1.0) * (1.0/dt) * ((alpha * S_v / Kdr + Ss * p )) * phip(ip,0);
        }
        return;
        
    }
    
    REAL c = fk/feta;
    TPZFNMatrix<3,STATE> phi_q_i(3,1), phi_q_j(3,1);
    int s_i, s_j;
    int v_i, v_j;
    REAL jac_det = datavec[q_b].detjac;
    TPZFNMatrix<40,STATE> div_on_master;
    REAL div_qq;
    this->ComputeDivergenceOnMaster(datavec, div_on_master,div_qq);
    
    for (int iq = 0; iq < nphi_q; iq++)
    {
        
        v_i = datavec[q_b].fVecShapeIndex[iq].first;
        s_i = datavec[q_b].fVecShapeIndex[iq].second;
        
        REAL Kl_inv_dot_q = 0.0;
        for (int k = 0; k < q.size(); k++) {
            phi_q_i(k,0) = phiq(s_i,0) * datavec[q_b].fNormalVec(k,v_i);
            Kl_inv_dot_q        += (1.0/c)*q[k]*phi_q_i(k,0);
        }
        
        ef(iq + first_q) +=  weight * ( Kl_inv_dot_q - (1.0/jac_det) * p_n * div_on_master(iq,0));
        
        for (int jq = 0; jq < nphi_q; jq++)
        {
            
            v_j = datavec[q_b].fVecShapeIndex[jq].first;
            s_j = datavec[q_b].fVecShapeIndex[jq].second;
            
            STATE Kl_inv_phi_u_j_dot_phi_u_j = 0.0;
            for (int k = 0; k < q.size(); k++) {
                phi_q_j(k,0) = phiq(s_j,0) * datavec[q_b].fNormalVec(k,v_j);
                Kl_inv_phi_u_j_dot_phi_u_j += (1.0/c)*phi_q_j(k,0)*phi_q_i(k,0);
            }
            
            ek(iq + first_q,jq + first_q) +=  weight * Kl_inv_phi_u_j_dot_phi_u_j;
        }
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(iq + first_q, jp + first_p) +=  weight * ( - (1.0/jac_det) * phip(jp,0) * div_on_master(iq,0) ) ;
        }
        
    }
    
    
    REAL div_q = (dq(0,0) + dq(1,1) + dq(2,2));
    // Darcy mixed mono-phasic flow restriction equation
    for (int ip = 0; ip < nphi_p; ip++) {
        
        
        ef(ip + first_p, 0)		+= (fCsymetric) * weight * (div_q + (1.0/dt) * ((alpha * S_v / Kdr + Ss * p_n))) * phip(ip,0);
        
        for (int jq = 0; jq < nphi_q; jq++)
        {
            ek(ip + first_p, jq + first_q) += (fCsymetric) * weight * (1.0/jac_det) * div_on_master(jq,0) * phip(ip,0);
        }
        
        for (int jp = 0; jp < nphi_p; jp++) {
            
            ek(ip + first_p, jp + first_p)  += (fCsymetric) * weight * ( (1.0/dt) * Ss * phip(jp,0) * phip(ip,0) );
        }
        
    }
    
    
}

void TPZDarcyFlow::ContributeMF(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->ContributeMF(datavec, weight, ek_fake, ef);
}

void TPZDarcyFlow::ContributeMFBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int q_b = 0;
    int p_b = 1;
    
    TPZFMatrix<REAL>  &phiq = datavec[q_b].phi;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,3> q = datavec[q_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    int phrq = phiq.Rows();
    
    int firstq = 0;
    
    REAL v[1];
    
    v[0] = bc.Val2()(0,0);	//	Pressure or flux
    
    REAL time = this->SimulationData()->t();
    REAL dt  = this->SimulationData()->dt();
    REAL Value = bc.Val2()(0,0);
    if (bc.HasTimedependentBCForcingFunction()) {
        TPZManVector<REAL,3> f(3);
        TPZFMatrix<REAL> gradf;
        bc.TimedependentBCForcingFunction()->Execute(datavec[p_b].x, time, f, gradf);
        v[0] = f[0];	//	Pressure
    }
    else{
        Value = bc.Val2()(0,0);
    }
    
    // Dirichlet in Pressure
    switch (bc.Type())
    {
        case 0 :
        {
            
            //	Diffusion Equation Pressure data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) +=  weight * v[0] * phiq(iq,0);
            }
            
            break;
        }
            
        case 1 :
        {
            
            //	Diffusion Equation Neumann data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += weight * gBigNumber * (q[0] - v[0]) * phiq(iq,0);
                
                for (int jq = 0; jq < phrq; jq++)
                {
                    ek(iq + firstq,jq + firstq) += weight * gBigNumber * phiq(jq,0) * phiq(iq,0);
                }
                
            }
            break;
        }
            
        default:
        {
            DebugStop();
        }
            break;
    }
    
}

/** returns the solution associated with the var index based on the finite element approximation */
void TPZDarcyFlow::SolutionMF(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    
    int q_b = 0;
    int p_b = 1;
    
    // Getting the space functions
    
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> q = datavec[q_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    
    REAL to_Mpa     = 1.0;//1.0e-6;
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_p(3,1,0.0),Grad_u(3,3,0.0),Grad_u_n(3,3,0.0),e_e(3,3,0.0),e_p(3,3,0.0),S;
    
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0); // dp/dx
    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1); // dp/dy
    

    //	pressure
    if(var == 0) {
        Solout[0] = p[0]*to_Mpa;
        return;
    }
    
    //	darcy
    if(var == 1) {
        for (int d = 0; d <fdimension; d++) {
            Solout[d] = q[d];
        }
        return;
    }
    
    // Pressure exact
    if(var == 2) {
        
        TPZManVector<STATE,5> f(1,0.0);
        TPZFNMatrix<4,STATE> df(4,1,0.0);
        if (this->HasTimedependentForcingFunction()) {
            REAL time = fSimulationData->t();
            this->fTimeDependentForcingFunction->Execute(datavec[q_b].x, time, f, df);
        }
        
        Solout[0] = f[2]*to_Mpa;
        
        return;
    }
    
    //	darcy exact
    if(var == 3) {
        TPZManVector<STATE,5> f(1,0.0);
        TPZFNMatrix<4,STATE> df(4,1,0.0);
        if (this->HasTimedependentForcingFunction()) {
            REAL time = fSimulationData->t();
            this->fTimeDependentForcingFunction->Execute(datavec[q_b].x, time, f, df);
        }
        
        Solout[0] = f[3];
        Solout[1] = f[4];
        return;
    }
    
    
    
    std::cout  << "not implemented. " << std::endl;
    DebugStop();
    
}


/** @brief Unique identifier for serialization purposes */
int TPZDarcyFlow::ClassId() const {
    return -6378637866;
}

/** @brief Save object data to a stream */
void TPZDarcyFlow::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
}

/** @brief Read object data from a stream */
void TPZDarcyFlow::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}
