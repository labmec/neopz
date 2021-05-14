//
//  TPZDualPoisson.cpp
//  PZ
//
//  Created by omar on 04/07/2016.
//
//

#include "TPZDualPoisson.h"


TPZDualPoisson::TPZDualPoisson(){
    
}

TPZDualPoisson::~TPZDualPoisson(){
    
}

TPZDualPoisson::TPZDualPoisson(int mat_id): TPZMaterial(mat_id){
    
}

TPZDualPoisson::TPZDualPoisson(const TPZDualPoisson &copy) : TPZMaterial(copy){
    
}

TPZMaterial * TPZDualPoisson::NewMaterial(){
    return new TPZDualPoisson(*this);
}

TPZDualPoisson & TPZDualPoisson::operator=(const TPZDualPoisson &other){
    
    if (this != & other) // prevent self-assignment
    {
        TPZMaterial::operator=(other);
    }
    return *this;
}

int TPZDualPoisson::Dimension() const { return 3;}

int TPZDualPoisson::NStateVariables() const
{
    return 1;
}

void TPZDualPoisson::Print(std::ostream & out){
    TPZMaterial::Print(out);
}

std::string TPZDualPoisson::Name() { return "TPZDualPoisson"; }


void TPZDualPoisson::FillDataRequirements(TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
}

void TPZDualPoisson::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
}

void TPZDualPoisson::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

void TPZDualPoisson::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

int TPZDualPoisson::ClassId() const{
    return Hash("TPZDualPoisson") ^ TPZMaterial::ClassId() << 1;
}

void TPZDualPoisson::Write(TPZStream &buf, int withclassid) const{
    DebugStop();
}


void TPZDualPoisson::Read(TPZStream &buf, void *context){
    DebugStop();
}


/** @} */

/**
 * @name Contribute methods (weak formulation) no multiphysics mesh
 * @{
 */

/** @brief Volumetric contribute with jacobian matrix */
void TPZDualPoisson::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    DebugStop();
    
}


/** @brief Volumetric contribute */
void TPZDualPoisson::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef){
    
    DebugStop();
}


/** @brief Boundary contribute with jacobian matrix */
void TPZDualPoisson::ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    DebugStop();
    
}


/** @brief Boundary contribute */
void TPZDualPoisson::ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    DebugStop();
    
}


/** @} */


/**
 * @name Contribute methods (weak formulation) multiphysics mesh
 * @{
 */

/** @brief Volumetric contribute with jacobian matrix */
void TPZDualPoisson::Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    int ub = 0;
    int pb = 1;
    
    TPZFNMatrix<100,REAL> phi_us       = datavec[ub].phi;
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,REAL> dphi_us      = datavec[ub].dphix;
    TPZFNMatrix<100,REAL> dphi_ps      = datavec[pb].dphix;
    
   
	TPZFNMatrix<40, REAL> div_on_master = datavec[ub].divphi;
	STATE divu = datavec[ub].divsol[0][0];
    STATE divflux;
    REAL jac_det = datavec[ub].detjac;
    
    int nphiu       = datavec[ub].fVecShapeIndex.NElements();
    int nphip       = phi_ps.Rows();
    int firstu      = 0;
    int firstp      = nphiu + firstu;
    
    TPZManVector<STATE,3> u  = datavec[ub].sol[0];
    STATE p                  = datavec[pb].sol[0][0];
    
    TPZFNMatrix<10,STATE> Graduaxes = datavec[ub].dsol[0];
    
    TPZFNMatrix<3,STATE> phi_u_i(3,1), phi_u_j(3,1);
    
    int s_i, s_j;
    int v_i, v_j;
    
    for (int iu = 0; iu < nphiu; iu++)
    {
        
        v_i = datavec[ub].fVecShapeIndex[iu].first;
        s_i = datavec[ub].fVecShapeIndex[iu].second;
        
        STATE u_dot_phi_u_i = 0.0;
        for (int i = 0; i < u.size(); i++) {
            phi_u_i(i,0) = phi_us(s_i,0) * datavec[ub].fDeformedDirections(i,v_i);
            u_dot_phi_u_i        += u[i]*phi_u_i(i,0);
        }
        
        ef(iu + firstu) += weight * ( u_dot_phi_u_i - p * div_on_master(iu,0));
        
        for (int ju = 0; ju < nphiu; ju++)
        {
            
            v_j = datavec[ub].fVecShapeIndex[ju].first;
            s_j = datavec[ub].fVecShapeIndex[ju].second;
            
            STATE phi_u_j_dot_phi_u_i = 0.0;
            for (int j = 0; j < u.size(); j++) {
                phi_u_j(j,0) = phi_us(s_j,0) * datavec[ub].fDeformedDirections(j,v_j);
                phi_u_j_dot_phi_u_i += phi_u_j(j,0)*phi_u_i(j,0);
            }
            
            ek(iu + firstu,ju + firstu) += weight * phi_u_j_dot_phi_u_i;
        }
        
        for (int jp = 0; jp < nphip; jp++)
        {
            ek(iu + firstu, jp + firstp) += weight * ( - div_on_master(iu,0) ) * phi_ps(jp,0);
        }
        
    }
    
    TPZManVector<STATE,1> f(1,0.0);
    if (this->HasForcingFunction()) {
        this->fForcingFunction->Execute(datavec[ub].x, f);
    }
    
    for (int ip = 0; ip < nphip; ip++)
    {
        
        ef(ip + firstp) += -1.0 * weight * (divu - f[0]) * phi_ps(ip,0);
        
        for (int ju = 0; ju < nphiu; ju++)
        {
            ek(ip + firstp, ju + firstu) += -1.0 * weight * div_on_master(ju,0) * phi_ps(ip,0);
        }
        
    }
    
}

/** @brief Volumetric contribute */
void TPZDualPoisson::Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ekfake, ef);
}


/** @brief Boundary contribute without jacobian matrix */
void TPZDualPoisson::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->ContributeBC(datavec, weight, ekfake, ef, bc);
}


/** @brief Boundary contribute */
void TPZDualPoisson::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    int ub = 0;
    TPZFNMatrix<100,REAL> phi_us       = datavec[ub].phi;
    
    int nphiu       = phi_us.Rows();
    int firstu      = 0;
    
    TPZManVector<STATE,3> u  = datavec[ub].sol[0];
    
    TPZManVector<STATE,1> bc_data(1,0.0);
    bc_data[0] = bc.Val2()(0,0);
    if (bc.HasForcingFunction()) {
        bc.ForcingFunction()->Execute(datavec[ub].x, bc_data);
    }
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD
        {
            STATE p_D = bc_data[0];
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
                STATE un_N = bc_data[0], un = u[0];
                ef(iu + firstu) += weight * gBigNumber * (un - un_N) * phi_us(iu,0);
                
                for (int ju = 0; ju < nphiu; ju++)
                {
                    
                    ek(iu + firstu,ju + firstu) += weight * gBigNumber * phi_us(ju,0) * phi_us(iu,0);
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


/** @} */

/**
 * @name Post-processing methods
 * @{
 */


int TPZDualPoisson::NEvalErrors() {return 3;}

int TPZDualPoisson::VariableIndex(const std::string &name){
    
    if(!strcmp("q",name.c_str()))               return  1;
    if(!strcmp("p",name.c_str()))               return  2;
    if(!strcmp("q_exact",name.c_str()))         return  3;
    if(!strcmp("p_exact",name.c_str()))         return  4;
    if(!strcmp("f_exact",name.c_str()))         return  5;
    if(!strcmp("div_q",name.c_str()))           return  6;
    
    return TPZMaterial::VariableIndex(name);
}

int TPZDualPoisson::NSolutionVariables(int var){
    if(var == 1) return this->Dimension();
    if(var == 2) return 1;
    if(var == 3) return this->Dimension();
    if(var == 4) return 1;
    if(var == 5) return 1;
    if(var == 6) return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}

void TPZDualPoisson::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    DebugStop();
    
}

void TPZDualPoisson::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    int ub = 0;
    int pb = 1;
    
    Solout.Resize( this->NSolutionVariables(var));
    TPZManVector<STATE,3> p, u, f;
    
    u = datavec[ub].sol[0];
    p = datavec[pb].sol[0];
    
    STATE div_u = datavec[ub].divsol[0][0];
    
    if(var == 1){
        for (int i=0; i < this->Dimension(); i++)
        {
            Solout[i] = u[i];
        }
        return;
    }
    
    if(var == 2){
        Solout[0] = p[0];
        return;
    }
    
    if(var == 3){
        TPZManVector<STATE,1> f(1,0.0);
        TPZFNMatrix<4,STATE> df(4,1,0.0);
        if (this->HasExactSol()) {
            this->fExactSol->Execute(datavec[ub].x, f, df);
        }

        for (int i=0; i < this->Dimension(); i++)
        {
            Solout[i] = df(i,0);
        }
        return;
    }
    
    if(var == 4){
        TPZManVector<STATE,1> f(1,0.0);
        TPZFNMatrix<4,STATE> df(4,1,0.0);
        if (this->HasExactSol()) {
            this->fExactSol->Execute(datavec[ub].x, f, df);
        }
        Solout[0] = f[0];
        return;
    }
    
    if(var == 5){
        TPZManVector<STATE,1> f(1,0.0);
        TPZFNMatrix<4,STATE> df(4,1,0.0);
        if (this->HasExactSol()) {
            this->fExactSol->Execute(datavec[ub].x, f, df);
        }
        Solout[0] = df(3,0);
        return;
    }
    
    if(var == 6){
        Solout[0] = div_u;
        return;
    }
    
    
    DebugStop();
}

void TPZDualPoisson::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{

    errors.Fill(0.0);
    
    int ub = 0;
    int pb = 1;
    TPZManVector<STATE,3> p, u, f;
    
    u = data[ub].sol[0];
    p = data[pb].sol[0];
    
    TPZFMatrix<STATE> dudx = data[ub].dsol[0];
    
    STATE div_u = 0.0;
    for(int i = 0; i < this->Dimension(); i++) {
        div_u  += dudx(i,i);
    }
    
    STATE div_error = div_u - du_exact(3,0); // using f source term on the fourth position of du_exact
    
    /** @brief   error[0] : primal error using L2 norm */
    STATE p_error = p[0] - u_exact[0];
    errors[0]  = p_error*p_error;
    
    /** @brief   error[1] : dual error using L2 norm */
    for(int i = 0; i < this->Dimension(); i++) {
        STATE d_error = u[i] - du_exact(i,0);
        errors[1]  += d_error*d_error;
    }
    
    /** @brief   error[2] : div error norm */
    errors[2]= div_error * div_error;
}

/** @} */
