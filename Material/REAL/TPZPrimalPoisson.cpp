//
//  TPZPrimalPoisson.cpp
//  PZ
//
//  Created by omar on 04/07/2016.
//
//

#include "TPZPrimalPoisson.h"

TPZPrimalPoisson::TPZPrimalPoisson(){
    
}

TPZPrimalPoisson::~TPZPrimalPoisson(){
    
}

TPZPrimalPoisson::TPZPrimalPoisson(int mat_id): TPZMaterial(mat_id){
    
}

TPZPrimalPoisson::TPZPrimalPoisson(const TPZPrimalPoisson &copy){

}

TPZMaterial * TPZPrimalPoisson::NewMaterial(){
    return new TPZPrimalPoisson(*this);
}

TPZPrimalPoisson & TPZPrimalPoisson::operator=(const TPZPrimalPoisson &other){
    if (this != & other) // prevent self-assignment
    {
    }
    return *this;
}

int TPZPrimalPoisson::Dimension() const { return 3;}

int TPZPrimalPoisson::NStateVariables() {return 1;}

void TPZPrimalPoisson::Print(std::ostream & out){
    TPZMaterial::Print(out);
}

std::string TPZPrimalPoisson::Name() { return "TPZPrimalPoisson"; }


void TPZPrimalPoisson::FillDataRequirements(TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
}

void TPZPrimalPoisson::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
}

void TPZPrimalPoisson::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

void TPZPrimalPoisson::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

int TPZPrimalPoisson::ClassId() const{
    return Hash("TPZPrimalPoisson") ^ TPZMaterial::ClassId() << 1;
}

void TPZPrimalPoisson::Write(TPZStream &buf, int withclassid) const{
    DebugStop();
}


void TPZPrimalPoisson::Read(TPZStream &buf, void *context){
    DebugStop();
}


/** @} */

/**
 * @name Contribute methods (weak formulation) no multiphysics mesh
 * @{
 */

/** @brief Volumetric contribute with jacobian matrix */
void TPZPrimalPoisson::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){

    TPZFNMatrix<220,REAL> &phi     = data.phi;
    TPZFNMatrix<660,REAL> &dphix   = data.dphix;
    
    TPZFNMatrix<15,STATE> &dpdx    = data.dsol[0];
    
    int nphi_p = phi.Rows();

    TPZManVector<STATE,1> f(1,0.0);
    TPZFMatrix<STATE> df;
    if (this->HasForcingFunction()) {
        this->fForcingFunction->Execute(data.x, f, df);
    }
   
    int dim = this->Dimension();
    
    for (int ip = 0; ip < nphi_p; ip++) {
        
        STATE dp_dot_dphi_i = 0.0;
        for (int i = 0; i < dim; i++) {
            dp_dot_dphi_i += dpdx[i]*dphix(i,ip);
        }
        
        ef(ip,0) += weight * (dp_dot_dphi_i - f[0] * phi(ip,0));
        
        for (int jp = 0; jp < nphi_p; jp++) {
            
            STATE dphi_j_dot_dphi_i = 0.0;
            for (int i = 0; i < this->Dimension(); i++) {
                dphi_j_dot_dphi_i += dphix(i,jp)*dphix(i,ip);
            }
            
            ek(ip,jp) += weight * dphi_j_dot_dphi_i;
        }
        
    }
    
}


/** @brief Volumetric contribute */
void TPZPrimalPoisson::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef){

    TPZFNMatrix<220,REAL> &phi     = data.phi;
    TPZFNMatrix<660,REAL> &dphix   = data.dphix;
    
    TPZFNMatrix<15,STATE> &dpdx    = data.dsol[0];
    
    int nphi_p = phi.Rows();
    
    TPZManVector<STATE,1> f(1,0.0);
    TPZFMatrix<STATE> df;
    if (this->HasForcingFunction()) {
        this->fForcingFunction->Execute(data.x, f, df);
    }
    
    for (int ip = 0; ip < nphi_p; ip++) {
        
        STATE dp_dot_dphi_i = 0.0;
        for (int i = 0; i < this->Dimension(); i++) {
            dp_dot_dphi_i += dpdx[i]*dphix(i,ip);
        }
        
        ef(ip,0) += weight * (dp_dot_dphi_i - f[0] * phi(ip,0));
        
    }
    
}


/** @brief Boundary contribute with jacobian matrix */
void TPZPrimalPoisson::ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    TPZFMatrix<REAL> &phi =     data.phi;
    TPZVec<STATE>    &p   =     data.sol[0];
    
    int nphi_p = phi.Rows();
    
    TPZManVector<STATE,1> bc_data(1,0.0);
    bc_data[0] = bc.Val2()(0,0);
    if (bc.HasForcingFunction()) {
        bc.ForcingFunction()->Execute(data.x, bc_data);
    }
    
    
    switch (bc.Type()) {
        case 0 : {      // Dirichlet condition
            STATE p_D = bc_data[0];
            for(int ip = 0 ; ip < nphi_p; ip++) {
                ef(ip,0) += weight * gBigNumber * ( p[0] -  p_D ) * phi(ip,0);
            }
        }
            break;
            
        case 1 : {      // Neumann condition
            STATE q_N = bc_data[0];
            for(int ip = 0 ; ip < nphi_p; ip++) {
                ef(ip,0) += weight * q_N * phi(ip,0);
            }
        }
            break;
        default :{
            PZError << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " - Error! boundary condition not implemented\n";
            DebugStop();
        }
            break;
    }
}


/** @brief Boundary contribute */
void TPZPrimalPoisson::ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    TPZFMatrix<REAL> &phi =     data.phi;
    TPZVec<STATE>    &p   =     data.sol[0];
    
    int nphi_p = phi.Rows();

    TPZManVector<STATE,1> bc_data(1,0.0);
    bc_data[0] = bc.Val2()(0,0);
    if (bc.HasForcingFunction()) {
        //TPZFMatrix<STATE> df;
        bc.ForcingFunction()->Execute(data.x, bc_data);   ///Jorge  2017 It is not used: , df);
    }

    
    switch (bc.Type()) {
        case 0 : {      // Dirichlet condition
            STATE p_D = bc_data[0];
            for(int ip = 0 ; ip < nphi_p; ip++) {
                ef(ip,0) += weight * gBigNumber * ( p[0] -  p_D ) * phi(ip,0);
                for (int jp = 0 ; jp < nphi_p; jp++) {
                    ek(ip,jp) += gBigNumber * phi(ip,0) * phi(jp,0) * weight;
                }
            }
        }
            break;
            
        case 1 : {      // Neumann condition
            STATE q_N = bc_data[0];
            for(int ip = 0 ; ip < nphi_p; ip++) {
                ef(ip,0) += weight * q_N * phi(ip,0);
            }
        }
            break;
        default :{
            PZError << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " - Error! boundary condition not implemented\n";
            DebugStop();
        }
            break;
    }
    
}


/** @} */


/**
 * @name Contribute methods (weak formulation) multiphysics mesh
 * @{
 */

/** @brief Volumetric contribute with jacobian matrix */
void TPZPrimalPoisson::Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    DebugStop();
}

/** @brief Volumetric contribute */
void TPZPrimalPoisson::Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef){
    DebugStop();
}


/** @brief Boundary contribute with jacobian matrix */
void TPZPrimalPoisson::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    DebugStop();
}


/** @brief Boundary contribute */
void TPZPrimalPoisson::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    DebugStop();
}


/** @} */

/**
 * @name Post-processing methods
 * @{
 */


int TPZPrimalPoisson::NEvalErrors() {return 3;}

int TPZPrimalPoisson::VariableIndex(const std::string &name){
    
    if(!strcmp("q",name.c_str()))               return  1;
    if(!strcmp("p",name.c_str()))               return  2;
    if(!strcmp("q_exact",name.c_str()))         return  3;
    if(!strcmp("p_exact",name.c_str()))         return  4;
    if(!strcmp("f_exact",name.c_str()))         return  5;
    
    return TPZMaterial::VariableIndex(name);
}

int TPZPrimalPoisson::NSolutionVariables(int var){
    if(var == 1) return this->Dimension();
    if(var == 2) return 1;
    if(var == 3) return this->Dimension();
    if(var == 4) return 1;
    if(var == 5) return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}

void TPZPrimalPoisson::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    TPZVec<STATE> p, q, f;
    
    if(var == 1){
        for (int i=0; i < this->Dimension(); i++)
        {
            Solout[i] = -data.dsol[0][i];
        }
        return;
    }
    
    if(var == 2){
            Solout[0] = data.sol[0][0];
        return;
    }
    
    if(var == 3){
        TPZManVector<STATE,1> f(1,0.0);
        TPZFNMatrix<4,STATE> df(4,1,0.0);
        if (this->HasForcingFunctionExact()) {
            this->fForcingFunctionExact->Execute(data.x, f, df);
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
        if (this->HasForcingFunctionExact()) {
            this->fForcingFunctionExact->Execute(data.x, f, df);
        }
        Solout[0] = f[0];
        return;
    }
    
    if(var == 5){
        TPZManVector<STATE,1> f(1,0.0);
        TPZFNMatrix<4,STATE> df(4,1,0.0);
        if (this->HasForcingFunctionExact()) {
            this->fForcingFunctionExact->Execute(data.x, f, df);
        }
        Solout[0] = df(3,0);
        return;
    }
    
    DebugStop();
    
}

void TPZPrimalPoisson::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    DebugStop();
}

void TPZPrimalPoisson::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &du, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &error){
    
    error.Fill(0.0);
    //  q = - grad (p)
    du *= -1.0;
    
    /** @brief   error[0] : primal error using L2 norm */
    STATE p_error = u[0] - u_exact[0];
    error[0]  = p_error*p_error;

    /** @brief   error[1] : dual error using L2 norm */
    for(int i = 0; i < this->Dimension(); i++) {
        STATE d_error = du(i,0) - du_exact(i,0);
        error[1]  += d_error*d_error;
    }
    
    /** @brief   error[2] : H1 error norm */
    error[2]= error[1];
    
}

/** @} */
