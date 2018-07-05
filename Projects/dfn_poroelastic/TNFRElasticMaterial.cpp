//
//  TNFRElasticMaterial.cpp
//  pz
//
//  Created by Omar Durán on 7/2/18.
//

#include "TNFRElasticMaterial.h"


/// Constructor
TNFRElasticMaterial::TNFRElasticMaterial(){
    
}

/// Constructor based on material identifier
TNFRElasticMaterial::TNFRElasticMaterial(int id) : TPZMatWithMem<TNRFElasticMemory,TPZMaterial>(id) {
    
}

/// Destructor
TNFRElasticMaterial::~TNFRElasticMaterial(){
    
}

/// material dimension
int TNFRElasticMaterial::Dimension() const {
    return 2;
}

int TNFRElasticMaterial::NStateVariables(){
    return Dimension();
}

/// setting data on domain
void TNFRElasticMaterial::FillDataRequirements(TPZMaterialData &data){
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
}

///  setting data on boundary
void TNFRElasticMaterial::FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data){
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
    
}

/// jacobian contribution
void TNFRElasticMaterial::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    unsigned int m_dim =  Dimension();
    TPZFMatrix<REAL> & phi        = data.phi;
    TPZFNMatrix<20,REAL> grad_phi;
    TPZAxesTools<STATE>::Axes2XYZ(data.dphix, grad_phi, data.axes);
    
    TPZFNMatrix<3,REAL> grad_u;
    TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], grad_u, data.axes);
    int n_phi = phi.Rows();
    
    TPZFNMatrix<9,REAL> Grad_vx_i(m_dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_i(m_dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vx_j(m_dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_j(m_dim,1,0.0);
    
    TPZFNMatrix<9,REAL> S(3,3);
    grad_u.Resize(3, 3);
    
    //  Computing sigma
    REAL trace;
    for (int i = 0; i < 3; i++) {
        trace = 0.0;
        for (int j = 0; j < 3; j++) {
            S(i,j) = m_mu * (grad_u(i,j) + grad_u(j,i));
            trace +=  grad_u(j,j);
        }
        S(i,i) += m_lambda * trace;
    }
    
    REAL dvxdx, dvxdy;
    REAL dvydx, dvydy;
    
    REAL duxdx, duxdy;
    REAL duydx, duydy;
    
    for (unsigned int i = 0; i < n_phi; i++) {
        
        for (int d = 0; d < m_dim; d++) {
            Grad_vx_i(d,0) = grad_phi(d,i);
            Grad_vy_i(d,0) = grad_phi(d,i);
        }
        
        dvxdx = Grad_vx_i(0,0);
        dvxdy = Grad_vx_i(1,0);
        
        dvydx = Grad_vy_i(0,0);
        dvydy = Grad_vy_i(1,0);
        
        ef(2*i, 0)   += weight * (S(0,0) * Grad_vx_i(0,0) + S(0,1) * Grad_vx_i(1,0));
        ef(2*i+1, 0) += weight * (S(1,0) * Grad_vy_i(0,0) + S(1,1) * Grad_vy_i(1,0));
        
        for (int j = 0; j < n_phi; j++) {
            
            for (int d = 0; d < m_dim; d++) {
                Grad_vx_j(d,0) = grad_phi(d,j);
                Grad_vy_j(d,0) = grad_phi(d,j);
            }
            
            duxdx = Grad_vx_j(0,0);
            duxdy = Grad_vx_j(1,0);
            
            duydx = Grad_vy_j(0,0);
            duydy = Grad_vy_j(1,0);
            
            ek(2*i, 2*j)      += weight * ( (2.0*m_mu + m_lambda)*duxdx*dvxdx + m_mu*duxdy*dvxdy);
            ek(2*i, 2*j+1)    += weight * ( m_lambda*duydy*dvxdx + m_mu*duydx*dvxdy);
            ek(2*i+1, 2*j)    += weight * ( m_lambda*duxdx*dvydy + m_mu*duxdy*dvydx);
            ek(2*i+1, 2*j+1)    += weight * ( (2.0*m_mu + m_lambda)*duydy*dvydy + m_mu*duydx*dvydx);
            
        }
    }
}

/// residual contribution
void TNFRElasticMaterial::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    Contribute(data, weight, ek_fake, ef);
}

void TNFRElasticMaterial::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    // Get the data at the integrations points
    TPZFMatrix<REAL>  &phi = data.phi;
    TPZManVector<REAL,10>  &u = data.sol[0];
    TPZManVector<REAL,3>  &n = data.normal;
    
    int n_phi = phi.Rows();
    short i,j;

    switch (bc.Type())
    {
        case 0 : // Displacement case
        {
            TPZManVector<REAL,2> v(2);
            v[0] = bc.Val2()(0,0);
            v[1] = bc.Val2()(1,0);
            
            REAL dux = u[0] - v[0];
            REAL duy = u[1] - v[1];
            for(i = 0 ; i < n_phi; i++)
            {
                ef(2*i,0)      += weight * gBigNumber * dux * phi(i,0); // x direction
                ef(2*i+1,0)    += weight * gBigNumber * duy * phi(i,0); // y direction
                
                for (j = 0 ; j < n_phi; j++)
                {
                    ek(2*i,2*j)        += weight * gBigNumber * phi(i,0) * phi(j,0);
                    ek(2*i+1,2*j+1)    += weight * gBigNumber * phi(i,0) * phi(j,0);
                }
            }
            break;
        }
            
        case 1 : // Normal displacement case
        {
            
            REAL u_dot_n    = bc.Val2()(0,0);
            REAL uh_dot_n   = u[0]*n[0] + u[1]*n[1];
            
            REAL dun = uh_dot_n - u_dot_n;
            for(i = 0 ; i < n_phi; i++)
            {
                ef(2*i,0)      += weight * gBigNumber * dun * phi(i,0) * n[0]; // x direction
                ef(2*i+1,0)    += weight * gBigNumber * dun * phi(i,0) * n[1]; // y direction
                
                for (j = 0 ; j < n_phi; j++)
                {
                    ek(2*i,2*j)        += weight * gBigNumber * phi(i,0) * phi(j,0) * n[0];
                    ek(2*i+1,2*j+1)    += weight * gBigNumber * phi(i,0) * phi(j,0) * n[1];
                }
            }
            break;
        }
            
        case 2 : // Normal stress case
        {
            TPZTensor<REAL> S;
            S.Zero();
            S.XX() = bc.Val2()(0,0);
            S.XY() = bc.Val2()(1,0);
            S.XZ() = bc.Val2()(2,0);
            S.YY() = bc.Val2()(3,0);
            S.YZ() = bc.Val2()(4,0);
            S.ZZ() = bc.Val2()(5,0);
            
            TPZManVector<REAL,3> S_dot_n(3);
            
            S_dot_n[0] = S.XX()*n[0] + S.XY()*n[1] + S.XZ()*n[2];
            S_dot_n[1] = S.XY()*n[0] + S.YY()*n[1] + S.YZ()*n[2];
            S_dot_n[2] = S.XZ()*n[0] + S.YZ()*n[1] + S.ZZ()*n[2];
            
            for(i = 0 ; i < n_phi; i++)
            {
                ef(2*i,0)      += -1.0 * weight * S_dot_n[0] * phi(i,0); //    S_dot_n_x
                ef(2*i+1,0)    += -1.0 * weight * S_dot_n[1] * phi(i,0); //    S_dot_n_y
            }
            break;
        }
        case 3 : // Displacement x case
        {
            TPZManVector<REAL,2> v(1);
            v[0] = bc.Val2()(0,0);
            
            REAL dux = u[0] - v[0];
            for(i = 0 ; i < n_phi; i++)
            {
                ef(2*i,0)      += weight * gBigNumber * dux * phi(i,0); // x direction
                
                for (j = 0 ; j < n_phi; j++)
                {
                    ek(2*i,2*j)        += weight * gBigNumber * phi(i,0) * phi(j,0);
                }
            }
            break;
        }
        case 4 : // Displacement y case
        {
            TPZManVector<REAL,2> v(1);
            v[0] = bc.Val2()(0,0);
            
            REAL duy = u[1] - v[0];
            for(i = 0 ; i < n_phi; i++)
            {
                ef(2*i+1,0)    += weight * gBigNumber * duy * phi(i,0); // y direction
                
                for (j = 0 ; j < n_phi; j++)
                {
                    ek(2*i+1,2*j+1)    += weight * gBigNumber * phi(i,0) * phi(j,0);
                }
            }
            break;
        }
        case 5 : // Null Normal displacement case
        {
            
            for(i = 0 ; i < n_phi; i++)
            {
                for (j = 0 ; j < n_phi; j++)
                {
                    ek(2*i,2*j)        += weight * gBigNumber * phi(i,0) * phi(j,0) * n[0];
                    ek(2*i+1,2*j+1)    += weight * gBigNumber * phi(i,0) * phi(j,0) * n[1];
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

/// Returns the variable index associated with the name
int TNFRElasticMaterial::VariableIndex(const std::string &name){
    
    if (!strcmp("u", name.c_str()))   return 0;
    if (!strcmp("s_xx", name.c_str())) return 1;
    if (!strcmp("s_xy", name.c_str())) return 2;
    if (!strcmp("s_xz", name.c_str())) return 3;
    if (!strcmp("s_yy", name.c_str())) return 4;
    if (!strcmp("s_yz", name.c_str())) return 5;
    if (!strcmp("s_zz", name.c_str())) return 6;
    if (!strcmp("s_1", name.c_str())) return 7;
    if (!strcmp("s_2", name.c_str())) return 8;
    if (!strcmp("s_3", name.c_str())) return 9;
    if (!strcmp("i_1", name.c_str())) return 10;
    if (!strcmp("j_2", name.c_str())) return 11;
    
    return TPZMatWithMem::VariableIndex(name);
}

/// Returns the number of variables associated with provided index var
int TNFRElasticMaterial::NSolutionVariables(int var){
    switch(var) {
        case 0:
            return Dimension(); // vector
        case 1:
            return 1; // scalar
        case 2:
            return 1; // scalar
        case 3:
            return 1; // scalar
        case 4:
            return 1; // scalar
        case 5:
            return 1; // scalar
        case 6:
            return 1; // scalar
        case 7:
            return 1; // scalar
        case 8:
            return 1; // scalar
        case 9:
            return 1; // scalar
        case 10:
            return 1; // scalar
        case 11:
            return 1; // scalar
    }
    
    DebugStop();
    return -1;
}

/// post-processing
void TNFRElasticMaterial::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    
    TPZManVector<REAL,3> u = data.sol[0];
    TPZFNMatrix<3,REAL> grad_u;
    TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], grad_u, data.axes);
    
    TPZFNMatrix<9,REAL> S(3,3);
    grad_u.Resize(3, 3);
    
    //  Computing sigma
    REAL trace;
    for (int i = 0; i < 3; i++) {
        trace = 0.0;
        for (int j = 0; j < 3; j++) {
            S(i,j) = m_mu * (grad_u(i,j) + grad_u(j,i));
            trace +=  grad_u(j,j);
        }
        S(i,i) += m_lambda * trace;
    }
    
    TPZTensor<REAL> T;
    T.XX() = S(0,0);
    T.XY() = S(0,1);
    T.XZ() = S(0,2);
    T.YY() = S(1,1);
    T.YZ() = S(1,2);
    T.ZZ() = S(2,2);
    
    switch (var) {
        case 0:
        {
            for (unsigned int d = 0; d < Dimension(); d++) {
                Solout[d] = u[d];
            }
        }
            break;
        case 1:
        {
            Solout[0] = T.XX();
        }
            break;
        case 2:
        {
            Solout[0] = T.XY();
        }
            break;
        case 3:
        {
            Solout[0] = T.XZ();
        }
            break;
        case 4:
        {
            Solout[0] = T.YY();
        }
            break;
        case 5:
        {
            Solout[0] = T.YZ();
        }
            break;
        case 6:
        {
            Solout[0] = T.ZZ();
        }
            break;
        case 7:
        {
            TPZTensor<REAL>::TPZDecomposed sig_eigen_system;
            T.EigenSystem(sig_eigen_system);
            Solout[0] = sig_eigen_system.fEigenvalues[0];
        }
            break;
        case 8:
        {
            TPZTensor<REAL>::TPZDecomposed sig_eigen_system;
            T.EigenSystem(sig_eigen_system);
            Solout[0] = sig_eigen_system.fEigenvalues[1];
        }
            break;
        case 9:
        {
            TPZTensor<REAL>::TPZDecomposed sig_eigen_system;
            T.EigenSystem(sig_eigen_system);
            Solout[0] = sig_eigen_system.fEigenvalues[2];
        }
            break;
        case 10:
        {
            Solout[0] = T.I1();
        }
            break;
        case 11:
        {
            Solout[0] = T.I2();
        }
            break;
        default:
            std::cout  << " not implemented. " << std::endl;
            DebugStop();
            break;
    }
}

void TNFRElasticMaterial::SetLameParameters(REAL lambda, REAL mu){
    m_lambda = lambda;
    m_mu     = mu;
}

