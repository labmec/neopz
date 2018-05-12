//
//  TPZPoroPermCoupling3D.cpp
//  PZ
//
//  Created by Manouchehr on 4/16/18.
//
//

#include "TPZPoroPermCoupling3D.h"
#include <iostream>
#include <string>
#include "pzbndcond.h"
#include "pzaxestools.h"
#include <algorithm>

#include "pzfmatrix.h"
#include "TPZTensor.h"


#include "TPZSandlerExtended.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZElasticCriterion.h"


#include "TPZSandlerDimaggio.h"

/** @brief default constructor */
TPZPoroPermCoupling3D::TPZPoroPermCoupling3D():TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin>()
{

m_Dim = 3;
    
}

/** @brief constructor based on a material id */
TPZPoroPermCoupling3D::TPZPoroPermCoupling3D(int matid, int dim):TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin>(matid) {

    m_Dim = dim;
    
}

/** @brief default desconstructor */
TPZPoroPermCoupling3D::~TPZPoroPermCoupling3D(){
}


/** @brief copy constructor $ */
TPZPoroPermCoupling3D::TPZPoroPermCoupling3D(const TPZPoroPermCoupling3D& other){
    this->m_Dim    = other.m_Dim;
    this->m_SimulationData    = other.m_SimulationData;
}


/** @brief copy assignemnt operator $ */
TPZPoroPermCoupling3D& TPZPoroPermCoupling3D::operator = (const TPZPoroPermCoupling3D& other){
    
    if (this != & other) // prevent self-assignment
    {
        this->m_Dim    = other.m_Dim;
        this->m_SimulationData    = other.m_SimulationData;
    }
    return *this;
}


 /** @brief Set the required data at each integration point */
void TPZPoroPermCoupling3D::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)

{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNeighborSol = true;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = true;
    }
}


/** @brief Set the required data at each integration point; BoundaryConditionData */
void TPZPoroPermCoupling3D::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
        datavec[i].fNeedsNeighborSol = true;
    }
}


/** returns the number of state variables associated with the material */
int TPZPoroPermCoupling3D::NStateVariables() {
    return m_Dim;
}


/** @brief permeability model */
REAL TPZPoroPermCoupling3D::k_permeability(REAL &phi, REAL &k){
    
    
    k = 0.0;
    REAL tom2 = 9.869233e-16;
    switch (m_k_model) {
        case 0:
        {
            k = m_k;
        }
            break;
            
        case 1:
        {
            k = m_k*pow((phi/m_porosity_0),4.0);
        }
            break;
            
        case 2:
        {
            k = 0.136*(pow(phi,1.4))*tom2;
        }
            break;
            
        case 3:
        {
            k = (100.0*pow(phi,2.25))*(100.0*pow(phi,2.25))*tom2;
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
    
    return k;
}

/** @brief Poroelastic porosity correction */
REAL TPZPoroPermCoupling3D::porosoty_corrected(TPZVec<TPZMaterialData> &datavec){
    
    int u_b = 0;
    int p_b = 1;
    
    // Getting the space functions
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    
    TPZFNMatrix<6,REAL> Grad_u(2,2,0.0);
    
    // Computing Gradient of the Solution
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(1,0) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    
    Grad_u(0,1) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    
    REAL div_u = Grad_u(0,0) + Grad_u(1,1);
    REAL phi = m_porosity_0 + m_alpha * div_u + m_Se * p[0];
    
    return phi;
    
}



/** @brief the print function */
void TPZPoroPermCoupling3D::Print(std::ostream &out) {
    out << "\t Base on the class print:\n";
    out << " Name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}



/** Returns the variable index associated with the name */
int TPZPoroPermCoupling3D::VariableIndex(const std::string &name)
{
    //	Elasticity Variables
    if(!strcmp("u",name.c_str()))				return	1;
    if(!strcmp("s_x",name.c_str()))             return	2;
    if(!strcmp("s_y",name.c_str()))             return	3;
    if(!strcmp("s_z",name.c_str()))             return	4;
    if(!strcmp("t_xy",name.c_str()))            return	5;
    if(!strcmp("t_xz",name.c_str()))            return	6;
    if(!strcmp("t_yz",name.c_str()))            return	7;
    
    //	Diffusion Variables
    if(!strcmp("p",name.c_str()))				return	8;
    if(!strcmp("v",name.c_str()))					return	9;
    if(!strcmp("k_x",name.c_str()))					return	10;
    if(!strcmp("k_y",name.c_str()))					return	11;
    if(!strcmp("k_z",name.c_str()))					return	12;
    if(!strcmp("phi",name.c_str()))					return	13;
    
    //	Defformation Variables

    if(!strcmp("e_x",name.c_str()))             return	14;
    if(!strcmp("e_y",name.c_str()))             return	15;
    if(!strcmp("e_z",name.c_str()))             return	16;
    if(!strcmp("e_xy",name.c_str()))            return	17;
    if(!strcmp("e_xz",name.c_str()))            return	18;
    if(!strcmp("e_yz",name.c_str()))            return	19;
    
    if(!strcmp("ep_x",name.c_str()))             return	20;
    if(!strcmp("ep_y",name.c_str()))             return	21;
    if(!strcmp("ep_z",name.c_str()))             return	22;
    if(!strcmp("ep_xy",name.c_str()))            return	23;
    if(!strcmp("ep_xz",name.c_str()))            return	24;
    if(!strcmp("ep_yz",name.c_str()))            return	25;
    
    
    if(!strcmp("K_0",name.c_str()))            return	26;
    
    return TPZMaterial::VariableIndex(name);
}



/** returns the number of variables associated with the variable indexed  */
int TPZPoroPermCoupling3D::NSolutionVariables(int var){
    if(var == 1)	return m_Dim;
    if(var == 2)	return 1;
    if(var == 3)	return 1;
    if(var == 4)	return 1;
    if(var == 5)	return 1;
    if(var == 6)	return 1;
    if(var == 7)	return 1;
    if(var == 8)	return 1;
    if(var == 9)	return m_Dim;
    if(var == 10)	return 1;
    if(var == 11)	return 1;
    if(var == 12)	return 1;
    if(var == 13)	return 1;
    if(var == 14)	return 1;
    if(var == 15)	return 1;
    if(var == 16)	return 1;
    if(var == 17)	return 1;
    if(var == 18)	return 1;
    if(var == 19)	return 1;
    if(var == 20)	return 1;
    if(var == 21)	return 1;
    if(var == 22)	return 1;
    if(var == 23)	return 1;
    if(var == 24)	return 1;
    if(var == 25)	return 1;
    if(var == 26)	return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}


/** @brief Compute effective stress */
void TPZPoroPermCoupling3D::Compute_Sigma(REAL & l, REAL & mu, TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_u){
    
    
    REAL trace;
    for (int i = 0; i < 3; i++) {
        trace = 0.0;
        for (int j = 0; j < 3; j++) {
            S(i,j) = mu * (Grad_u(i,j) + Grad_u(j,i));
            trace +=  Grad_u(j,j);
        }
        S(i,i) += l * trace;
    }
    
    return;
}



// @brief the inner product function
REAL TPZPoroPermCoupling3D::Inner_Product(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & T){
    REAL inner_product = S(0,0) * T(0,0) + S(0,1) * T(0,1) + S(0,2) * T(0,2) + S(1,0) * T(1,0) + S(1,1) * T(1,1) + S(1,2) * T(1,2) + S(2,0) * T(2,0) + S(2,1) * T(2,1) + S(2,2) * T(2,2); //     S11 T11 + S12 T12 + S13 T13 + S21 T21 + S22 T22 + S23 T23 + S31 T31 + S32 T32 + S33 T33
    return inner_product;
}


// Contribute Methods being used

void TPZPoroPermCoupling3D::Contribute_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int u_b = 0;
    int p_b = 1;
    
    // Compute porosity poroelastic correction
    REAL phi_poro = porosoty_corrected(datavec);
    
    REAL dt = m_SimulationData->dt();
    
    if (!m_SimulationData->IsCurrentStateQ()) {
        return;
    }
    
    
    REAL rho_avg = (1.0-phi_poro)*m_rho_s+phi_poro*m_rho_f;
    m_b[0] = rho_avg*m_SimulationData->Gravity()[0];
    m_b[1] = rho_avg*m_SimulationData->Gravity()[1];
    m_b[2] = rho_avg*m_SimulationData->Gravity()[2];
    
    
    // Get the solution at the integrations points
    long global_point_index = datavec[0].intGlobPtIndex;
    TPZPoroPermMemory &point_memory = GetMemory()[global_point_index];
    
    
    // Getting the space functions
    TPZFMatrix<REAL>    & phi_u   =   datavec[u_b].phi;
    TPZFMatrix<REAL>    & phi_p   =   datavec[p_b].phi;
    
    TPZFMatrix<REAL>    & grad_phi_u   =   datavec[u_b].dphix;
    TPZFMatrix<REAL>    & grad_phi_p   =   datavec[p_b].dphix;

    
    TPZFMatrix<REAL> & grad_u_0 = point_memory.grad_u_0();
    TPZFMatrix<REAL> & grad_u   = point_memory.grad_u();
    TPZFMatrix<REAL> & grad_u_n = point_memory.grad_u_n();
    
    TPZFMatrix<REAL> Grad_phi_u,Grad_u;
    TPZAxesTools<STATE>::Axes2XYZ(datavec[u_b].dphix, Grad_phi_u, datavec[u_b].axes);
    TPZAxesTools<STATE>::Axes2XYZ(datavec[u_b].dsol[0], Grad_u, datavec[u_b].axes);
    
    TPZFMatrix<REAL> Grad_phi_p,Grad_p;
    TPZAxesTools<STATE>::Axes2XYZ(datavec[p_b].dphix, Grad_phi_p, datavec[p_b].axes);
    TPZAxesTools<STATE>::Axes2XYZ(datavec[p_b].dsol[0], Grad_p, datavec[p_b].axes);
    
    
    
    
    REAL & p_n  = point_memory.p_n();
    REAL & p_0  = point_memory.p_0();
    
    
    
    int nphi_u = phi_u.Rows();
    int nphi_p = phi_p.Rows();
    
    int first_u = 0;
    int first_p = 2*nphi_u;
    
    TPZFNMatrix<9,REAL> Grad_vx_i(m_Dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_i(m_Dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vz_i(m_Dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vx_j(m_Dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_j(m_Dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vz_j(m_Dim,1,0.0);
    
    
    TPZFNMatrix<9,REAL> Grad_wx_i(m_Dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_wy_i(m_Dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_wz_i(m_Dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_wx_j(m_Dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_wy_j(m_Dim,1,0.0);
    TPZFNMatrix<9,REAL> Grad_wz_j(m_Dim,1,0.0);
    
    
    
    TPZFNMatrix<9,REAL> S_0(3,3),S(3,3),S_n(3,3);
    Compute_Sigma(m_lambda, m_mu, S_0, grad_u_0);
    Compute_Sigma(m_lambda, m_mu, S, grad_u);
    Compute_Sigma(m_lambda, m_mu, S_n, grad_u_n);
    
    if (m_SimulationData->IsInitialStateQ()) {
        S_0.Zero();
    }
    else{
        S_n += S_0; // S_0 is the total initial stress
    }
    
    
    REAL source = m_alpha * (p_n - p_0);
    
    REAL dvxdx, dvxdy, dvxdz;
    REAL dvydx, dvydy, dvydz;
    REAL dvzdx, dvzdy, dvzdz;
    
    REAL duxdx, duxdy, duxdz;
    REAL duydx, duydy, duydz;
    REAL duzdx, duzdy, duzdz;
    
    
    
    for (int iu = 0; iu < nphi_u; iu++)
    {
        
        // Computing Gradient of the test function for each component
        for (int d = 0; d < m_Dim; d++)
        {
            Grad_vx_i(d,0) = grad_phi_u(d,iu);
            Grad_vy_i(d,0) = grad_phi_u(d,iu);
            Grad_vz_i(d,0) = grad_phi_u(d,iu);
        }
    
        
        ef(3*iu + first_u, 0)   += weight * ((S_n(0,0) - source) * Grad_vx_i(0,0) + S_n(0,1) * Grad_vx_i(1,0) + S_n(0,2) * Grad_vx_i(2,0) - m_b[0] * phi_u(iu,0));
        ef(3*iu+1 + first_u, 0)	+= weight * (S_n(1,0) * Grad_vy_i(0,0) + (S_n(1,1) - source) * Grad_vy_i(1,0) + S_n(1,2) * Grad_vy_i(2,0) - m_b[1] * phi_u(iu,0));
        ef(3*iu+2 + first_u, 0)	+= weight * (S_n(2,0) * Grad_vz_i(0,0) + S_n(2,1) * Grad_vz_i(1,0) + (S_n(2,2) - source) * Grad_vz_i(2,0) - m_b[2] * phi_u(iu,0));
        
        //x
        dvxdx = Grad_vx_i(0,0);
        dvxdy = Grad_vx_i(1,0);
        dvxdz = Grad_vx_i(2,0);
        
        //y
        dvydx = Grad_vy_i(0,0);
        dvydy = Grad_vy_i(1,0);
        dvydz = Grad_vy_i(2,0);
        
        //z
        dvzdx = Grad_vz_i(0,0);
        dvzdy = Grad_vz_i(1,0);
        dvzdz = Grad_vz_i(2,0);
    
        
        for (int ju = 0; ju < nphi_u; ju++)
        {
            
            // Computing Gradient of the test function for each component
            for (int d = 0; d < m_Dim; d++)
            {
                Grad_vx_j(d,0) = grad_phi_u(d,ju);
                Grad_vy_j(d,0) = grad_phi_u(d,ju);
                Grad_vz_j(d,0) = grad_phi_u(d,ju);
            }
            
            //x
            duxdx = Grad_vx_j(0,0);
            duxdy = Grad_vx_j(1,0);
            duxdz = Grad_vx_j(2,0);
            
            //y
            duydx = Grad_vy_j(0,0);
            duydy = Grad_vy_j(1,0);
            duydz = Grad_vy_j(2,0);
            
            //z
            duzdx = Grad_vz_j(0,0);
            duzdy = Grad_vz_j(1,0);
            duzdz = Grad_vz_j(2,0);
        
        
            
            // Gradient 1
            ek(3*iu + first_u, 3*ju + first_u)      += weight * ((m_lambda + 2.*m_mu)*duxdx*dvxdx + m_mu*duxdy*dvxdy + m_mu*duxdz*dvxdz);
            ek(3*iu + first_u, 3*ju+1  + first_u)   += weight * (m_lambda*duydy*dvxdx + m_mu*duydx*dvxdy);
            ek(3*iu + first_u, 3*ju+2  + first_u)   += weight * (m_lambda*duzdz*dvxdx + m_mu*duzdx*dvxdz);
            
            // Gradient 2
            ek(3*iu+1 + first_u, 3*ju  + first_u)   += weight * (m_lambda*duxdx*dvydy + m_mu*duxdy*dvydx);
            ek(3*iu+1 + first_u, 3*ju+1  + first_u) += weight * ((m_lambda + 2.*m_mu)*duydy*dvydy + m_mu*duydx*dvydx + m_mu*duydz*dvydz);
            ek(3*iu+1 + first_u, 3*ju+2  + first_u) += weight * (m_lambda*duzdz*dvydy + m_mu*duzdy*dvydz);
            
            // Gradient 3
            ek(3*iu+2 + first_u, 3*ju  + first_u)   += weight * (m_lambda*duxdx*dvzdz + m_mu*duxdz*dvzdx);
            ek(3*iu+2 + first_u, 3*ju+1  + first_u) += weight * (m_lambda*duydy*dvzdz + m_mu*duydz*dvzdy);
            ek(3*iu+2 + first_u, 3*ju+2  + first_u) += weight * ((m_lambda + 2.*m_mu)*duzdz*dvzdz + m_mu*duzdx*dvzdx + m_mu*duzdy*dvzdy);
            
        }

    }
    
    
    TPZFNMatrix<6,REAL> dv(3,1,0.0);
    
    //	Matrix -Qc
    //	Coupling matrix
    for(int iu = 0; iu < nphi_u; iu++ )
    {
        
        // Computing Gradient of the test function for each component
        for (int d = 0; d < m_Dim; d++)
        {
            Grad_vx_i(d,0) = grad_phi_u(d,iu); // dvx/dx
            Grad_vy_i(d,0) = grad_phi_u(d,iu); // dvy/dy
            Grad_vz_i(d,0) = grad_phi_u(d,iu); // dvz/dz
        }
        
        for(int jp = 0; jp < nphi_p; jp++)
        {
            
            ek(2*iu,first_p+jp) += (-1.)* weight * m_alpha * phi_p(jp,0) * Grad_vx_i(0,0);
            ek(2*iu+1,first_p+jp) += (-1.)* weight * m_alpha * phi_p(jp,0) * Grad_vy_i(1,0);
            ek(2*iu+2,first_p+jp) += (-1.)* weight * m_alpha * phi_p(jp,0) * Grad_vz_i(2,0);
        }
    }
    
    //	Matrix QcˆT
    //	Coupling matrix transpose
    for(int ip = 0; ip < nphi_p; ip++ )
    {
        
        
        for(int ju = 0; ju < nphi_u; ju++)
        {
            
            dv(0,0) = Grad_vx_i(0,0);
            dv(1,0) = Grad_vy_i(1,0);
            dv(2,0) = Grad_vz_i(2,0);
            
            ek(first_p+ip,2*ju) += (-1.) * weight * m_alpha * dv(0,0) * phi_p(ip,0);
            ek(first_p+ip,2*ju+1) += (-1.) * weight * m_alpha * dv(1,0) * phi_p(ip,0);
            ek(first_p+ip,2*ju+2) += (-1.) * weight * m_alpha * dv(2,0) * phi_p(ip,0);
            
        }
    }
    

    
    /** @brief Rudnicki diffusion coefficient */
    /** J. W. Rudnicki. Fluid mass sources and point forces in linear elastic di usive solids. Journal of Mechanics of Materials, 5:383–393, 1986. */
    REAL k = 0.0;
    k_permeability(phi_poro,k);
    m_lambdau *=1.1;
    REAL c = 1.0;//(k/feta)*(flambdau-flambda)*(flambda + 2.0*fmu)/(falpha*falpha*(flambdau + 2.0*fmu));

    // Darcy mono-phascis flow
    for (int ip = 0; ip < nphi_p; ip++)
    {
        
        // Computing Gradient of the test function for each component
        for (int d = 0; d < m_Dim; d++)
        {
            Grad_wx_i(d,0) = grad_phi_p(d,ip); // dvx/dx
            Grad_wx_i(d,0) = grad_phi_p(d,ip); // dvy/dy
            Grad_wx_i(d,0) = grad_phi_p(d,ip); // dvz/dz
        }
        
        
        
        REAL dot = 0.0;
        for (int i = 0;  i < m_Dim; i++) {
            dot += Grad_p(i,0) * Grad_wx_i(i,0);
        }
        
        ef(ip + first_p, 0)		+= -1.0 * weight * (dt * c * dot + (phi_poro) * phi_p(ip,0));
        
        for (int jp = 0; jp < nphi_p; jp++) {
            
            Grad_phi_p(0,0) = grad_phi_p(0,jp);
            Grad_phi_p(1,0) = grad_phi_p(1,jp);
            Grad_phi_p(2,0) = grad_phi_p(2,jp);
            
            REAL dot = 0.0;
            for (int i = 0;  i < m_Dim; i++) {
                dot += Grad_wx_i(i,0) * Grad_wx_i(i,0);
            }
            
            ek(ip + first_p, jp + first_p)		+= -1.0 * weight * (dt * c * dot + m_Se * phi_p(jp,0) * phi_p(ip,0) );
        }
        
    }
    
    
}


void TPZPoroPermCoupling3D::Contribute_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute_3D(datavec, weight, ek_fake, ef);
    
}


void TPZPoroPermCoupling3D::ContributeBC_3D(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    if (!m_SimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int u_b = 0;
    int p_b = 1;
    
    REAL dt = m_SimulationData->dt();
    
    // Get the data at the integrations points
    TPZMatWithMem<TPZPoroPermMemory,TPZBndCond>  & material_mem = dynamic_cast<TPZMatWithMem<TPZPoroPermMemory,TPZBndCond > & >(bc);
    long global_point_index = datavec[u_b].intGlobPtIndex;
    TPZPoroPermMemory & point_memory = material_mem.GetMemory()[global_point_index];
    
    TPZFMatrix<REAL>  & phiu = datavec[u_b].phi;
    TPZFMatrix<REAL>  &phip = datavec[p_b].phi;
    
    TPZFNMatrix<3,REAL> u_m = point_memory.u_n();
    TPZManVector<REAL,1> p_m = point_memory.p_n();
    
    TPZManVector<REAL,3> u(3,0.0);
    u[0] = u_m(0,0);
    u[1] = u_m(0,1);
    u[2] = u_m(0,2);
    
    int phru = phiu.Rows();
    int phrp = phip.Rows();
    short in,jn;
    
    
    switch (bc.Type())
    {
        case 0 : // Du_Dp
            // Dirichlet in displacement and Pressure
        {
            REAL v[4];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);	  //	Uz displacement
            v[3] = bc.Val2()(3,0);    //    Pressure
            
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[2] - v[2])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in,3*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            
            //    Diffusion Equation
            REAL p_s = p_m[0];
            REAL d_p = (p_s-v[3]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }
            
        case 1 : // Duxy_Dp
            // Dirichlet in x and y direction of displacement and Pressure
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in,3*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            //    Diffusion Equation
            REAL p_s = p_m[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }
            
        case 2 : // Duxz_Dp
            // Dirichlet in x and z direction of displacement and Pressure
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uz displacement
            v[2] = bc.Val2()(2,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in,3*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// z displacement
                }
            }
            //    Diffusion Equation
            REAL p_s = p_m[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }

        case 3 : // Duyz_Dp
            // Dirichlet in y and z direction of displacement and Pressure
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            v[1] = bc.Val2()(1,0);    //    Uz displacement
            v[2] = bc.Val2()(2,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+1,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// y displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// y displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// z displacement
                }
            }
            //    Diffusion Equation
            REAL p_s = p_m[0];
            REAL d_p = (p_s-v[2]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }

        case 4 : // Dux_Dp
            // Dirichlet in x direction of displacement and Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
 
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in,3*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
            
                }
            }
            //    Diffusion Equation
            REAL p_s = p_m[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }
       
            
        case 5 : // Duy_Dp
            // Dirichlet in y direction of displacement and Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            v[1] = bc.Val2()(1,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+1,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// y displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// y displacement
                    
                }
            }
            //    Diffusion Equation
            REAL p_s = p_m[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }

            
        case 6 : // Duz_Dp
            // Dirichlet in z direction of displacement and Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Uz displacement
            v[1] = bc.Val2()(1,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+2,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// z displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+2,3*jn+2)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// z displacement
                    
                }
            }
            //    Diffusion Equation
            REAL p_s = p_m[0];
            REAL d_p = (p_s-v[1]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
            
        }
            
         
        case 7 :
            // Neumann condition in displacement and Dirichlet in Pressure
        {
            REAL v[4];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);	  //	Uz displacement
            v[3] = bc.Val2()(3,0);    //    Pressure
            
 
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(3*in,0)		+= -1.0 * weight * v[0] * phiu(in,0);		//	Tnx
                ef(3*in+1,0)	+= -1.0 * weight * v[1] * phiu(in,0);		//	Tny
                ef(3*in+2,0)	+= -1.0 * weight * v[2] * phiu(in,0);		//	Tnz
            }
            
            //    Diffusion Equation
            REAL p_s = p_m[0];
            REAL d_p = (p_s-v[3]);
            for(in = 0 ; in < phrp; in++)
            {
                //    Contribution for load Vector
                ef(in+3*phru,0)        += gBigNumber*(d_p)*phip(in,0)*weight;    // P Pressure
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(in+3*phru,jn+3*phru)        += gBigNumber*phip(in,0)*phip(jn,0)*weight;    // P Pressure
                }
            }
            break;
        }
            
        case 8 : // Du_Dq
            // Dirichlet in displacement and Neumann in Pressure
        {
            REAL v[4];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);	  //	Uz displacement
            v[3] = bc.Val2()(3,0);    //    Pressure
            
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[2] - v[2])*phiu(in,0)*weight;	// Z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in,3*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Z displacement
                }
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += -1.0 * dt * v[3]*phip(in,0)*weight;    // Qnormal
            }
            break;
            
        }
            
        case 9 : // Duxy_Dq
            // Dirichlet in x and y direction of displacement and Neumann in Pressure
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in,3*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+1,3*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += -1.0 * dt * v[2]*phip(in,0)*weight;    // Qnormal
            }
            break;
            
        }
            
        case 10 : // Duxz_Dq
            // Dirichlet in x and z direction of displacement and Neumann in Pressure
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uz displacement
            v[2] = bc.Val2()(2,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in,3*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// z displacement
                }
            }
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += -1.0 * dt * v[2]*phip(in,0)*weight;    // Qnormal
            }
            break;
            
        }
            
        case 11 : // Duyz_Dq
            // Dirichlet in y and z direction of displacement and Neumann in Pressure
        {
            REAL v[3];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            v[1] = bc.Val2()(1,0);    //    Uz displacement
            v[2] = bc.Val2()(2,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+1,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// y displacement Value
                ef(3*in+2,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// z displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// y displacement
                    ek(3*in+2,3*jn+2)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// z displacement
                }
            }
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += -1.0 * dt * v[2]*phip(in,0)*weight;    // Qnormal
            }
            break;
            
        }
            
        case 12 : // Dux_Dq
            // Dirichlet in x direction of displacement and Neumann in Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in,3*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    
                }
            }
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += -1.0 * dt * v[1]*phip(in,0)*weight;    // Qnormal
            }
            break;
            
        }
            
            
        case 13 : // Duy_Dq
            // Dirichlet in y direction of displacement and Neumann in Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Uy displacement
            v[1] = bc.Val2()(1,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+1,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// y displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+1,3*jn+1)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// y displacement
                    
                }
            }
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += -1.0 * dt * v[1]*phip(in,0)*weight;    // Qnormal
            }
            break;
            
        }
            
            
        case 14 : // Duz_Dq
            // Dirichlet in z direction of displacement and Neumann in Pressure
        {
            REAL v[2];
            v[0] = bc.Val2()(0,0);    //    Uz displacement
            v[1] = bc.Val2()(1,0);    //    Pressure
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(3*in+2,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// z displacement Value
                
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(3*in+2,3*jn+2)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// z displacement
                    
                }
            }
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += -1.0 * dt * v[1]*phip(in,0)*weight;    // Qnormal
            }
            break;
            
        }
            
            
        case 15 : // Nt_Nq
            // Neumann condition in displacement and Pressure
        {
            REAL v[4];
            v[0] = bc.Val2()(0,0);    //    Ux displacement
            v[1] = bc.Val2()(1,0);    //    Uy displacement
            v[2] = bc.Val2()(2,0);	  //	Uz displacement
            v[3] = bc.Val2()(3,0);    //    Pressure
            
            
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(3*in,0)		+= -1.0 * weight * v[0] * phiu(in,0);		//	Tnx
                ef(3*in+1,0)	+= -1.0 * weight * v[1] * phiu(in,0);		//	Tny
                ef(3*in+2,0)	+= -1.0 * weight * v[2] * phiu(in,0);		//	Tnz
            }
            
            //    Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //    Normal Flux on neumman boundary
                ef(in+3*phru,0)    += -1.0 * dt * v[3]*phip(in,0)*weight;    // Qnormal
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






//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
void TPZPoroPermCoupling3D::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    
    int u_b = 0;
    int p_b = 1;
    
    
    // Getting the space functions
    
    TPZFNMatrix <9,REAL> axes_u	=	datavec[u_b].axes;
    TPZFNMatrix <9,REAL> axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,3> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    
    REAL to_Mpa     = 1.0e-6;
    REAL to_Darcy   = 1.013249966e+12;
    
    // Computing Gradient of the Solution
    TPZFNMatrix<9,REAL> Grad_p(3,1,0.0),Grad_u(3,3,0.0),Grad_u_n(3,3,0.0),e_e(3,3,0.0),e_p(3,3,0.0),S(3,3,0.0);
    
    
    TPZAxesTools<STATE>::Axes2XYZ(dp, Grad_p, axes_p);
    TPZAxesTools<STATE>::Axes2XYZ(du, Grad_u, axes_u);
    
  
    Grad_u.Resize(3, 3);
    this->Compute_Sigma(m_lambda, m_mu, S, Grad_u);
    
    
    
    
    
    //	Displacements
    if(var == 1){
        Solout[0] = u[0];
        Solout[1] = u[1];
        Solout[2] = u[2];
        return;
    }

    //	sigma_x
    if(var == 2) {
        Solout[0] = S(0,0)*to_Mpa - m_alpha * p[0]*to_Mpa;
        return;
    }
    
    //	sigma_y
    if(var == 3) {
        Solout[0] = S(1,1)*to_Mpa - m_alpha * p[0]*to_Mpa;
        return;
    }
    
    //	sigma_z
    if(var == 4) {
        Solout[0] = S(2,2)*to_Mpa - m_alpha * p[0]*to_Mpa;
        return;
    }
    
    //	tau_xy
    if(var == 5) {
        Solout[0] = S(0,1)*to_Mpa;
        return;
    }
    
    //	tau_xz
    if(var == 6) {
        Solout[0] = S(0,2)*to_Mpa;
        return;
    }
    
    //	tau_yz
    if(var == 7) {
        Solout[0] = S(1,2)*to_Mpa;
        return;
    }
    
    //	Pore pressure excess
    if(var == 8) {
        Solout[0] = p[0]*to_Mpa;
        return;
    }
    
    //	v
    if(var == 9) {
        REAL phi = porosoty_corrected(datavec);
        REAL k;
        k_permeability(phi, k);
        Solout[0] = -(k/m_eta) * Grad_p(0,0);
        Solout[1] = -(k/m_eta) * Grad_p(1,0);
        Solout[2] = -(k/m_eta) * Grad_p(2,0);
        return;
    }
    
    //	k_x
    if(var == 10) {
        REAL phi = porosoty_corrected(datavec);
        REAL k = 0.0;
        k_permeability(phi, k);
        Solout[0] = k*to_Darcy;
        return;
    }
    
    //	k_y
    if(var == 11) {
        REAL phi = porosoty_corrected(datavec);
        REAL k = 0.0;
        k_permeability(phi, k);
        Solout[0] = k*to_Darcy;
        return;
    }
    
    //	k_z
    if(var == 12) {
        REAL phi = porosoty_corrected(datavec);
        REAL k = 0.0;
        k_permeability(phi, k);
        Solout[0] = k*to_Darcy;
        return;
    }
    
    //	Porosity form poroelastic correction
    if(var == 13) {
        Solout[0] = porosoty_corrected(datavec);
        return;
    }
    
    //	epsilon_x
    if(var == 14) {
        Solout[0] = e_e(0,0);
        return;
    }
    
    //	epsilon_y
    if(var == 15) {
        Solout[0] = e_e(1,1);
        return;
    }
    
    //	epsilon_z
    if(var == 16) {
        Solout[0] = e_e(2,2);
        return;
    }
    
    //	epsilon_xy
    if(var == 17) {
        Solout[0] = e_e(0,1);
        return;
    }
    
    //	epsilon_xz
    if(var == 18) {
        Solout[0] = e_e(0,2);
        return;
    }
    
    //	epsilon_yz
    if(var == 19) {
        Solout[0] = e_e(1,2);
        return;
    }
    
    //	epsilon_p_x
    if(var == 20) {
        Solout[0] = e_p(0,0);
        return;
    }
    
    //	epsilon_p_y
    if(var == 21) {
        Solout[0] = e_p(1,1);
        return;
    }
    
    //	epsilon_p_z
    if(var == 21) {
        Solout[0] = e_p(2,2);
        return;
    }
    
    //	epsilon_p_xy
    if(var == 22) {
        Solout[0] = e_p(0,1);
        return;
    }
    
    //	epsilon_p_xz
    if(var == 23) {
        Solout[0] = e_p(0,2);
        return;
    }
    
    //	epsilon_p_yz
    if(var == 24) {
        Solout[0] = e_p(1,2);
        return;
    }
    
    //	K_0
    if(var == 25) {
        Solout[0] = S(0,0)/S(1,1);
        return;
    }
    
}


