//
//  TPZLinearElliptic.cpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#include "TPZLinearElliptic.h"
#include <iostream>
#include <string>
#include "pzbndcond.h"
#include "pzaxestools.h"

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.permeabilityc"));
#endif


TPZLinearElliptic::TPZLinearElliptic():TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin>(), fnu(0.), falpha(0.), fk(0.), feta(0.), fPlaneStress(0) {
    
    fDim = 2;
    fb.resize(2);
    fb[0]=0.;
    fb[1]=0.;
    fPlaneStress = 1.;
    
    frho_s = 2700.0; // @omar:: put a method for the right set up of these values
    frho_f = 1000.0;
    
    feta_dp = 0.0;
    fxi_dp = 0.0;
    
}

TPZLinearElliptic::TPZLinearElliptic(int matid, int dim):TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin>(matid), fnu(0.), falpha(0.), fk(0.), feta(0.),fPlaneStress(0) {
    
    fDim = dim;
    fb.resize(2);
    fb[0]=0.;
    fb[1]=0.;
    fPlaneStress = 1;
    
    frho_s = 2700.0; // @omar:: put a method for the right set up of these values
    frho_f = 1000.0;
    fk_model = 0;
    feta_dp = 0.0;
    fxi_dp = 0.0;
    
}

TPZLinearElliptic::~TPZLinearElliptic(){
}


int TPZLinearElliptic::NStateVariables() {
    return 1;
}


REAL TPZLinearElliptic::k_permeability(REAL &phi, REAL &k){
    
    
    k = 0.0;
    REAL tom2 = 9.869233e-16;
    switch (fk_model) {
        case 0:
        {
            k = fk;
        }
            break;
            
        case 1:
        {
            k = fk*pow((phi/fporosity_0),4.0);
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
REAL TPZLinearElliptic::porosoty_corrected(TPZVec<TPZMaterialData> &datavec){
    
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
    REAL phi = fporosity_0 + falpha * div_u + fSe * p[0];
    
    return phi;
    
}

void TPZLinearElliptic::Compute_Sigma(TPZFMatrix<REAL> & S_eff,TPZFMatrix<REAL> & Grad_u, REAL p_ex){
    
    TPZFNMatrix<6,REAL> Grad_ut(2,2,0.0), epsilon(2,2,0.0), I(2,2,0.0);
    Grad_u.Transpose(&Grad_ut);
    
    epsilon = Grad_u + Grad_ut;
    epsilon *= 0.5;
    
    I(0,0) = 1.0;
    I(1,1) = 1.0;
    
    REAL trace = (epsilon(0,0) + epsilon(1,1));
    
    S_eff = 2.0 * fmu * epsilon + flambda * trace * I - 0.0*falpha * p_ex * I;
    
}

void TPZLinearElliptic::Compute_Sigma(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_v){
    
    TPZFNMatrix<6,REAL> Grad_vt(3,3,0.0), epsilon(3,3,0.0), I(3,3,0.0);
    Grad_v.Transpose(&Grad_vt);
    
    epsilon = Grad_v + Grad_vt;
    epsilon *= 0.5;
    
    I.Identity();
    
    REAL trace = (epsilon(0,0) + epsilon(1,1));
    
    S = 2.0 * fmu * epsilon + flambda * trace * I;
    
}

REAL TPZLinearElliptic::Inner_Product(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & T){
    REAL inner_product = S(0,0) * T(0,0) + S(0,1) * T(0,1) + S(1,0) * T(1,0) + S(1,1) * T(1,1); //     S11 T11 + S12 T12 + S21 T21 + S22 T22
    return inner_product;
}

void TPZLinearElliptic::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef){
    
    
    int u_b = 0;
    int p_b = 1;

    TPZMaterialData::MShapeFunctionType shapetype = datavec[u_b].fShapeType;
    if(shapetype == datavec[u_b].EVecShape){
        this->ContributeVec(datavec, weight, ek, ef);
    }

    
    // Getting the space functions
    TPZFMatrix<REAL>    &phiu   =   datavec[u_b].phi;
    TPZFMatrix<REAL>    &phip   =   datavec[p_b].phi;
    
    TPZFMatrix<REAL>    &dphiu   =   datavec[u_b].dphix;
    TPZFMatrix<REAL>    &dphip   =   datavec[p_b].dphix;
    
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    TPZFNMatrix<6,REAL> Grad_p(2,1,0.0),Grad_phi_i(2,1,0.0),Grad_phi_j(2,1,0.0);
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0);
    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1);
    
    int nphi_u = phiu.Rows();
    int nphi_p = phip.Rows();
    
    int first_u = 0;
    int first_p = 2*nphi_u;
    
    
    // Compute porosity poroelastic correction
    REAL phi_poro = porosoty_corrected(datavec);
    
    REAL dt = fSimulationData->dt();
    if (!fSimulationData->IsCurrentStateQ()) {
        
        return;
    }
    

    REAL rho_avg = (1.0-phi_poro)*frho_s+phi_poro*frho_f;
    fb[0] = rho_avg*fSimulationData->Gravity()[0];
    fb[1] = rho_avg*fSimulationData->Gravity()[1];
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_u(3,3,0.0),Grad_u_n,e_e,e_p,S;
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    
    Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    
    // Get the pressure at the integrations points
    int64_t global_point_index = datavec[0].intGlobPtIndex;
    TPZPoroPermMemory &point_memory = GetMemory()[global_point_index];
    e_e = point_memory.epsilon_e_n();
    e_p = point_memory.epsilon_p_n();
    Grad_u_n = point_memory.grad_u_n();
    
    Compute_Sigma(S, Grad_u);
    
    
    TPZFNMatrix<6,REAL> Grad_vx_i(2,1,0.0),Si_x;
    TPZFNMatrix<6,REAL> Grad_vy_i(2,1,0.0),Si_y;
    
    TPZFNMatrix<6,REAL> Grad_v(2,2,0.0),T(2,2,0.0);
    TPZFNMatrix<6,REAL> Grad_vx_j(2,1,0.0),Tj_x;
    TPZFNMatrix<6,REAL> Grad_vy_j(2,1,0.0),Tj_y;
    
    
    for (int iu = 0; iu < nphi_u; iu++) {
        
        // Computing Gradient of the test function for each component
        Grad_vx_i(0,0) = dphiu(0,iu)*axes_u(0,0)+dphiu(1,iu)*axes_u(1,0); // dvx/dx
        Grad_vx_i(1,0) = dphiu(0,iu)*axes_u(0,1)+dphiu(1,iu)*axes_u(1,1); // dvx/dy
        
        Grad_vy_i(0,0) = dphiu(0,iu)*axes_u(0,0)+dphiu(1,iu)*axes_u(1,0); // dvy/dx
        Grad_vy_i(1,0) = dphiu(0,iu)*axes_u(0,1)+dphiu(1,iu)*axes_u(1,1); // dvy/dy
        
        ef(2*iu + first_u, 0)   += weight * ((S(0,0) + falpha * p[0]) * Grad_vx_i(0,0) + S(0,1) * Grad_vx_i(1,0) - fb[0] * phiu(iu, 0));
        ef(2*iu+1 + first_u, 0)	+= weight * (S(1,0) * Grad_vy_i(0,0) + (S(1,1) + falpha * p[0]) * Grad_vy_i(1,0) - fb[1] * phiu(iu, 0));
        
        
        for (int ju = 0; ju < nphi_u; ju++) {
            
            
            // Computing Gradient of the test function
            Grad_vx_j(0,0) = dphiu(0,ju)*axes_u(0,0)+dphiu(1,ju)*axes_u(1,0); // dvx/dx
            Grad_vx_j(1,0) = dphiu(0,ju)*axes_u(0,1)+dphiu(1,ju)*axes_u(1,1); // dvx/dy
            
            Grad_vy_j(0,0) = dphiu(0,ju)*axes_u(0,0)+dphiu(1,ju)*axes_u(1,0); // dvy/dx
            Grad_vy_j(1,0) = dphiu(0,ju)*axes_u(0,1)+dphiu(1,ju)*axes_u(1,1); // dvy/dy
            
            
            ek(2*iu + first_u, 2*ju + first_u)      += weight * ( ( (2.0*fmu + flambda) * Grad_vx_j(0,0) ) * Grad_vx_i(0,0) + fmu * Grad_vx_j(1,0) * Grad_vx_i(1,0) );
            ek(2*iu + first_u, 2*ju+1 + first_u)    += weight * ( (flambda * Grad_vy_j(1,0) ) * Grad_vx_i(0,0) + fmu * Grad_vy_j(0,0) * Grad_vx_i(1,0)  );
            ek(2*iu+1 + first_u, 2*ju + first_u)	+= weight * ( fmu * Grad_vx_j(1,0) * Grad_vy_i(0,0) + flambda * Grad_vx_j(0,0) * Grad_vy_i(1,0));
            ek(2*iu+1 + first_u, 2*ju+1 + first_u)	+= weight * ( (2.0*fmu + flambda) * Grad_vy_j(1,0) * Grad_vy_i(1,0) + fmu * Grad_vy_j(0,0) * Grad_vy_i(0,0) );
            
        }
        
    }
    

    for (int ip = 0; ip < nphi_p; ip++) {
        
        ef(ip + first_p, 0)		+= weight *  ( 0.0 * phip(ip,0));
        
        for (int jp = 0; jp < nphi_p; jp++) {
            
            ek(ip + first_p, jp + first_p)		+= weight * phip(jp,0) * phip(ip,0);
        }
        
    }
    
    
    
}

void TPZLinearElliptic::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
    
}

void TPZLinearElliptic::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int u_b = 0;
    int p_b = 1;
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[u_b].fShapeType;
    if(shapetype == datavec[u_b].EVecShape){
        this->ContributeVecBC(datavec, weight, ek, ef, bc);
        return;
    }
    
    TPZFMatrix<REAL>  &phiu = datavec[u_b].phi;
    TPZFMatrix<REAL>  &phip = datavec[p_b].phi;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    int phru = phiu.Rows();
//    int phrp = phip.Rows();
    short in,jn;
    REAL v[3];
    v[0] = bc.Val2()(0,0);	//	Ux displacement
    v[1] = bc.Val2()(1,0);	//	Uy displacement
    v[2] = bc.Val2()(2,0);	//	Pressure
    
    REAL time = this->SimulationData()->t();
    REAL Value = bc.Val2()(0,0);
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<REAL,3> f(3);
        TPZFMatrix<REAL> gradf;
        bc.TimedependentBCForcingFunction()->Execute(datavec[p_b].x, time, f, gradf);
        v[0] = f[0];	//	Ux displacement or Tx
        v[1] = f[1];	//	Uy displacement or Ty
        v[2] = f[2];	//	Pressure
    }
    else{
        Value = bc.Val2()(0,0);
    }
    
    // Dirichlet in Pressure
    switch (bc.Type())
    {
        case 0 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(2*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(2*in+1,2*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            
            break;
        }
            
        case 1 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                }
            }
            
            break;
        }
            
        case 2 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            
            break;
        }
            
        case 3 :
        {
            //	Neumann condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(2*in,0)		+= -1.0*v[0]*phiu(in,0)*weight;		//	Tnx
                ef(2*in+1,0)	+= -1.0*v[1]*phiu(in,0)*weight;		//	Tny
            }
            
            break;
        }
            
        case 4 :
        {
            //	Neumann condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(2*in,0)		+= -1.0*v[0]*phiu(in,0)*weight;		//	Tnx
            }
            
            break;
        }
            
        case 5 :
        {
            //	Neumann condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(2*in+1,0)	+= -1.0*v[1]*phiu(in,0)*weight;		//	Tny
            }
            
            break;
        }
            
        case 6 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(2*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(2*in+1,2*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            
            break;
        }
            
        case 7 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)		+= gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                }
            }
            
            break;
        }
            
        case 8 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in+1,0)	+= gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)	+= gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            

            break;
        }
            
        case 9 :
        {
            //	Neumann condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(2*in,0)		+= -1.0*v[0]*phiu(in,0)*weight;		//	Tnx
                ef(2*in+1,0)	+= -1.0*v[1]*phiu(in,0)*weight;		//	Tny
            }
            
            break;
        }
            
        case 10 :
        {
            //	Neumann condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(2*in,0)		+= -1.0*v[0]*phiu(in,0)*weight;		//	Tnx
            }

            break;
        }
            
        case 11 :
        {
            //	Neumann condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(2*in+1,0)	+= -1.0*v[1]*phiu(in,0)*weight;		//	Tny
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


void TPZLinearElliptic::ContributeVec(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){

    int u_b = 0;
    int p_b = 1;
    
    // Getting the space functions
    TPZFMatrix<REAL>    &phiu   =   datavec[u_b].phi;
    TPZFMatrix<REAL>    &phip   =   datavec[p_b].phi;
    
    TPZFMatrix<REAL>    &dphiu   =   datavec[u_b].dphix;
    TPZFMatrix<REAL>    &dphip   =   datavec[p_b].dphix;
    
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    TPZFNMatrix<6,REAL> Grad_p(2,1,0.0),Grad_phi_i(2,1,0.0),Grad_phi_j(2,1,0.0);
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0);
    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1);
    
    int nphi_u = phiu.Rows();
    int nphi_p = phip.Rows();
    
    int first_u = 0;
    int first_p = 2*nphi_u;
    
    
    // Compute porosity poroelastic correction
    REAL phi_poro = porosoty_corrected(datavec);
    
    REAL dt = fSimulationData->dt();
    if (!fSimulationData->IsCurrentStateQ()) {
        
        return;
    }
    
    REAL rho_avg = (1.0-phi_poro)*frho_s+phi_poro*frho_f;
    fb[0] = rho_avg*fSimulationData->Gravity()[0];
    fb[1] = rho_avg*fSimulationData->Gravity()[1];
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_u(3,3,0.0),Grad_u_n,e_e,e_p,S;
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    
    Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    
    // Get the pressure at the integrations points
    int64_t global_point_index = datavec[0].intGlobPtIndex;
    TPZPoroPermMemory &point_memory = GetMemory()[global_point_index];
    e_e = point_memory.epsilon_e_n();
    e_p = point_memory.epsilon_p_n();
    Grad_u_n = point_memory.grad_u_n();
    
    Compute_Sigma(S, Grad_u);
    
    
    TPZFNMatrix<6,REAL> Grad_vx_i(2,1,0.0),Si_x;
    TPZFNMatrix<6,REAL> Grad_vy_i(2,1,0.0),Si_y;
    
    TPZFNMatrix<6,REAL> Grad_v(2,2,0.0),T(2,2,0.0);
    TPZFNMatrix<6,REAL> Grad_vx_j(2,1,0.0),Tj_x;
    TPZFNMatrix<6,REAL> Grad_vy_j(2,1,0.0),Tj_y;
    
    
    for (int iu = 0; iu < nphi_u; iu++) {
        
        // Computing Gradient of the test function for each component
        Grad_vx_i(0,0) = dphiu(0,iu)*axes_u(0,0)+dphiu(1,iu)*axes_u(1,0); // dvx/dx
        Grad_vx_i(1,0) = dphiu(0,iu)*axes_u(0,1)+dphiu(1,iu)*axes_u(1,1); // dvx/dy
        
        Grad_vy_i(0,0) = dphiu(2,iu)*axes_u(0,0)+dphiu(3,iu)*axes_u(1,0); // dvy/dx
        Grad_vy_i(1,0) = dphiu(2,iu)*axes_u(0,1)+dphiu(3,iu)*axes_u(1,1); // dvy/dy

        
        ef(2*iu + first_u, 0)   += weight * ((S(0,0) + falpha * p[0]) * Grad_vx_i(0,0) + S(0,1) * Grad_vx_i(1,0) - fb[0] * phiu(iu, 0));
        ef(2*iu+1 + first_u, 0)	+= weight * (S(1,0) * Grad_vy_i(0,0) + (S(1,1) + falpha * p[0]) * Grad_vy_i(1,0) - fb[1] * phiu(iu, 0));
        
        
        for (int ju = 0; ju < nphi_u; ju++) {
            
            
            // Computing Gradient of the test function for each component
            Grad_vx_j(0,0) = dphiu(0,ju)*axes_u(0,0)+dphiu(1,ju)*axes_u(1,0); // dvx/dx
            Grad_vx_j(1,0) = dphiu(0,ju)*axes_u(0,1)+dphiu(1,ju)*axes_u(1,1); // dvx/dy
            
            Grad_vy_j(0,0) = dphiu(2,ju)*axes_u(0,0)+dphiu(3,ju)*axes_u(1,0); // dvy/dx
            Grad_vy_j(1,0) = dphiu(2,ju)*axes_u(0,1)+dphiu(3,ju)*axes_u(1,1); // dvy/dy
            
            
            ek(2*iu + first_u, 2*ju + first_u)      += weight * ( ( (2.0*fmu + flambda) * Grad_vx_j(0,0) ) * Grad_vx_i(0,0) + fmu * Grad_vx_j(1,0) * Grad_vx_i(1,0) );
            ek(2*iu + first_u, 2*ju+1 + first_u)    += weight * ( (flambda * Grad_vy_j(1,0) ) * Grad_vx_i(0,0) + fmu * Grad_vy_j(0,0) * Grad_vx_i(1,0)  );
            ek(2*iu+1 + first_u, 2*ju + first_u)	+= weight * ( fmu * Grad_vx_j(1,0) * Grad_vy_i(0,0) + flambda * Grad_vx_j(0,0) * Grad_vy_i(1,0));
            ek(2*iu+1 + first_u, 2*ju+1 + first_u)	+= weight * ( (2.0*fmu + flambda) * Grad_vy_j(1,0) * Grad_vy_i(1,0) + fmu * Grad_vy_j(0,0) * Grad_vy_i(0,0) );
            
        }
        
    }
    
    
    for (int ip = 0; ip < nphi_p; ip++) {
        
        ef(ip + first_p, 0)		+= weight *  ( 0.0 * phip(ip,0));
        
        for (int jp = 0; jp < nphi_p; jp++) {
            
            ek(ip + first_p, jp + first_p)		+= weight * phip(jp,0) * phip(ip,0);
        }
        
    }

}

void TPZLinearElliptic::ContributeVec(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    DebugStop();
}

void TPZLinearElliptic::ContributeVecBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    return;
}

void TPZLinearElliptic::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)

{
//    datavec[0].fShapeType = TPZMaterialData::EVecShape;
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

void TPZLinearElliptic::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){

//    datavec[0].fShapeType = TPZMaterialData::EVecShape;
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
        datavec[i].fNeedsNeighborSol = true;
    }
}

void TPZLinearElliptic::Print(std::ostream &out)
{
    out << "Material Name : " << Name() << "\n";
    out << "Plane Problem (fPlaneStress = 0, for Plane Strain conditions) " << fPlaneStress << std::endl;
    out << "Properties for TPZLinearElliptic: \n";
    out << "\t Poisson Ratio   = "											<< fnu		<< std::endl;
    out << "\t Undarined Poisson Ratio   = "								<< fnuu		<< std::endl;
    out << "\t First Lamé Parameter   = "									<< flambda	<< std::endl;
    out << "\t Second Lamé Parameter   = "									<< fmu		<< std::endl;
    out << "\t Undrained First Lamé Parameter   = "							<< flambdau	<< std::endl;
    out << "\t Biot coefficient   = "										<< falpha	<< std::endl;
    out << "\t Body force vector B {X-direction, Y-direction}   = "			<< fb[0] << ' ' << fb[1]   << std::endl;
    out << "Properties for Diffusion: \n";
    out << "\t Permeability   = "											<< fk		<< std::endl;
    out << "\t Fluid Viscosity   = "										<< feta	<< std::endl;
    out << "\t Constrained specific storage at constant strain Se = "		<< fSe		<< std::endl;
    out << "Class properties :";
    TPZMaterial::Print(out);
    out << "\n";
    
}

/** Returns the variable index associated with the name */
int TPZLinearElliptic::VariableIndex(const std::string &name)
{
    //	Elasticity Variables
    if(!strcmp("u",name.c_str()))				return	1;
    if(!strcmp("s_x",name.c_str()))             return	2;
    if(!strcmp("s_y",name.c_str()))             return	3;
    if(!strcmp("s_z",name.c_str()))             return	4;
    if(!strcmp("t_xy",name.c_str()))            return	5;
    
    //	Diffusion Variables
    if(!strcmp("p_ex",name.c_str()))				return	6;
    if(!strcmp("v",name.c_str()))					return	7;
    if(!strcmp("k_x",name.c_str()))					return	8;
    if(!strcmp("k_y",name.c_str()))					return	9;
    if(!strcmp("phi",name.c_str()))					return	10;
    
    if(!strcmp("e_x",name.c_str()))             return	11;
    if(!strcmp("e_y",name.c_str()))             return	12;
    if(!strcmp("e_xy",name.c_str()))            return	13;
    if(!strcmp("ep_x",name.c_str()))             return	14;
    if(!strcmp("ep_y",name.c_str()))             return	15;
    if(!strcmp("ep_xy",name.c_str()))            return	16;
    
    if(!strcmp("K_0",name.c_str()))            return	17;
    
    return TPZMaterial::VariableIndex(name);
}

int TPZLinearElliptic::NSolutionVariables(int var){
    if(var == 1)	return fDim;
    if(var == 2)	return 1;
    if(var == 3)	return 1;
    if(var == 4)	return 1;
    if(var == 5)	return 1;
    if(var == 6)	return 1;
    if(var == 7)	return fDim;
    if(var == 8)	return 1;
    if(var == 9)	return 1;
    if(var == 10)	return 1;
    if(var == 11)	return 1;
    if(var == 12)	return 1;
    if(var == 13)	return 1;
    if(var == 14)	return 1;
    if(var == 15)	return 1;
    if(var == 16)	return 1;
    if(var == 17)	return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}

//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
void TPZLinearElliptic::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    
    int u_b = 0;
    int p_b = 1;
    
    
    // Getting the space functions
    
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    
    REAL to_Mpa     = 1.0e-6;
    REAL to_Darcy   = 1.013249966e+12;
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_p(3,1,0.0),Grad_u(3,3,0.0),Grad_u_n(3,3,0.0),e_e(3,3,0.0),e_p(3,3,0.0),S;
    
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0); // dp/dx
    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1); // dp/dy
    
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    
    Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    
    Compute_Sigma(S, Grad_u);
    
    //	Displacements
    if(var == 1){
        Solout[0] = u[0];
        Solout[1] = u[1];
        return;
    }
    
    //	sigma_x
    if(var == 2) {
        Solout[0] = S(0,0)*to_Mpa;
        return;
    }
    
    //	sigma_y
    if(var == 3) {
        Solout[0] = S(1,1)*to_Mpa;
        return;
    }
    
    //	sigma_z
    if(var == 4) {
        Solout[0] = S(2,2)*to_Mpa;
        return;
    }
    
    //	tau_xy
    if(var == 5) {
        Solout[0] = S(0,1)*to_Mpa;
        return;
    }
    
    //	Pore pressure excess
    if(var == 6) {
        Solout[0] = p[0]*to_Mpa;
        return;
    }
    
    //	v
    if(var == 7) {
        REAL phi = porosoty_corrected(datavec);
        REAL k;
        k_permeability(phi, k);
        Solout[0] = -(k/feta) * Grad_p(0,0);
        Solout[1] = -(k/feta) * Grad_p(1,0);
        return;
    }
    
    //	k_x
    if(var == 8) {
        REAL phi = porosoty_corrected(datavec);
        REAL k = 0.0;
        k_permeability(phi, k);
        Solout[0] = k*to_Darcy;
        return;
    }
    
    //	k_y
    if(var == 9) {
        REAL phi = porosoty_corrected(datavec);
        REAL k = 0.0;
        k_permeability(phi, k);
        Solout[0] = k*to_Darcy;
        return;
    }
    
    //	Porosity form poroelastic correction
    if(var == 10) {
        Solout[0] = porosoty_corrected(datavec);
        return;
    }
    
    //	epsilon_x
    if(var == 11) {
        Solout[0] = e_e(0,0);
        return;
    }
    
    //	epsilon_y
    if(var == 12) {
        Solout[0] = e_e(1,1);
        return;
    }
    
    //	epsilon_xy
    if(var == 13) {
        Solout[0] = e_e(0,1);
        return;
    }
    
    //	epsilon_p_x
    if(var == 14) {
        Solout[0] = e_p(0,0);
        return;
    }
    
    //	epsilon_p_y
    if(var == 15) {
        Solout[0] = e_p(1,1);
        return;
    }
    
    //	epsilon_p_xy
    if(var == 16) {
        Solout[0] = e_p(0,1);
        return;
    }
    
    //	K_0
    if(var == 17) {
        Solout[0] = S(0,0)/S(1,1);
        return;
    }
    
}

