//
//  TPZPoroPermCoupling.cpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#include "TPZPoroPermCoupling.h"
#include <iostream>
#include <string>
#include "pzbndcond.h"
#include "pzaxestools.h"


#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.permeabilityc"));
#endif


TPZPoroPermCoupling::TPZPoroPermCoupling():TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin>(), fnu(0.), falpha(0.), fk(0.), feta(0.), fPlaneStress(0) {

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

TPZPoroPermCoupling::TPZPoroPermCoupling(int matid, int dim):TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin>(matid), fnu(0.), falpha(0.), fk(0.), feta(0.),fPlaneStress(0) {

    fDim = dim;
    fb.resize(2);
    fb[0]=0.;
    fb[1]=0.;
    fPlaneStress = 1;
    
    frho_s = 2700.0; // @omar:: put a method for the right set up of these values
    frho_f = 1000.0;
    
    feta_dp = 0.0;
    fxi_dp = 0.0;
    
}

TPZPoroPermCoupling::~TPZPoroPermCoupling(){
}


int TPZPoroPermCoupling::NStateVariables() {
    return 1;
}


REAL TPZPoroPermCoupling::c_diffusion(REAL phi){

    REAL c = 0.0;
    REAL k = fk*(phi/fporosity_0);
    c = (k/feta)*(flambdau-flambda)*(flambdau + 2.0*fmu)/(falpha*falpha*(flambdau + 2.0*fmu));
    
    return c;
}

/** @brief Poroelastic porosity correction */
REAL TPZPoroPermCoupling::porosoty_corrected(TPZVec<TPZMaterialData> &datavec){
    
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

void TPZPoroPermCoupling::Compute_Sigma(TPZFMatrix<REAL> & S_eff,TPZFMatrix<REAL> & Grad_u, REAL p_ex){
    
    TPZFNMatrix<6,REAL> Grad_ut(2,2,0.0), epsilon(2,2,0.0), I(2,2,0.0);
    Grad_u.Transpose(&Grad_ut);
    
    epsilon = Grad_u + Grad_ut;
    epsilon *= 0.5;
    
    I(0,0) = 1.0;
    I(1,1) = 1.0;
    
    REAL trace = (epsilon(0,0) + epsilon(1,1));
    
    S_eff = 2.0 * fmu * epsilon + flambda * trace * I - 0.0*falpha * p_ex * I;
    
}

void TPZPoroPermCoupling::Compute_Sigma(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_v){
    
    TPZFNMatrix<6,REAL> Grad_vt(2,2,0.0), epsilon(2,2,0.0), I(2,2,0.0);
    Grad_v.Transpose(&Grad_vt);
    
    epsilon = Grad_v + Grad_vt;
    epsilon *= 0.5;
    
    I(0,0) = 1.0;
    I(1,1) = 1.0;
    
    REAL trace = (epsilon(0,0) + epsilon(1,1));
    
    S = 2.0 * fmu * epsilon + flambda * trace * I;
    
}

REAL TPZPoroPermCoupling::Inner_Product(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & T){
    REAL inner_product = S(0,0) * T(0,0) + S(0,1) * T(0,1) + S(1,0) * T(1,0) + S(1,1) * T(1,1); //     S11 T11 + S12 T12 + S21 T21 + S22 T22
    return inner_product;
}

void TPZPoroPermCoupling::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef){
    
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
        

        // Darcy mono-phascis flow
        for (int ip = 0; ip < nphi_p; ip++) {
            
            ef(ip + first_p, 0)		+= weight *  (-1.0/dt) * (phi_poro)  * phip(ip,0);
        }
        
        return;
    }
    
    REAL rho_avg = (1.0-phi_poro)*frho_s+phi_poro*frho_f;
    fb[0] = rho_avg*fSimulationData->Gravity()[0];
    fb[1] = rho_avg*fSimulationData->Gravity()[1];

    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_u(2,2,0.0),S(2,2,0.0);
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    
    Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
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
        
        ef(2*iu + first_u, 0)   += weight * (S(0,0) * Grad_vx_i(0,0) + S(0,1) * Grad_vx_i(1,0) - (-1.0*falpha * Grad_p(0,0) + fb[0])*phiu(iu, 0));
        ef(2*iu+1 + first_u, 0)	+= weight * (S(1,0) * Grad_vy_i(0,0) + S(1,1) * Grad_vy_i(1,0) - (-1.0*falpha * Grad_p(1,0) + fb[1])*phiu(iu, 0));
        
        
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
    
    TPZFNMatrix<6,REAL> dv(2,1,0.0);
    
    //	Matrix Qc
    //	Coupling matrix
    for(int iu = 0; iu < nphi_u; iu++ )
    {
//        dv(0,0) = dphiu(0,iu)*axes_u(0,0)+dphiu(1,iu)*axes_u(1,0);
//        dv(1,0) = dphiu(0,iu)*axes_u(0,1)+dphiu(1,iu)*axes_u(1,1);
        
        for(int jp = 0; jp < nphi_p; jp++)
        {
            Grad_phi_j(0,0) = dphip(0,jp)*axes_p(0,0)+dphip(1,jp)*axes_p(1,0);
            Grad_phi_j(1,0) = dphip(0,jp)*axes_p(0,1)+dphip(1,jp)*axes_p(1,1);
            
            ek(2*iu,first_p+jp) += (+1.)* weight * 1.0*falpha * Grad_phi_j(0,0) * phiu(iu,0);
            ek(2*iu+1,first_p+jp) += (+1.)* weight * 1.0*falpha *Grad_phi_j(1,0) * phiu(iu,0);
        }
    }
    
    //	Matrix QcˆT
    //	Coupling matrix transpose
    for(int ip = 0; ip < nphi_p; ip++ )
    {
        
        
        for(int ju = 0; ju < nphi_u; ju++)
        {
            
            dv(0,0) = dphiu(0,ju)*axes_u(0,0)+dphiu(1,ju)*axes_u(1,0);
            dv(1,0) = dphiu(0,ju)*axes_u(0,1)+dphiu(1,ju)*axes_u(1,1);
            
            ek(first_p+ip,2*ju) += (1.) * weight * (1.0/dt) * falpha * dv(0,0) * phip(ip,0);
            ek(first_p+ip,2*ju+1) += (1.) * weight * (1.0/dt) * falpha * dv(1,0) * phip(ip,0);
            
        }
    }
    
    REAL c = c_diffusion(phi_poro);
    
    // Darcy mono-phascis flow
    for (int ip = 0; ip < nphi_p; ip++) {
        
        Grad_phi_i(0,0) = dphip(0,ip)*axes_p(0,0)+dphip(1,ip)*axes_p(1,0);
        Grad_phi_i(1,0) = dphip(0,ip)*axes_p(0,1)+dphip(1,ip)*axes_p(1,1);
        
        REAL dot = 0.0;
        for (int i = 0;  i < fDim; i++) {
            dot += Grad_p(i,0) * Grad_phi_i(i,0);
        }
        
        ef(ip + first_p, 0)		+= weight *  (c * dot + (1.0/dt) * (phi_poro) * phip(ip,0));
        
        for (int jp = 0; jp < nphi_p; jp++) {
            
            Grad_phi_j(0,0) = dphip(0,jp)*axes_p(0,0)+dphip(1,jp)*axes_p(1,0);
            Grad_phi_j(1,0) = dphip(0,jp)*axes_p(0,1)+dphip(1,jp)*axes_p(1,1);
            
            REAL dot = 0.0;
            for (int i = 0;  i < fDim; i++) {
                dot += Grad_phi_j(i,0) * Grad_phi_i(i,0);
            }
            
            ek(ip + first_p, jp + first_p)		+= weight * ( c * dot + (1.0/dt) * (fSe * phip(jp,0)) * phip(ip,0) );
        }
        
    }
    
    
//    // Darcy mono-phascis flow
//    for (int ip = 0; ip < nphi_p; ip++) {
//        
//       
//        ef(ip + first_p, 0)		+= weight * (p[0] - M_PI)  * phip(ip,0);
//        
//        for (int jp = 0; jp < nphi_p; jp++) {
//            
//            ek(ip + first_p, jp + first_p)		+= weight * ( phip(jp,0) * phip(ip,0) );
//        }
//        
//    }
    
}

void TPZPoroPermCoupling::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
    
}

void TPZPoroPermCoupling::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int u_b = 0;
    int p_b = 1;
    
    TPZFMatrix<REAL>  &phiu = datavec[u_b].phi;
    TPZFMatrix<REAL>  &phip = datavec[p_b].phi;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    int phru = phiu.Rows();
    int phrp = phip.Rows();
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
        DebugStop();

//        bc.TimedependentBCForcingFunction()->Execute(datavec[p_b].x, time, f, gradf);
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
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in+2*phru,0)		+= gBigNumber*(p[0]-v[2])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
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
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in+2*phru,0)		+= gBigNumber*(p[0]-v[2])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
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
            
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in+2*phru,0)		+= gBigNumber*(p[0]-v[2])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
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
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in+2*phru,0)		+= gBigNumber*(p[0]-v[2])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
                }
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
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in+2*phru,0)		+= gBigNumber*(p[0]-v[2])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
                }
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
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in+2*phru,0)		+= gBigNumber*(p[0]-v[2])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in+2*phru,jn+2*phru)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
                }
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
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0*v[2]*phip(in,0)*weight;	// Qnormal
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
        
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0*v[2]*phip(in,0)*weight;	// Qnormal
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

            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0*v[2]*phip(in,0)*weight;	// Qnormal
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
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0*v[2]*phip(in,0)*weight;	// Qnormal
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
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0*v[2]*phip(in,0)*weight;	// Qnormal
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
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0*v[2]*phip(in,0)*weight;	// Qnormal
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

void TPZPoroPermCoupling::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)

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

void TPZPoroPermCoupling::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
        datavec[i].fNeedsNeighborSol = true;
    }
}

void TPZPoroPermCoupling::Print(std::ostream &out)
{
    out << "Material Name : " << Name() << "\n";
    out << "Plane Problem (fPlaneStress = 0, for Plane Strain conditions) " << fPlaneStress << std::endl;
    out << "Properties for TPZPoroPermCoupling: \n";
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
int TPZPoroPermCoupling::VariableIndex(const std::string &name)
{
    //	Elasticity Variables
    if(!strcmp("u",name.c_str()))				return	1;
    if(!strcmp("s_x",name.c_str()))             return	2;
    if(!strcmp("s_y",name.c_str()))             return	3;
    if(!strcmp("t_xy",name.c_str()))            return	4;
    
    //	Diffusion Variables
    if(!strcmp("p_ex",name.c_str()))				return	5;
    if(!strcmp("v",name.c_str()))					return	6;
    if(!strcmp("k_x",name.c_str()))					return	7;
    if(!strcmp("k_y",name.c_str()))					return	8;
    if(!strcmp("phi",name.c_str()))					return	9;
    
    if(!strcmp("e_x",name.c_str()))             return	10;
    if(!strcmp("e_y",name.c_str()))             return	11;
    if(!strcmp("e_xy",name.c_str()))            return	12;
    
    return TPZMaterial::VariableIndex(name);
}

int TPZPoroPermCoupling::NSolutionVariables(int var){
    if(var == 1)	return fDim;
    if(var == 2)	return 1;
    if(var == 3)	return 1;
    if(var == 4)	return 1;
    if(var == 5)	return 1;
    if(var == 6)	return fDim;
    if(var == 7)	return 1;
    if(var == 8)	return 1;
    if(var == 9)	return 1;
    if(var == 10)	return 1;
    if(var == 11)	return 1;
    if(var == 12)	return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}

//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
void TPZPoroPermCoupling::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    
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
    
    
    int nphi_u = phiu.Rows();
    int nphi_p = phip.Rows();
    
    int first_u = 0;
    int first_p = 2*nphi_u;
    
    TPZFNMatrix<6,REAL> Grad_u(2,2,0.0),Grad_u_t,S,e;
    REAL to_Mpa     = 1.0e-6;
    REAL to_Darcy   = 1.013249966e+12;
    // Computing Gradient of the Solution
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(1,0) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    Grad_u(0,1) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    Grad_u.Transpose(&Grad_u_t);
    e = 0.5 * (Grad_u + Grad_u_t);
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
    
    //	tau_xy
    if(var == 4) {
        Solout[0] = S(0,1)*to_Mpa;
        return;
    }
    
    //	Pore pressure excess
    if(var == 5) {
        Solout[0] = p[0]*to_Mpa;
        return;
    }
    
    //	k_x
    if(var == 7) {
        REAL phi = porosoty_corrected(datavec);
        Solout[0] = fk*(phi/fporosity_0)*to_Darcy;
        return;
    }
    
    //	k_y
    if(var == 6) {
        REAL phi = porosoty_corrected(datavec);
        Solout[0] = fk*(phi/fporosity_0)*to_Darcy;
        return;
    }
    
    //	Porosity form poroelastic correction
    if(var == 9) {
        Solout[0] = porosoty_corrected(datavec);
        return;
    }
    
    //	epsilon_x
    if(var == 10) {
        Solout[0] = e(0,0);
        return;
    }
    
    //	epsilon_y
    if(var == 11) {
        Solout[0] = e(1,1);
        return;
    }
    
    //	epsilon_xy
    if(var == 12) {
        Solout[0] = e(0,1);
        return;
    }
    
    
    //	Darcy's velocity
//    if (var == 7)
//    { 
//        int id;
//        TPZManVector<STATE> dsolp(2,0);
//        dsolp[0] = datavec[1].dsol[0](0,0)*datavec[1].axes(0,0)+datavec[1].dsol[0](1,0)*datavec[1].axes(1,0);
//        dsolp[1] = datavec[1].dsol[0](0,0)*datavec[1].axes(0,1)+datavec[1].dsol[0](1,0)*datavec[1].axes(1,1);			
//        for(id=0 ; id<fDim; id++) 
//        {
//            Solout[id] = -1. * this->fK * dsolp[id];
//        }
//        Solout[2] = 0.0;
//        return;
//    }
    
    
    
}

/** @brief mean stress */
REAL TPZPoroPermCoupling::p(TPZFMatrix<REAL> T){
    
#ifdef PZDEBUG
    if (T.Rows() != 2 && T.Cols() != 2) {
        DebugStop();
    }
#endif
    
    REAL mean_stress = 0.0;
    mean_stress = (T(0,0) + T(1,1) + T(1,1))/3;
    return mean_stress;
}

/** @brief mean stress */
TPZFMatrix<REAL> TPZPoroPermCoupling::s(TPZFMatrix<REAL> T){

#ifdef PZDEBUG
    if (T.Rows() != 2 && T.Cols() != 2) {
        DebugStop();
    }
#endif
    
    REAL  mean_stress = p(T);
    TPZFMatrix<REAL> H = T;
    TPZFNMatrix<6,REAL> I(2,2,0.0);
    I.Identity();
    H = T - mean_stress * I;
    return H;
}

/** @brief J2 invariant stress */
REAL TPZPoroPermCoupling::J2(TPZFMatrix<REAL> T){
    
#ifdef PZDEBUG
    if (T.Rows() != 2 && T.Cols() != 2) {
        DebugStop();
    }
#endif
    TPZFMatrix<REAL> S = T;
    TPZFMatrix<REAL> S_inner = T;
    S.Transpose(&S);
    S.Multiply(T,S_inner);
    
    REAL j2 = 0.5*(S_inner(0,0) + S_inner(1,1) + S_inner(1,1));
    
    return j2
    ;
}

/** @brief J3 invariant stress */
REAL TPZPoroPermCoupling::J3(TPZFMatrix<REAL> T){
    
#ifdef PZDEBUG
    if (T.Rows() != 2 && T.Cols() != 2) {
        DebugStop();
    }
#endif
    TPZFMatrix<REAL> Tinv = T;
    REAL det = T(0,0)*T(1,1)-T(1,0)*T(0,1);
    Tinv.Resize(3, 3);
    Tinv(2,2) = Tinv(1,1);
    Tinv.DeterminantInverse(det, Tinv);
    
    return det;
    
}

/** @brief theta */
REAL TPZPoroPermCoupling::theta(TPZFMatrix<REAL> T){
    
#ifdef PZDEBUG
    if (T.Rows() != 2 && T.Cols() != 2) {
        DebugStop();
    }
#endif
    
    REAL theta;
    REAL arg = -3.0*sqrt(3.0)*J3(s(T))/(2.0*pow(J2(s(T)), 1.5)) + 1.0e-14;
    theta = (1.0/3.0)*asin(arg);
    
    return theta;
}

/** @brief Phi Mohr-Coulomb */
REAL TPZPoroPermCoupling::Phi_MC(TPZFMatrix<REAL> T){

#ifdef PZDEBUG
    if (T.Rows() != 2 && T.Cols() != 2) {
        DebugStop();
    }
#endif
    
    REAL theta_v = theta(T);
    
    REAL phi = (cos(theta_v) - (1.0/sqrt(3.0))*sin(theta_v)*sin(fphi_f)) * sqrt(J2(s(T))) + p(T)*sin(fphi_f) - fc * cos(fphi_f);
    return phi;
}

/** @brief Phi Drucker-Prager */
REAL TPZPoroPermCoupling::Phi_DP(TPZFMatrix<REAL> T){
    
#ifdef PZDEBUG
    if (T.Rows() != 2 && T.Cols() != 2) {
        DebugStop();
    }
#endif
    REAL eta = 6.0*(sin(fphi_f))/(sqrt(3.0)*(3.0-sin(fphi_f)));
    REAL xi = 6.0*(cos(fphi_f))/(sqrt(3.0)*(3.0-sin(fphi_f)));
    REAL phi = sqrt(J2(s(T))) + eta *  p(T) - xi * fc ;
    return phi;

}

/** @brief plasticity multiplier delta_gamma */
REAL TPZPoroPermCoupling::Phi_tilde_DP(TPZFMatrix<REAL> T, REAL d_gamma_guest){
    
#ifdef PZDEBUG
    if (T.Rows() != 2 && T.Cols() != 2) {
        DebugStop();
    }
#endif
    REAL phi = sqrt(J2(s(T))) - fmu * d_gamma_guest + feta_dp *  (p(T) - fK * feta_dp * d_gamma_guest) - fxi_dp * fc ;
    return phi;
    
}

/** @brief plasticity multiplier delta_gamma */
REAL TPZPoroPermCoupling::Phi_tilde_DP_delta_gamma(TPZFMatrix<REAL> T, REAL d_gamma_guest){
    
#ifdef PZDEBUG
    if (T.Rows() != 2 && T.Cols() != 2) {
        DebugStop();
    }
#endif

    REAL d_phi_d_delta_gamma = - fmu - fK * feta_dp * feta_dp;
    
    if (fabs(d_phi_d_delta_gamma) <= 1.0e-18) {
        d_phi_d_delta_gamma = 1.0e-18;
    }
    
    return d_phi_d_delta_gamma;
    
}

/** @brief plasticity multiplier delta_gamma using newton iterations */
REAL TPZPoroPermCoupling::delta_gamma_finder(TPZFMatrix<REAL> T, REAL d_gamma_guest){
    
    REAL tol = 1.0e-10;
    REAL error = 1.0;
    int n_iter = 20;
    REAL d_gamma_converged = d_gamma_guest;
    
    for (int i = 0; i < n_iter; i++) {
        d_gamma_converged = d_gamma_converged - Phi_tilde_DP(T,d_gamma_converged) / Phi_tilde_DP_delta_gamma(T,d_gamma_converged);
        error = Phi_tilde_DP(T,d_gamma_converged);
        if (error <= tol) {
            break;
        }
        
    }
    
    return d_gamma_converged;
    
}

/** @brief Drucker prager strain update */
TPZFMatrix<REAL> TPZPoroPermCoupling::strain_DP(TPZFMatrix<REAL> T){
    DebugStop();
}

/** @brief Drucker prager stress update */
TPZFMatrix<REAL> TPZPoroPermCoupling::stress_DP(TPZFMatrix<REAL> T){
    DebugStop();
}

/** @brief Drucker prager elastoplastic corrector  */
void TPZPoroPermCoupling::corrector_DP(TPZFMatrix<REAL> Grad_u_n, TPZFMatrix<REAL> Grad_u, TPZFMatrix<REAL> &e_e, TPZFMatrix<REAL> &e_p, TPZFMatrix<REAL> &S){
    
#ifdef PZDEBUG
    if (Grad_u_n.Rows() != 2 && Grad_u_n.Cols() != 2) {
        DebugStop();
    }

    if (Grad_u.Rows() != 2 && Grad_u.Cols() != 2) {
        DebugStop();
    }

    if (e_p.Rows() != 2 && e_p.Cols() != 2) {
        DebugStop();
    }

    if (e_e.Rows() != 2 && e_e.Cols() != 2) {
        DebugStop();
    }
#endif
    
   
    TPZFMatrix<REAL> Grad_du, Grad_du_Transpose, delta_e;
    
    Grad_du = Grad_u - Grad_u_n;
    Grad_du.Transpose(&Grad_du_Transpose);
    delta_e = Grad_du + Grad_du_Transpose;
    delta_e *= 0.5;
    
    
    TPZFMatrix<REAL> e_t, e_trial;
    TPZFMatrix<REAL> s_trial, I(e_t.Rows(),e_t.Cols(),0.0);
    I.Identity();
    
    /** Trial strain */
    e_t = e_e + e_p;
    e_trial = e_t + delta_e;

    /** Trial stress */
    REAL trace = (e_trial(0,0) + e_trial(1,1)); // @omar:: this is ok since in plane strain modelling  e_trial(2,2) = 0.0
    s_trial = 2.0 * fmu * e_trial + flambda * trace * I;
    
    if (Phi_DP(s_trial) <= 0.0) {
        /** Elastic update */
        e_e = e_trial;
        S = s_trial;
    }
    else{
        
        /** Plastic update */
        DebugStop();
    }
    
    
}