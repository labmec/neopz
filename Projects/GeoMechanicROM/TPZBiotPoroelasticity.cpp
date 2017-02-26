//
//  TPZBiotPoroelasticity.cpp
//  PZ
//
//  Created by Omar on 2/25/17.
//
//

#include "TPZBiotPoroelasticity.h"


/** @brief Default constructor */
TPZBiotPoroelasticity::TPZBiotPoroelasticity() : TPZMatWithMem<TPZPoroPermMemory, TPZDiscontinuousGalerkin>(){
    
    fdimension = 0;
    fnon_symetric = 0.0;
    
}

/** @brief Constructor based on a material id */
TPZBiotPoroelasticity::TPZBiotPoroelasticity(int matid, int dimension) : TPZMatWithMem<TPZPoroPermMemory, TPZDiscontinuousGalerkin>(matid){
    
    fdimension = dimension;
    fnon_symetric = 0.0;
}

/** @brief Constructor based on a Biot Poroelasticity  object */
TPZBiotPoroelasticity::TPZBiotPoroelasticity(const TPZBiotPoroelasticity &mat) : TPZMatWithMem<TPZPoroPermMemory, TPZDiscontinuousGalerkin>(mat){
    
    this->fdimension    = mat.fdimension;
    this->fnon_symetric = mat.fnon_symetric;
    
}


/** @brief Default destructor */
TPZBiotPoroelasticity::~TPZBiotPoroelasticity(){
    
}

/** @brief Set the required data at each integration point */
void TPZBiotPoroelasticity::FillDataRequirements(TPZVec<TPZMaterialData> &datavec){
    
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

/** @brief Set the required data at each integration point */
void TPZBiotPoroelasticity::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec){
    
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}


/** print out the data associated with the material */
void TPZBiotPoroelasticity::Print(std::ostream &out){
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

/** returns the variable index associated with the name */
int TPZBiotPoroelasticity::VariableIndex(const std::string &name){
    
    //	Elasticity Variables
    if(!strcmp("u",name.c_str()))              return	0;
    
    //	Diffusion Variables
    if(!strcmp("p_ex",name.c_str()))			return	1;
    if(!strcmp("v",name.c_str()))				return	2;
    
    //	Elasticity Variables
    if(!strcmp("ue",name.c_str()))              return	3;
    
    //	Diffusion Variables
    if(!strcmp("pe_ex",name.c_str()))			return	4;
    if(!strcmp("ve",name.c_str()))				return	5;
    
}

/** returns the number of variables associated with the variable
 indexed by var.  var is obtained by calling VariableIndex */
int TPZBiotPoroelasticity::NSolutionVariables(int var){
    
    switch(var) {
        case 0:
            return fdimension; // Vector
        case 1:
            return 1; // Scalar
        case 2:
            return fdimension; // Vector
        case 3:
            return fdimension; // Vector
        case 4:
            return 1; // Scalar
        case 5:
            return fdimension; // Vector
    }
    return TPZMatWithMem::NSolutionVariables(var);
    
}

void TPZBiotPoroelasticity::Compute_Sigma(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_u){
    
    TPZFNMatrix<6,REAL> Grad_ut(3,3,0.0), epsilon(3,3,0.0), I(3,3,0.0);
    Grad_u.Transpose(&Grad_ut);
    
    epsilon = Grad_u + Grad_ut;
    epsilon *= 0.5;
    
    I.Identity();
    
    REAL trace = (epsilon(0,0) + epsilon(1,1));
    S = 2.0 * fmu * epsilon + flambda * trace * I;
    
}

// Contribute Methods being used
void TPZBiotPoroelasticity::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){

    int u_b = 0;
    int p_b = 1;
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[u_b].fShapeType;
    if(shapetype == datavec[u_b].EVecShape)
    {
        DebugStop();
//        this->ContributeVec(datavec, weight, ek, ef);
        return;
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
    
    
    // Transformations
    TPZFNMatrix<27,REAL> grad_phi_u;
    TPZAxesTools<STATE>::Axes2XYZ(dphiu, grad_phi_u, axes_u);
    
    TPZFNMatrix<3,REAL> grad_u;
    TPZAxesTools<STATE>::Axes2XYZ(du, grad_u, axes_u);
    
    TPZFNMatrix<9,REAL> grad_phi_p;
    TPZAxesTools<STATE>::Axes2XYZ(dphip, grad_phi_p, axes_p);
    
    TPZFNMatrix<3,REAL> grad_p;
    TPZAxesTools<STATE>::Axes2XYZ(dp, grad_p, axes_p);
    
    int nphi_u = phiu.Rows();
    int nphi_p = phip.Rows();
    
    int first_u = 0;
    int first_p = 2*nphi_u;
    
    REAL dt = fSimulationData->dt();
    REAL div_u = grad_u(0,0) + grad_u(1,1);
    
    TPZFNMatrix<6,REAL> Grad_u(3,3,0.0);
    TPZFNMatrix<9,REAL> S;
    grad_u.Resize(3, 3);
    this->Compute_Sigma(S,grad_u);
    
    TPZFNMatrix<9,REAL> Grad_vx_i(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_i(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vx_j(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_j(fdimension,1,0.0);
    
    if (!fSimulationData->IsCurrentStateQ()) {
        
        // Darcy mono-phascis flow
        for (int ip = 0; ip < nphi_p; ip++) {
            ef(ip + first_p, 0)		+= weight *  (-1.0) * (1.0/dt) * (falpha * div_u + fSe * p[0]) * phip(ip,0);
        }
        return;
    }
    
    for (int iu = 0; iu < nphi_u; iu++) {
        
        // Computing Gradient of the test function for each component
        for (int d = 0; d < fdimension; d++) {
            Grad_vx_i(d,0) = grad_phi_u(d,iu);
            Grad_vy_i(d,0) = grad_phi_u(d,iu);
        }
        
        ef(2*iu + first_u, 0)   += weight * ((S(0,0) - falpha * p[0]) * Grad_vx_i(0,0) + S(0,1) * Grad_vx_i(1,0));
        ef(2*iu+1 + first_u, 0)	+= weight * (S(1,0) * Grad_vy_i(0,0) + (S(1,1) - falpha * p[0]) * Grad_vy_i(1,0));
        
        for (int ju = 0; ju < nphi_u; ju++) {
            
            // Computing Gradient of the test function for each component
            for (int d = 0; d < fdimension; d++) {
                Grad_vx_j(d,0) = grad_phi_u(d,ju);
                Grad_vy_j(d,0) = grad_phi_u(d,ju);
            }
            
            ek(2*iu + first_u, 2*ju + first_u)      += weight * ( ( (2.0*fmu + flambda) * Grad_vx_j(0,0) ) * Grad_vx_i(0,0) + fmu * Grad_vx_j(1,0) * Grad_vx_i(1,0));
            ek(2*iu + first_u, 2*ju+1 + first_u)    += weight * ( (flambda * Grad_vy_j(1,0) ) * Grad_vx_i(0,0) + fmu * Grad_vy_j(0,0) * Grad_vx_i(1,0)  );
            ek(2*iu+1 + first_u, 2*ju + first_u)	+= weight * ( fmu * Grad_vx_j(1,0) * Grad_vy_i(0,0) + flambda * Grad_vx_j(0,0) * Grad_vy_i(1,0));
            ek(2*iu+1 + first_u, 2*ju+1 + first_u)	+= weight * ( (2.0*fmu + flambda) * Grad_vy_j(1,0) * Grad_vy_i(1,0) + fmu * Grad_vy_j(0,0) * Grad_vy_i(0,0) );
        }
        
        for(int jp = 0; jp < nphi_p; jp++)
        {
            ek(2*iu,first_p+jp) +=      (-1.0)* weight * falpha * phip(jp,0) * Grad_vx_i(0,0);
            ek(2*iu + 1,first_p+jp) +=  (-1.0)* weight * falpha * phip(jp,0) * Grad_vy_i(1,0);
        }
        
    }
    
    REAL c = fk/feta;
    // Darcy mono-phasic flow
    for (int ip = 0; ip < nphi_p; ip++) {
        
        REAL dot = 0.0;
        for (int i = 0;  i < fdimension; i++) {
            dot += grad_p(i,0) * grad_phi_p(i,ip);
        }
        
        ef(ip + first_p, 0)		+= 1.0 * weight * (c * dot + (1.0/dt) * (falpha * div_u + fSe * p[0]) * phip(ip,0) );
        
        //	Coupling matrix
        for(int ju = 0; ju < nphi_u; ju++ )
        {
            // Computing Gradient of the test function for each component
            for (int d = 0; d < fdimension; d++) {
                Grad_vx_j(d,0) = grad_phi_u(d,ju);
                Grad_vy_j(d,0) = grad_phi_u(d,ju);
            }
            
            ek(ip + first_p,2*ju)     += (1.0)* (1.0/dt)  * weight * falpha * phip(ip,0) * Grad_vx_j(0,0);
            ek(ip + first_p,2*ju+1)   += (1.0)* (1.0/dt)  * weight * falpha * phip(ip,0) * Grad_vy_j(1,0);
        }
        
        for (int jp = 0; jp < nphi_p; jp++) {
            
            REAL dot = 0.0;
            for (int i = 0;  i < fdimension; i++) {
                dot += grad_phi_p(i,jp) * grad_phi_p(i,ip);
            }
            
            ek(ip + first_p, jp + first_p)  += 1.0 * weight * ( c * dot + (1.0/dt) * (fSe * phip(jp,0)) * phip(ip,0) );
        }
        
    }
    
}


void TPZBiotPoroelasticity::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){

    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
}

void TPZBiotPoroelasticity::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){

    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int u_b = 0;
    int p_b = 1;
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[u_b].fShapeType;
    if(shapetype == datavec[u_b].EVecShape)
    {
        DebugStop();
        //ContributeVecBC(datavec, weight, ek, ef, bc);
        return;
    }
    
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
    REAL dt  = this->SimulationData()->dt();
    REAL Value = bc.Val2()(0,0);
    if (bc.HasTimedependentBCForcingFunction()) {
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
                ef(2*in,0)		+= -1.0* v[0]*phiu(in,0)*weight;		//	Tnx
                ef(2*in+1,0)	+= -1.0* v[1]*phiu(in,0)*weight;		//	Tny
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
                ef(2*in,0)		+= -1.0 * v[0]*phiu(in,0)*weight;		//	Tnx
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
                ef(2*in+1,0)	+= -1.0 *  v[1]*phiu(in,0)*weight;		//	Tny
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
                ef(2*in,0)		+= -1.0 * v[0]*phiu(in,0)*weight;		//	Tnx
                ef(2*in+1,0)	+= -1.0 * v[1]*phiu(in,0)*weight;		//	Tny
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
                ef(2*in,0)		+= -1.0 * v[0]*phiu(in,0)*weight;		//	Tnx
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0 * v[2]*phip(in,0)*weight;	// Qnormal
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
                ef(2*in+1,0)	+= -1.0 * v[1]*phiu(in,0)*weight;		//	Tny
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


/** returns the solution associated with the var index based on the finite element approximation */
void TPZBiotPoroelasticity::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
    
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
    REAL to_Darcy   = 1.013249966e+15; //md
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_p(3,1,0.0),Grad_u(3,3,0.0),Grad_u_n(3,3,0.0),e_e(3,3,0.0),e_p(3,3,0.0),S;
    
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0); // dp/dx
    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1); // dp/dy
    
    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
    
    Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    
    //	Displacements
    if(var == 0){
        Solout[0] = u[0];
        Solout[1] = u[1];
        return;
    }
    
    //	pressure
    if(var == 1) {
        Solout[0] = p[0]*to_Mpa;
        return;
    }
    
    //	darcy
    if(var == 2) {
        Solout[0] = -1.0*(fk/feta)*Grad_p(0,0);
        Solout[0] = -1.0*(fk/feta)*Grad_p(1,0);
        return;
    }
    
    // displacement exact
    if(var == 3) {
        
        TPZManVector<STATE,5> f(1,0.0);
        TPZFNMatrix<4,STATE> df(4,1,0.0);
        if (this->HasTimedependentForcingFunction()) {
            REAL time = fSimulationData->t();
            this->fTimeDependentForcingFunction->Execute(datavec[u_b].x, time, f, df);
        }
        
        Solout[0] = f[0];
        Solout[1] = f[1];
        
        return;
    }
    
    // Pressure exact
    if(var == 4) {
        
        TPZManVector<STATE,5> f(1,0.0);
        TPZFNMatrix<4,STATE> df(4,1,0.0);
        if (this->HasTimedependentForcingFunction()) {
            REAL time = fSimulationData->t();
            this->fTimeDependentForcingFunction->Execute(datavec[u_b].x, time, f, df);
        }
        
        Solout[0] = f[2]*to_Mpa;
        
        return;
    }
    
    //	darcy exact
    if(var == 5) {
        TPZManVector<STATE,5> f(1,0.0);
        TPZFNMatrix<4,STATE> df(4,1,0.0);
        if (this->HasTimedependentForcingFunction()) {
            REAL time = fSimulationData->t();
            this->fTimeDependentForcingFunction->Execute(datavec[u_b].x, time, f, df);
        }
        
        Solout[0] = f[3];
        Solout[1] = f[4];
        return;
    }
    

    
    std::cout  << "not implemented. " << std::endl;
    DebugStop();
    
}

/** @brief Unique identifier for serialization purposes */
int TPZBiotPoroelasticity::ClassId() const {
    return -6378637866;
}

/** @brief Save object data to a stream */
void TPZBiotPoroelasticity::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
}

/** @brief Read object data from a stream */
void TPZBiotPoroelasticity::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}
