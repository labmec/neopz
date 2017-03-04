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
    fMFsymetric = 1.0;
    fCsymetric = 1.0;
    fIsMixedQ = false;
    fIsSymmetricQ = false;
    
}

/** @brief Constructor based on a material id */
TPZBiotPoroelasticity::TPZBiotPoroelasticity(int matid, int dimension) : TPZMatWithMem<TPZPoroPermMemory, TPZDiscontinuousGalerkin>(matid){
    
    fdimension = dimension;
    fMFsymetric = 1.0;
    fCsymetric = 1.0;
    fIsMixedQ = false;
    fIsSymmetricQ = false;
}

/** @brief Constructor based on a Biot Poroelasticity  object */
TPZBiotPoroelasticity::TPZBiotPoroelasticity(const TPZBiotPoroelasticity &mat) : TPZMatWithMem<TPZPoroPermMemory, TPZDiscontinuousGalerkin>(mat){
    
    this->fdimension    = mat.fdimension;
    this->fMFsymetric = mat.fMFsymetric;
    this->fCsymetric = mat.fCsymetric;
    this->fIsMixedQ     = mat.fIsMixedQ;
    this->fIsSymmetricQ     = mat.fIsSymmetricQ;
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
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[0].fShapeType;
    if(shapetype == datavec[0].EVecShape)
    {
        // RB case
        datavec[0].fNeedsBasis = false;
        datavec[0].fNeedsSol = false;
        return;
    }
    
}

/** @brief Set the required data at each integration point */
void TPZBiotPoroelasticity::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec){
    
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsBasis = true;
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
    
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[0].fShapeType;
    if(shapetype == datavec[0].EVecShape)
    {
        // RB case
        datavec[0].fNeedsBasis = false;
        datavec[0].fNeedsSol = false;
        datavec[0].fNeedsNormal = false;
        return;
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
    
    return -1;
    
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

void TPZBiotPoroelasticity::Compute_Sigma_fast(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_u){
    
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

// Contribute Methods being used
void TPZBiotPoroelasticity::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){

    if (fIsSymmetricQ) {
        this->fMFsymetric = fSimulationData->dt();
        this->fCsymetric  = -1.0*fSimulationData->dt();
    }
    
    if (fIsMixedQ) {
        ContributeMF(datavec, weight, ek, ef);
        return;
    }
    
    int u_b = 0;
    int p_b = 1;
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[u_b].fShapeType;
    if(shapetype == datavec[u_b].EVecShape){
        this->ContributeRB(datavec, weight, ek, ef);
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
    
    TPZFNMatrix <15,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    // Transformations
    TPZFNMatrix<27,REAL> grad_phi_u;
    TPZFNMatrix<3,REAL> grad_u;
    
    TPZAxesTools<STATE>::Axes2XYZ(dphiu, grad_phi_u, axes_u);
    TPZAxesTools<STATE>::Axes2XYZ(du, grad_u, axes_u);
    
    TPZFNMatrix<9,REAL> grad_phi_p;
    TPZAxesTools<STATE>::Axes2XYZ(dphip, grad_phi_p, axes_p);
    
    TPZFNMatrix<3,REAL> grad_p;
    TPZAxesTools<STATE>::Axes2XYZ(dp, grad_p, axes_p);
    
    int nphi_u = phiu.Rows();
    int nphi_p = phip.Rows();
    
    int first_u = 0;
    int first_p = first_u + 2*nphi_u;
    
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
            ef(ip + first_p, 0)		+= (fCsymetric) * weight *  (-1.0) * (1.0/dt) * (falpha * div_u + fSe * p[0]) * phip(ip,0);
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
        
        ef(ip + first_p, 0)		+= (fCsymetric) * weight * (c * dot + (1.0/dt) * (falpha * div_u + fSe * p[0]) * phip(ip,0) );
        
        //	Coupling matrix
        for(int ju = 0; ju < nphi_u; ju++ )
        {

            // Computing Gradient of the test function for each component
            for (int d = 0; d < fdimension; d++) {
                Grad_vx_j(d,0) = grad_phi_u(d,ju);
                Grad_vy_j(d,0) = grad_phi_u(d,ju);
            }
            
            ek(ip + first_p,2*ju)     += (fCsymetric) * (1.0/dt)  * weight * falpha * phip(ip,0) * Grad_vx_j(0,0);
            ek(ip + first_p,2*ju+1)   += (fCsymetric) * (1.0/dt)  * weight * falpha * phip(ip,0) * Grad_vy_j(1,0);
        }
        
        for (int jp = 0; jp < nphi_p; jp++) {
            
            REAL dot = 0.0;
            for (int i = 0;  i < fdimension; i++) {
                dot += grad_phi_p(i,jp) * grad_phi_p(i,ip);
            }
            
            ek(ip + first_p, jp + first_p)  += (fCsymetric) * weight * ( c * dot + (1.0/dt) * (fSe * phip(jp,0)) * phip(ip,0) );
        }
        
    }
    
}


void TPZBiotPoroelasticity::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){

    if (fIsMixedQ) {
        ContributeMF(datavec, weight,ef);
        return;
    }
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
}

void TPZBiotPoroelasticity::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){

   
    if (fIsSymmetricQ) {
        this->fMFsymetric = fSimulationData->dt();
        this->fCsymetric  = -1.0*fSimulationData->dt();
    }
    
    if (fIsMixedQ) {
        this->ContributeMFBC(datavec, weight, ek, ef, bc);
        return;
    }
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int u_b = 0;
    int p_b = 1;

    int ux_id = 0;
    int uy_id = 0;
    TPZMaterialData::MShapeFunctionType shapetype = datavec[u_b].fShapeType;
    REAL c_big = 1.0;
    if(shapetype == datavec[u_b].EVecShape)
    {
        this->ContributeRB_BC(datavec, weight, ek, ef, bc);
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
                ef(2*in,0)		+= c_big*gBigNumber*(u[0] - v[0])*phiu(in,ux_id)*weight;	// X displacement Value
                ef(2*in+1,0)	+= c_big*gBigNumber*(u[1] - v[1])*phiu(in,uy_id)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= c_big*gBigNumber*phiu(in,0)*phiu(jn,ux_id)*weight;	// X displacement
                    ek(2*in+1,2*jn+1)	+= c_big*gBigNumber*phiu(in,0)*phiu(jn,uy_id)*weight;	// Y displacement
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
                ef(2*in,0)		+= c_big*gBigNumber*(u[0] - v[0])*phiu(in,ux_id)*weight;	// X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= c_big*gBigNumber*phiu(in,0)*phiu(jn,ux_id)*weight;	// X displacement
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
                ef(2*in+1,0)	+= c_big*gBigNumber*(u[1] - v[1])*phiu(in,uy_id)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)	+= c_big*gBigNumber*phiu(in,0)*phiu(jn,uy_id)*weight;	// Y displacement
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
                ef(2*in,0)		+= -1.0* v[0]*phiu(in,ux_id)*weight;		//	Tnx
                ef(2*in+1,0)	+= -1.0* v[1]*phiu(in,uy_id)*weight;		//	Tny
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
                ef(2*in,0)		+= -1.0 * v[0]*phiu(in,ux_id)*weight;		//	Tnx
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
                ef(2*in+1,0)	+= -1.0 * v[1]*phiu(in,uy_id)*weight;		//	Tny
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
                ef(2*in,0)		+= c_big*gBigNumber*(u[0] - v[0])*phiu(in,ux_id)*weight;	// X displacement Value
                ef(2*in+1,0)	+= c_big*gBigNumber*(u[1] - v[1])*phiu(in,uy_id)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= c_big*gBigNumber*phiu(in,0)*phiu(jn,ux_id)*weight;	// X displacement
                    ek(2*in+1,2*jn+1)	+= c_big*gBigNumber*phiu(in,0)*phiu(jn,uy_id)*weight;	// Y displacement
                }
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0 * (fCsymetric) * v[2]*phip(in,0)*weight;	// Qnormal
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
                ef(2*in,0)		+= c_big*gBigNumber*(u[0] - v[0])*phiu(in,ux_id)*weight;	// X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= c_big*gBigNumber*phiu(in,0)*phiu(jn,ux_id)*weight;	// X displacement
                }
            }
            
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0 * (fCsymetric) * v[2]*phip(in,0)*weight;	// Qnormal
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
                ef(2*in+1,0)	+= c_big*gBigNumber*(u[1] - v[1])*phiu(in,uy_id)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)	+= c_big*gBigNumber*phiu(in,0)*phiu(jn,uy_id)*weight;	// Y displacement
                }
            }
            
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0 * (fCsymetric) * v[2]*phip(in,0)*weight;	// Qnormal
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
                ef(2*in,0)		+= -1.0 * v[0]*phiu(in,ux_id)*weight;		//	Tnx
                ef(2*in+1,0)	+= -1.0 * v[1]*phiu(in,uy_id)*weight;		//	Tny
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0 * (fCsymetric) * v[2]*phip(in,0)*weight;	// Qnormal
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
                ef(2*in,0)		+= -1.0 * v[0]*phiu(in,ux_id)*weight;		//	Tnx
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0  * (fCsymetric) * v[2]*phip(in,0)*weight;	// Qnormal
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
                ef(2*in+1,0)	+= -1.0 * v[1]*phiu(in,uy_id)*weight;		//	Tny
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in+2*phru,0)	+= -1.0 * (fCsymetric) * v[2]*phip(in,0)*weight;	// Qnormal
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


// Reduce basis methods
void TPZBiotPoroelasticity::ContributeRB(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int u_b = 0;
    int p_b = 1;
    
    // RB functions must to include axes transformations .
    
    // Getting RB functions and solution form integration points
    // Get the data at the integrations points

    long global_point_index = datavec[u_b].intGlobPtIndex;
    TPZPoroPermMemory &point_memory = GetMemory()[global_point_index];
    TPZFMatrix<REAL> & phiu = point_memory.phi_u();
    TPZFMatrix<REAL> & dphiu = point_memory.grad_phi_u();
    
    TPZFNMatrix<3,REAL> & int_u = point_memory.u();
    TPZFNMatrix<3,REAL> & int_u_n = point_memory.u_n();
    
    TPZFNMatrix<9,REAL> & int_grad_u = point_memory.grad_u();
    TPZFNMatrix<9,REAL> & int_grad_u_n = point_memory.grad_u_n();
    
    // Getting the space functions
//    TPZFMatrix<REAL>    &phiu   =   datavec[u_b].phi;
    TPZFMatrix<REAL>    &phip   =   datavec[p_b].phi;
    
//    TPZFMatrix<REAL>    &dphiu   =   datavec[u_b].dphix;
    TPZFMatrix<REAL>    &dphip   =   datavec[p_b].dphix;
    
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,10> & u = datavec[u_b].sol[0];
    TPZManVector<REAL,10> & p = datavec[p_b].sol[0];
    
    TPZFNMatrix <15,REAL> & du = datavec[u_b].dsol[0];
    TPZFNMatrix <15,REAL> & dp = datavec[p_b].dsol[0];
    
//    // Transfering from integration points
//    int n_rb = int_phi_u.size();
//    phiu.Redim(n_rb, fdimension);
//    dphiu.Redim(n_rb, fdimension*fdimension);
//    for (int i = 0; i < n_rb; i++) {
//        int cd = 0;
//        for (int id = 0; id < fdimension; id++) {
//            phiu(i,id) = int_phi_u[i](id,0);
//            for (int jd = 0; jd < fdimension; jd++) {
//                dphiu(i,cd) = int_dphi_u[i](jd,id);
//                cd++;
//            }
//        }
//    }
    
    u.Resize(fdimension, 0.0);
    du.Resize(fdimension, fdimension);
    
    if (fSimulationData->IsCurrentStateQ()) {
        for (int i = 0; i <fdimension; i++) {
            u[i] = int_u_n(i,0);
            for (int j = 0; j <fdimension; j++) {
                du(i,j) = int_grad_u_n(i,j);
            }
        }
    }
    else{
        for (int i = 0; i <fdimension; i++) {
            u[i] = int_u(i,0);
            for (int j = 0; j <fdimension; j++) {
                du(i,j) = int_grad_u(i,j);
            }
        }
    }
    
    // Transformations
    TPZFNMatrix<27,REAL> grad_phi_u;
    TPZFNMatrix<9,REAL> grad_u(3,3,0.0);
    
    grad_u(0,0) = du(0,0); // dux/dx
    grad_u(0,1) = du(1,0); // dux/dy
    
    grad_u(1,0) = du(0,1); // duy/dx
    grad_u(1,1) = du(1,1); // duy/dy
    
    TPZFNMatrix<3,REAL> grad_p(fdimension,1,0.0);
    // Computing Gradient of the test function for each component
    grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0); // dp/dx
    grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1); // dp/dy
    
    int nphi_u = phiu.Rows();
    int nphi_p = phip.Rows();
    
    int first_u = 0;
    int first_p = first_u + nphi_u;
    
    REAL dt = fSimulationData->dt();
    REAL div_u = grad_u(0,0) + grad_u(1,1);
    
    TPZFNMatrix<9,REAL> S(3,3,0.0);
    this->Compute_Sigma_fast(S,grad_u);
    
    TPZFNMatrix<9,REAL> Grad_vx_i(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_i(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vx_j(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_j(fdimension,1,0.0);
    
    if (!fSimulationData->IsCurrentStateQ()) {
        
        // Darcy mono-phascis flow
        for (int ip = 0; ip < nphi_p; ip++) {
            ef(ip + first_p, 0)		+= (fCsymetric) * weight *  (-1.0) * (1.0/dt) * (falpha * div_u + fSe * p[0]) * phip(ip,0);
        }
        return;
    }
    
    
    for (int iu = 0; iu < nphi_u; iu++) {

//        // Computing Gradient of the test function for each component
//        Grad_vx_i(0,0) = dphiu(iu,0)*axes_u(0,0)+dphiu(iu,1)*axes_u(1,0); // dvx/dx
//        Grad_vx_i(1,0) = dphiu(iu,0)*axes_u(0,1)+dphiu(iu,1)*axes_u(1,1); // dvx/dy
//        
//        Grad_vy_i(0,0) = dphiu(iu,2)*axes_u(0,0)+dphiu(iu,3)*axes_u(1,0); // dvy/dx
//        Grad_vy_i(1,0) = dphiu(iu,2)*axes_u(0,1)+dphiu(iu,3)*axes_u(1,1); // dvy/dy
        
        // Computing Gradient of the test function for each component
        Grad_vx_i(0,0) = dphiu(iu,0); // dvx/dx
        Grad_vx_i(1,0) = dphiu(iu,1); // dvx/dy
        
        Grad_vy_i(0,0) = dphiu(iu,2); // dvy/dx
        Grad_vy_i(1,0) = dphiu(iu,3); // dvy/dy
        
        REAL iterm1 = weight * ((S(0,0) - falpha * p[0]) * Grad_vx_i(0,0) + S(0,1) * Grad_vx_i(1,0));
        REAL iterm2 = weight * (S(1,0) * Grad_vy_i(0,0) + (S(1,1) - falpha * p[0]) * Grad_vy_i(1,0));
        ef(iu + first_u, 0) += iterm1 + iterm2;
        
        for (int ju = 0; ju < nphi_u; ju++) {

//            // Computing Gradient of the test function
//            Grad_vx_j(0,0) = dphiu(ju,0)*axes_u(0,0)+dphiu(ju,1)*axes_u(1,0); // dvx/dx
//            Grad_vx_j(1,0) = dphiu(ju,0)*axes_u(0,1)+dphiu(ju,1)*axes_u(1,1); // dvx/dy
//            
//            Grad_vy_j(0,0) = dphiu(ju,2)*axes_u(0,0)+dphiu(ju,3)*axes_u(1,0); // dvy/dx
//            Grad_vy_j(1,0) = dphiu(ju,2)*axes_u(0,1)+dphiu(ju,3)*axes_u(1,1); // dvy/dy
            
            // Computing Gradient of the test function for each component
            Grad_vx_j(0,0) = dphiu(ju,0); // dvx/dx
            Grad_vx_j(1,0) = dphiu(ju,1); // dvx/dy
            
            Grad_vy_j(0,0) = dphiu(ju,2); // dvy/dx
            Grad_vy_j(1,0) = dphiu(ju,3); // dvy/dy
            
            REAL equ1_1 = weight * ( ( (2.0*fmu + flambda) * Grad_vx_j(0,0) ) * Grad_vx_i(0,0) + fmu * Grad_vx_j(1,0) * Grad_vx_i(1,0));
            REAL equ1_2 = weight * ( (flambda * Grad_vy_j(1,0) ) * Grad_vx_i(0,0) + fmu * Grad_vy_j(0,0) * Grad_vx_i(1,0)  );
            REAL equ2_1 = weight * ( fmu * Grad_vx_j(1,0) * Grad_vy_i(0,0) + flambda * Grad_vx_j(0,0) * Grad_vy_i(1,0));
            REAL equ2_2 = weight * ( (2.0*fmu + flambda) * Grad_vy_j(1,0) * Grad_vy_i(1,0) + fmu * Grad_vy_j(0,0) * Grad_vy_i(0,0) );
            
            ek(iu + first_u, ju + first_u) += equ1_1 + equ1_2 + equ2_1 + equ2_2;
        }
        
        for(int jp = 0; jp < nphi_p; jp++)
        {
            REAL jterm1     = (-1.0)* weight * falpha * phip(jp,0) * Grad_vx_i(0,0);
            REAL jterm2     =  (-1.0)* weight * falpha * phip(jp,0) * Grad_vy_i(1,0);
            ek(iu+first_u,jp+first_p) += jterm1 + jterm2;
        }
        
    }
    
    TPZFNMatrix<9,REAL> grad_phi_p(fdimension,nphi_p);
//    TPZAxesTools<STATE>::Axes2XYZ(dphip, grad_phi_p, axes_p);// @omar::Prohibitive at contribute
    for (int  ip = 0; ip < nphi_p; ip++) {
        grad_phi_p(0,ip) = dphip(0,ip)*axes_p(0,0)+dphip(1,ip)*axes_p(1,0); // dphi_p/dx
        grad_phi_p(1,ip) = dphip(0,ip)*axes_p(0,1)+dphip(1,ip)*axes_p(1,1); // dphi_p/dy
    }
    
    REAL c = fk/feta;
    // Darcy mono-phasic flow
    for (int ip = 0; ip < nphi_p; ip++) {
        
        REAL dot = 0.0;
        for (int i = 0;  i < fdimension; i++) {
            dot += grad_p(i,0) * grad_phi_p(i,ip);
        }
        
        ef(ip + first_p, 0)		+= (fCsymetric) * weight * (c * dot + (1.0/dt) * (falpha * div_u + fSe * p[0]) * phip(ip,0) );
        
        //	Coupling matrix
        for(int ju = 0; ju < nphi_u; ju++ )
        {
            
            // Computing Gradient of the test function for each component
            Grad_vx_j(0,0) = dphiu(ju,0); // dvx/dx
            Grad_vx_j(1,0) = dphiu(ju,1); // dvx/dy
            
            Grad_vy_j(0,0) = dphiu(ju,2); // dvy/dx
            Grad_vy_j(1,0) = dphiu(ju,3); // dvy/dy
            
            REAL jterm1 = (fCsymetric) * (1.0/dt)  * weight * falpha * phip(ip,0) * Grad_vx_j(0,0);
            REAL jterm2 = (fCsymetric) * (1.0/dt)  * weight * falpha * phip(ip,0) * Grad_vy_j(1,0);
            
            ek(ip + first_p,ju + first_u) += jterm1 + jterm2;
        }
        
        for (int jp = 0; jp < nphi_p; jp++) {
            
            REAL dot = 0.0;
            for (int i = 0;  i < fdimension; i++) {
                dot += grad_phi_p(i,jp) * grad_phi_p(i,ip);
            }
            
            ek(ip + first_p, jp + first_p)  += (fCsymetric) * weight * ( c * dot + (1.0/dt) * (fSe * phip(jp,0)) * phip(ip,0) );
        }
        
    }
    
}

void TPZBiotPoroelasticity::ContributeRB_BC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    int u_b = 0;
    int p_b = 1;
    
    int ux_id = 0;
    int uy_id = 1;
    REAL c_big = 1.0;
    
    
    // Getting RB functions and solution form integration points
    // Get the data at the integrations points
    TPZMatWithMem<TPZPoroPermMemory,TPZBndCond>  & material_bc_mem = dynamic_cast<TPZMatWithMem<TPZPoroPermMemory,TPZBndCond > & >(bc);

    TPZFNMatrix<3,REAL>  int_u, int_u_n;
    TPZFNMatrix<9,REAL>  int_grad_u, int_grad_u_n;
    long global_point_index = datavec[u_b].intGlobPtIndex;
    TPZPoroPermMemory &point_memory = material_bc_mem.GetMemory()[global_point_index];
    TPZFMatrix<REAL> & phiu = point_memory.phi_u();
    
    int_u = point_memory.u();
    int_u_n = point_memory.u_n();
    
    int_grad_u = point_memory.grad_u();
    int_grad_u_n = point_memory.grad_u_n();
    
    
//    TPZFMatrix<REAL>  &phiu = datavec[u_b].phi;
    TPZFMatrix<REAL>  &phip = datavec[p_b].phi;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    u.Resize(fdimension, 0.0);
    
    if (fSimulationData->IsCurrentStateQ()) {
        for (int i = 0; i <fdimension; i++) {
            u[i] = int_u_n(i,0);
        }
    }
    else{
        for (int i = 0; i <fdimension; i++) {
            u[i] = int_u(i,0);
        }
    }
    
    int phru = phiu.Rows();
    int phrp = phip.Rows();
    
    int first_u = 0;
    int first_p = first_u + phru;
    
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
                REAL iterm1		 = c_big*gBigNumber*(u[0] - v[0])*phiu(in,ux_id)*weight;	// X displacement Value
                REAL iterm2      = c_big*gBigNumber*(u[1] - v[1])*phiu(in,uy_id)*weight;	// y displacement Value
                ef(in + first_u,0) += iterm1 + iterm2;
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    REAL jterm1		= c_big*gBigNumber*phiu(in,0)*phiu(jn,ux_id)*weight;	// X displacement
                    REAL jterm2		= c_big*gBigNumber*phiu(in,0)*phiu(jn,uy_id)*weight;	// Y displacement
                    ek(in + first_u,jn + first_u)   += jterm1 + jterm2;
                }
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in + first_p,0)		+= gBigNumber*(p[0]-v[2])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in + first_p,jn + first_p)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
                }
            }
            break;
        }
            
        case 1 :
        {
            //	Dirichlet condition for x state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                REAL iterm1		 = c_big*gBigNumber*(u[0] - v[0])*phiu(in,ux_id)*weight;	// X displacement Value
                ef(in + first_u,0) += iterm1;
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    REAL jterm1		= c_big*gBigNumber*phiu(in,0)*phiu(jn,ux_id)*weight;	// X displacement
                    ek(in + first_u,jn + first_u)   += jterm1;
                }
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in + first_p,0)		+= gBigNumber*(p[0]-v[2])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in + first_p,jn + first_p)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
                }
            }
            break;
        }
            
        case 2 :
        {
            //	Dirichlet condition for y state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                REAL iterm2      = c_big*gBigNumber*(u[1] - v[1])*phiu(in,uy_id)*weight;	// y displacement Value
                ef(in + first_u,0) += iterm2;
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    REAL jterm2		= c_big*gBigNumber*phiu(in,0)*phiu(jn,uy_id)*weight;	// Y displacement
                    ek(in + first_u,jn + first_u)   += jterm2;
                }
            }
            
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in + first_p,0)		+= gBigNumber*(p[0]-v[2])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in + first_p,jn + first_p)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
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
                REAL iterm1		= -1.0* v[0]*phiu(in,ux_id)*weight;		//	Tnx
                REAL iterm2		= -1.0* v[1]*phiu(in,uy_id)*weight;		//	Tny
                ef(in + first_u,0) += iterm1 + iterm2;
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in + first_p,0)		+= gBigNumber*(p[0]-v[2])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in + first_p,jn + first_p)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
                }
            }
            break;
        }
            
        case 4 :
        {
            //	Neumann condition for x state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                REAL iterm1		= -1.0* v[0]*phiu(in,ux_id)*weight;		//	Tnx
                ef(in + first_u,0) += iterm1;
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in + first_p,0)		+= gBigNumber*(p[0]-v[2])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in + first_p,jn + first_p)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
                }
            }
            break;
        }
            
        case 5 :
        {
            //	Neumann condition for y state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                REAL iterm2		= -1.0* v[1]*phiu(in,uy_id)*weight;		//	Tny
                ef(in + first_u,0) += iterm2;
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Contribution for load Vector
                ef(in + first_p,0)		+= gBigNumber*(p[0]-v[2])*phip(in,0)*weight;	// P Pressure Value
                
                for (jn = 0 ; jn < phrp; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(in + first_p,jn + first_p)		+= gBigNumber*phip(in,0)*phip(jn,0)*weight;	// P Pressure
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
                REAL iterm1		 = c_big*gBigNumber*(u[0] - v[0])*phiu(in,ux_id)*weight;	// X displacement Value
                REAL iterm2      = c_big*gBigNumber*(u[1] - v[1])*phiu(in,uy_id)*weight;	// y displacement Value
                ef(in + first_u,0) += iterm1 + iterm2;
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    REAL jterm1		= c_big*gBigNumber*phiu(in,0)*phiu(jn,ux_id)*weight;	// X displacement
                    REAL jterm2		= c_big*gBigNumber*phiu(in,0)*phiu(jn,uy_id)*weight;	// Y displacement
                    ek(in + first_u,jn + first_u)   += jterm1 + jterm2;
                }
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in + first_p,0)	+= -1.0 * (fCsymetric) * v[2]*phip(in,0)*weight;	// Qnormal
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
                REAL iterm1		 = c_big*gBigNumber*(u[0] - v[0])*phiu(in,ux_id)*weight;	// X displacement Value
                ef(in + first_u,0) += iterm1;
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    REAL jterm1		= c_big*gBigNumber*phiu(in,0)*phiu(jn,ux_id)*weight;	// X displacement
                    ek(in + first_u,jn + first_u)   += jterm1;
                }
            }
            
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in + first_p,0)	+= -1.0 * (fCsymetric) * v[2]*phip(in,0)*weight;	// Qnormal
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
                REAL iterm2      = c_big*gBigNumber*(u[1] - v[1])*phiu(in,uy_id)*weight;	// y displacement Value
                ef(in + first_u,0) += iterm2;
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    REAL jterm2		= c_big*gBigNumber*phiu(in,0)*phiu(jn,uy_id)*weight;	// Y displacement
                    ek(in + first_u,jn + first_u)   += jterm2;
                }
            }
            
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in + first_p,0)	+= -1.0 * (fCsymetric) * v[2]*phip(in,0)*weight;	// Qnormal
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
                REAL iterm1		= -1.0* v[0]*phiu(in,ux_id)*weight;		//	Tnx
                REAL iterm2		= -1.0* v[1]*phiu(in,uy_id)*weight;		//	Tny
                ef(in + first_u,0) += iterm1 + iterm2;
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in + first_p,0)	+= -1.0 * (fCsymetric) * v[2]*phip(in,0)*weight;	// Qnormal
            }
            break;
        }
            
        case 10 :
        {
            //	Neumann condition for x state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                REAL iterm1		= -1.0* v[0]*phiu(in,ux_id)*weight;		//	Tnx
                ef(in + first_u,0) += iterm1;
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in + first_p,0)	+= -1.0  * (fCsymetric) * v[2]*phip(in,0)*weight;	// Qnormal
            }
            break;
        }
            
        case 11 :
        {
            //	Neumann condition for y state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                REAL iterm2		= -1.0* v[1]*phiu(in,uy_id)*weight;		//	Tny
                ef(in + first_u,0) += iterm2;
            }
            
            //	Diffusion Equation
            for(in = 0 ; in < phrp; in++)
            {
                //	Normal Flux on neumman boundary
                ef(in + first_p,0)	+= -1.0 * (fCsymetric) * v[2]*phip(in,0)*weight;	// Qnormal
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
    
    if (fIsMixedQ) {
        this->SolutionMF(datavec, var, Solout);
        return;
    }
    
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
        Solout[1] = -1.0*(fk/feta)*Grad_p(1,0);
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


///////////////////////////////////// monophasic mixed methods ////////////////////////////////////////////////


/** Computes the divergence over the parametric space */
void TPZBiotPoroelasticity::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &Divergence_of_q){

    int qblock = 1;
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

void TPZBiotPoroelasticity::ContributeMF(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    
    int u_b = 0;
    int q_b = 1;
    int p_b = 2;
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[u_b].fShapeType;
    
    // Getting the space functions
    TPZFMatrix<REAL>    &phiu   =   datavec[u_b].phi;
    TPZFMatrix<REAL>    &phiq   =   datavec[q_b].phi;
    TPZFMatrix<REAL>    &phip   =   datavec[p_b].phi;
    
    TPZFMatrix<REAL>    &dphiu   =   datavec[u_b].dphix;
    TPZFMatrix<REAL>    &dphip   =   datavec[p_b].dphix;
    
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,3> u = datavec[u_b].sol[0];
    TPZManVector<REAL,3> q = datavec[q_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix<10,STATE> dq = datavec[q_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    // Transformations
    TPZFNMatrix<27,REAL> grad_phi_u;
    if(!(shapetype == datavec[u_b].EVecShape))
    {
        TPZAxesTools<STATE>::Axes2XYZ(dphiu, grad_phi_u, axes_u);
        
    }
    
    TPZFNMatrix<3,REAL> grad_u;
    TPZAxesTools<STATE>::Axes2XYZ(du, grad_u, axes_u);
    
    
    // Size of matrix blocks
    int nphi_u = phiu.Rows();
    int nphi_q = datavec[q_b].fVecShapeIndex.NElements();
    int nphi_p = phip.Rows();
    
    int first_u = 0;
    int first_q = first_u + 2*nphi_u;
    int first_p = first_q + nphi_q;
    
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
            ef(ip + first_p, 0)		+= (fCsymetric) * weight *  (-1.0) * (1.0/dt) * (falpha * div_u + fSe * p[0]) * phip(ip,0);
        }
        return;
        
    }
    
    for (int iu = 0; iu < nphi_u; iu++) {
        
        if(shapetype == datavec[u_b].EVecShape) // RB functions
        {
            // Computing Gradient of the test function for each component
            Grad_vx_i(0,0) = dphiu(iu,0)*axes_u(0,0)+dphiu(iu,1)*axes_u(1,0); // dvx/dx
            Grad_vx_i(1,0) = dphiu(iu,0)*axes_u(0,1)+dphiu(iu,1)*axes_u(1,1); // dvx/dy
            
            Grad_vy_i(0,0) = dphiu(iu,2)*axes_u(0,0)+dphiu(iu,3)*axes_u(1,0); // dvy/dx
            Grad_vy_i(1,0) = dphiu(iu,2)*axes_u(0,1)+dphiu(iu,3)*axes_u(1,1); // dvy/dy
        }
        else{
            // Computing Gradient of the test function for each component
            for (int d = 0; d < fdimension; d++) {
                Grad_vx_i(d,0) = grad_phi_u(d,iu);
                Grad_vy_i(d,0) = grad_phi_u(d,iu);
            }
        }
        
        ef(2*iu + first_u, 0)   += weight * ((S(0,0) - falpha * p[0]) * Grad_vx_i(0,0) + S(0,1) * Grad_vx_i(1,0));
        ef(2*iu+1 + first_u, 0)	+= weight * (S(1,0) * Grad_vy_i(0,0) + (S(1,1) - falpha * p[0]) * Grad_vy_i(1,0));
        
        for (int ju = 0; ju < nphi_u; ju++) {
            
            if(shapetype == datavec[u_b].EVecShape) // RB functions
            {
                // Computing Gradient of the test function for each component
                Grad_vx_j(0,0) = dphiu(ju,0)*axes_u(0,0)+dphiu(ju,1)*axes_u(1,0); // dvx/dx
                Grad_vx_j(1,0) = dphiu(ju,0)*axes_u(0,1)+dphiu(ju,1)*axes_u(1,1); // dvx/dy
                
                Grad_vy_j(0,0) = dphiu(ju,2)*axes_u(0,0)+dphiu(ju,3)*axes_u(1,0); // dvy/dx
                Grad_vy_j(1,0) = dphiu(ju,2)*axes_u(0,1)+dphiu(ju,3)*axes_u(1,1); // dvy/dy
            }
            else{
                // Computing Gradient of the test function for each component
                for (int d = 0; d < fdimension; d++) {
                    Grad_vx_j(d,0) = grad_phi_u(d,ju);
                    Grad_vy_j(d,0) = grad_phi_u(d,ju);
                }
            }
            
            ek(2*iu + first_u, 2*ju + first_u)      += weight * ( ( (2.0*fmu + flambda) * Grad_vx_j(0,0) ) * Grad_vx_i(0,0) + fmu * Grad_vx_j(1,0) * Grad_vx_i(1,0));
            ek(2*iu + first_u, 2*ju+1 + first_u)    += weight * ( (flambda * Grad_vy_j(1,0) ) * Grad_vx_i(0,0) + fmu * Grad_vy_j(0,0) * Grad_vx_i(1,0)  );
            ek(2*iu+1 + first_u, 2*ju + first_u)	+= weight * ( fmu * Grad_vx_j(1,0) * Grad_vy_i(0,0) + flambda * Grad_vx_j(0,0) * Grad_vy_i(1,0));
            ek(2*iu+1 + first_u, 2*ju+1 + first_u)	+= weight * ( (2.0*fmu + flambda) * Grad_vy_j(1,0) * Grad_vy_i(1,0) + fmu * Grad_vy_j(0,0) * Grad_vy_i(0,0) );
        }
        
        for(int jp = 0; jp < nphi_p; jp++)
        {
            ek(2*iu + first_u,first_p+jp) +=      (-1.0)* weight * falpha * phip(jp,0) * Grad_vx_i(0,0);
            ek(2*iu + 1 + first_u,first_p+jp) +=  (-1.0)* weight * falpha * phip(jp,0) * Grad_vy_i(1,0);
        }
        
    }
    
    // Mixed part
   
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
        
        ef(iq + first_q) += (fMFsymetric) * weight * ( Kl_inv_dot_q - (1.0/jac_det) * p[0] * div_on_master(iq,0));
        
        for (int jq = 0; jq < nphi_q; jq++)
        {
            
            v_j = datavec[q_b].fVecShapeIndex[jq].first;
            s_j = datavec[q_b].fVecShapeIndex[jq].second;
            
            STATE Kl_inv_phi_u_j_dot_phi_u_j = 0.0;
            for (int k = 0; k < q.size(); k++) {
                phi_q_j(k,0) = phiq(s_j,0) * datavec[q_b].fNormalVec(k,v_j);
                Kl_inv_phi_u_j_dot_phi_u_j += (1.0/c)*phi_q_j(k,0)*phi_q_i(k,0);
            }
            
            ek(iq + first_q,jq + first_q) += (fMFsymetric) * weight * Kl_inv_phi_u_j_dot_phi_u_j;
        }
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(iq + first_q, jp + first_p) += (fMFsymetric) * weight * ( - (1.0/jac_det) * phip(jp,0) * div_on_master(iq,0) ) ;
        }
        
    }
    

    REAL div_q = (dq(0,0) + dq(1,1) + dq(2,2));///jac_det;
    // Darcy mixed mono-phasic flow restriction equation
    for (int ip = 0; ip < nphi_p; ip++) {
    
        
        ef(ip + first_p, 0)		+= (fCsymetric) * weight * (div_q + (1.0/dt) * (falpha * div_u + fSe * p[0])) * phip(ip,0);
        
        //	Coupling matrix
        for(int ju = 0; ju < nphi_u; ju++ )
        {
            if(shapetype == datavec[u_b].EVecShape) // RB functions
            {
                // Computing Gradient of the test function for each component
                Grad_vx_j(0,0) = dphiu(ju,0)*axes_u(0,0)+dphiu(ju,1)*axes_u(1,0); // dvx/dx
                Grad_vx_j(1,0) = dphiu(ju,0)*axes_u(0,1)+dphiu(ju,1)*axes_u(1,1); // dvx/dy
                
                Grad_vy_j(0,0) = dphiu(ju,2)*axes_u(0,0)+dphiu(ju,3)*axes_u(1,0); // dvy/dx
                Grad_vy_j(1,0) = dphiu(ju,2)*axes_u(0,1)+dphiu(ju,3)*axes_u(1,1); // dvy/dy
            }
            else{
                // Computing Gradient of the test function for each component
                for (int d = 0; d < fdimension; d++) {
                    Grad_vx_j(d,0) = grad_phi_u(d,ju);
                    Grad_vy_j(d,0) = grad_phi_u(d,ju);
                }
            }
            
            ek(ip + first_p,2*ju + first_u)     += (fCsymetric) * weight * (1.0/dt) * falpha * phip(ip,0) * Grad_vx_j(0,0);
            ek(ip + first_p,2*ju+1 + first_u)   += (fCsymetric) * weight * (1.0/dt) * falpha * phip(ip,0) * Grad_vy_j(1,0);
        }
        
        for (int jq = 0; jq < nphi_q; jq++)
        {
            ek(ip + first_p, jq + first_q) += (fCsymetric) * weight * (1.0/jac_det) * div_on_master(jq,0) * phip(ip,0);
        }
        
        for (int jp = 0; jp < nphi_p; jp++) {

            ek(ip + first_p, jp + first_p)  += (fCsymetric) * weight * ( (1.0/dt) * fSe * phip(jp,0) * phip(ip,0) );
        }
        
    }
}

void TPZBiotPoroelasticity::ContributeMF(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->ContributeMF(datavec, weight, ek_fake, ef);
}

void TPZBiotPoroelasticity::ContributeMFBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int u_b = 0;
    int q_b = 1;
    int p_b = 2;
    
    int ux_id = 0;
    int uy_id = 0;
    TPZMaterialData::MShapeFunctionType shapetype = datavec[u_b].fShapeType;
    REAL c_big = 1.0;
    if(shapetype == datavec[u_b].EVecShape)
    {
        ux_id = 0;
        uy_id = 1;
        c_big = 0.0;
    }
    
    TPZFMatrix<REAL>  &phiu = datavec[u_b].phi;
    TPZFMatrix<REAL>  &phiq = datavec[q_b].phi;
    TPZFMatrix<REAL>  &phip = datavec[p_b].phi;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,3> u = datavec[u_b].sol[0];
    TPZManVector<REAL,3> q = datavec[q_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    int phru = phiu.Rows();
    int phrq = phiq.Rows();
    int phrp = phip.Rows();
    
    int firstu      = 0;
    int firstq      = firstu + 2*phru;
    int firstp      = firstq + phrq;
    
    short in,jn;
    REAL v[3];
    
    v[0] = bc.Val2()(0,0);	//	Ux displacement
    v[1] = bc.Val2()(1,0);	//	Uy displacement
    v[2] = bc.Val2()(2,0);	//	Pressure or flux
    
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
                ef(2*in,0)		+= c_big*gBigNumber*(u[0] - v[0])*phiu(in,ux_id)*weight;	// X displacement Value
                ef(2*in+1,0)	+= c_big*gBigNumber*(u[1] - v[1])*phiu(in,uy_id)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= c_big*gBigNumber*phiu(in,0)*phiu(jn,ux_id)*weight;	// X displacement
                    ek(2*in+1,2*jn+1)	+= c_big*gBigNumber*phiu(in,0)*phiu(jn,uy_id)*weight;	// Y displacement
                }
            }
            
            //	Diffusion Equation Pressure data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += (fMFsymetric) * weight * v[2] * phiq(iq,0);
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
                ef(2*in,0)		+= c_big*gBigNumber*(u[0] - v[0])*phiu(in,ux_id)*weight;	// X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= c_big*gBigNumber*phiu(in,0)*phiu(jn,ux_id)*weight;	// X displacement
                }
            }
            
            //	Diffusion Equation Pressure data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += (fMFsymetric) * weight * v[2] * phiq(iq,0);
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
                ef(2*in+1,0)	+= c_big*gBigNumber*(u[1] - v[1])*phiu(in,uy_id)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)	+= c_big*gBigNumber*phiu(in,0)*phiu(jn,uy_id)*weight;	// Y displacement
                }
            }
            
            
            //	Diffusion Equation Pressure data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += (fMFsymetric) * weight * v[2] * phiq(iq,0);
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
                ef(2*in,0)		+= -1.0* v[0]*phiu(in,ux_id)*weight;		//	Tnx
                ef(2*in+1,0)	+= -1.0* v[1]*phiu(in,uy_id)*weight;		//	Tny
            }
            
            //	Diffusion Equation Pressure data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += (fMFsymetric) * weight * v[2] * phiq(iq,0);
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
                ef(2*in,0)		+= -1.0 * v[0]*phiu(in,ux_id)*weight;		//	Tnx
            }
            
            //	Diffusion Equation Pressure data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += (fMFsymetric) * weight * v[2] * phiq(iq,0);
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
                ef(2*in+1,0)	+= -1.0 * v[1]*phiu(in,uy_id)*weight;		//	Tny
            }
            
            //	Diffusion Equation Pressure data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += (fMFsymetric) * weight * v[2] * phiq(iq,0);
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
                ef(2*in,0)		+= c_big*gBigNumber*(u[0] - v[0])*phiu(in,ux_id)*weight;	// X displacement Value
                ef(2*in+1,0)	+= c_big*gBigNumber*(u[1] - v[1])*phiu(in,uy_id)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= c_big*gBigNumber*phiu(in,0)*phiu(jn,ux_id)*weight;	// X displacement
                    ek(2*in+1,2*jn+1)	+= c_big*gBigNumber*phiu(in,0)*phiu(jn,uy_id)*weight;	// Y displacement
                }
            }
            
            //	Diffusion Equation Neumann data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += weight * gBigNumber * (q[0] - v[2]) * phiq(iq,0);
                
                for (int jq = 0; jq < phrq; jq++)
                {
                    ek(iq + firstq,jq + firstq) += weight * gBigNumber * phiq(jq,0) * phiq(iq,0);
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
                ef(2*in,0)		+= c_big*gBigNumber*(u[0] - v[0])*phiu(in,ux_id)*weight;	// X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= c_big*gBigNumber*phiu(in,0)*phiu(jn,ux_id)*weight;	// X displacement
                }
            }
            
            
            //	Diffusion Equation Neumann data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += weight * gBigNumber * (q[0] - v[2]) * phiq(iq,0);
                
                for (int jq = 0; jq < phrq; jq++)
                {
                    ek(iq + firstq,jq + firstq) += weight * gBigNumber * phiq(jq,0) * phiq(iq,0);
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
                ef(2*in+1,0)	+= c_big*gBigNumber*(u[1] - v[1])*phiu(in,uy_id)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)	+= c_big*gBigNumber*phiu(in,0)*phiu(jn,uy_id)*weight;	// Y displacement
                }
            }
            
            
            //	Diffusion Equation Neumann data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += weight * gBigNumber * (q[0] - v[2]) * phiq(iq,0);
                
                for (int jq = 0; jq < phrq; jq++)
                {
                    ek(iq + firstq,jq + firstq) += weight * gBigNumber * phiq(jq,0) * phiq(iq,0);
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
                ef(2*in,0)		+= -1.0 * v[0]*phiu(in,ux_id)*weight;		//	Tnx
                ef(2*in+1,0)	+= -1.0 * v[1]*phiu(in,uy_id)*weight;		//	Tny
            }
            
            //	Diffusion Equation Neumann data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += weight * gBigNumber * (q[0] - v[2]) * phiq(iq,0);
                
                for (int jq = 0; jq < phrq; jq++)
                {
                    ek(iq + firstq,jq + firstq) += weight * gBigNumber * phiq(jq,0) * phiq(iq,0);
                }
                
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
                ef(2*in,0)		+= -1.0 * v[0]*phiu(in,ux_id)*weight;		//	Tnx
            }
            
            //	Diffusion Equation Neumann data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += weight * gBigNumber * (q[0] - v[2]) * phiq(iq,0);
                
                for (int jq = 0; jq < phrq; jq++)
                {
                    ek(iq + firstq,jq + firstq) += weight * gBigNumber * phiq(jq,0) * phiq(iq,0);
                }
                
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
                ef(2*in+1,0)	+= -1.0 * v[1]*phiu(in,uy_id)*weight;		//	Tny
            }
            
            //	Diffusion Equation Neumann data
            for (int iq = 0; iq < phrq; iq++)
            {
                ef(iq + firstq) += weight * gBigNumber * (q[0] - v[2]) * phiq(iq,0);
                
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
void TPZBiotPoroelasticity::SolutionMF(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    
    int u_b = 0;
    int q_b = 1;
    int p_b = 2;
    
    // Getting the space functions
    
    TPZFNMatrix <9,REAL>	&axes_p	=	datavec[p_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZManVector<REAL,2> q = datavec[q_b].sol[0];
    TPZManVector<REAL,1> p = datavec[p_b].sol[0];
    
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    TPZFNMatrix <6,REAL> dp = datavec[p_b].dsol[0];
    
    
    REAL to_Mpa     = 1.0e-6;
    REAL to_Darcy   = 1.013249966e+15; //md
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_p(3,1,0.0),Grad_u(3,3,0.0),Grad_u_n(3,3,0.0),e_e(3,3,0.0),e_p(3,3,0.0),S;
    
    Grad_p(0,0) = dp(0,0)*axes_p(0,0)+dp(1,0)*axes_p(1,0); // dp/dx
    Grad_p(1,0) = dp(0,0)*axes_p(0,1)+dp(1,0)*axes_p(1,1); // dp/dy
    
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
        for (int d = 0; d <fdimension; d++) {
            Solout[d] = q[d];
        }
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
