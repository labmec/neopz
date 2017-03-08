//
//  TPZElasticBiot.cpp
//  PZ
//
//  Created by Omar on 3/5/17.
//
//

#include "TPZElasticBiot.h"


/** @brief Default constructor */
TPZElasticBiot::TPZElasticBiot() : TPZMatWithMem<TPZElasticBiotMemory, TPZDiscontinuousGalerkin>(){
    
    fdimension = 0;
    
}

/** @brief Constructor based on a material id */
TPZElasticBiot::TPZElasticBiot(int matid, int dimension) : TPZMatWithMem<TPZElasticBiotMemory, TPZDiscontinuousGalerkin>(matid){
    
    fdimension = dimension;
}

/** @brief Constructor based on a Biot Poroelasticity  object */
TPZElasticBiot::TPZElasticBiot(const TPZElasticBiot &mat) : TPZMatWithMem<TPZElasticBiotMemory, TPZDiscontinuousGalerkin>(mat){
    
    this->fdimension    = mat.fdimension;
}


/** @brief Default destructor */
TPZElasticBiot::~TPZElasticBiot(){
    
}

/** @brief Set the required data at each integration point */
void TPZElasticBiot::FillDataRequirements(TPZVec<TPZMaterialData> &datavec){
    
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
        DebugStop();
        return;
    }
    
}

/** @brief Set the required data at each integration point */
void TPZElasticBiot::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec){
    
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
        DebugStop();
        return;
    }
    
}


/** print out the data associated with the material */
void TPZElasticBiot::Print(std::ostream &out){
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

/** returns the variable index associated with the name */
int TPZElasticBiot::VariableIndex(const std::string &name){
    
    //	Elasticity Variables
    if(!strcmp("u",name.c_str()))              return	0;
    
    //	Elasticity Variables
    if(!strcmp("ue",name.c_str()))              return	1;
    
    return -1;
    
}

/** returns the number of variables associated with the variable
 indexed by var.  var is obtained by calling VariableIndex */
int TPZElasticBiot::NSolutionVariables(int var){
    
    switch(var) {
        case 0:
            return fdimension; // Vector
        case 1:
            return fdimension; // Vector
    }
    return TPZMatWithMem::NSolutionVariables(var);
    
}


void TPZElasticBiot::Compute_Sigma(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_u){
    
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
void TPZElasticBiot::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){

    int u_b = 0;
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[u_b].fShapeType;
    if(shapetype == datavec[u_b].EVecShape){
        this->ContributeRB(datavec, weight, ek, ef);
        return;
    }
    
//    // Getting the space functions from memory
//    long global_point_index = datavec[u_b].intGlobPtIndex;
//    TPZElasticBiotMemory &point_memory = GetMemory()[global_point_index];
//    
//    TPZFMatrix<REAL> & phi_u = point_memory.phi_u();
//    TPZFMatrix<REAL> & grad_phi_u = point_memory.grad_phi_u();
//
//    TPZFMatrix<REAL> & u_n = point_memory.u_n();
//    TPZFMatrix<REAL> & grad_u_n = point_memory.grad_u_n();
    
    
    // Getting the space functions
    TPZFMatrix<REAL>    &phiu   =   datavec[u_b].phi;
    TPZFMatrix<REAL>    &dphiu   =   datavec[u_b].dphix;
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZFNMatrix <15,REAL> du = datavec[u_b].dsol[0];
    
    // Transformations
    TPZFNMatrix<27,REAL> grad_phi_u;
    TPZFNMatrix<3,REAL> grad_u;
    
    TPZAxesTools<STATE>::Axes2XYZ(dphiu, grad_phi_u, axes_u);
    TPZAxesTools<STATE>::Axes2XYZ(du, grad_u, axes_u);
    
//    grad_u_n.Print("grad_u_n omar = ");
//    grad_u.Print("grad_u pz = ");
    
//    u_n.Print("u_n = ");
//    std::cout << "u = " <<  u << std::endl;
//    std::cout << std::endl;
    
//    grad_u_n.Print("du_n = ");
//    du.Print("du  = ");
//    std::cout << std::endl;
    
//    grad_u_n.Print("grad_u_n omar = ");
//    grad_u.Print("grad_u pz = ");
    
    REAL p = 0.0;
    
    int nphi_u = phiu.Rows();
    int first_u = 0;
    
    REAL div_u = grad_u(0,0) + grad_u(1,1);
    
    TPZFNMatrix<6,REAL> Grad_u(3,3,0.0);
    TPZFNMatrix<9,REAL> S(3,3,0.0);
    grad_u.Resize(3, 3);
    this->Compute_Sigma(S,grad_u);
    
    TPZFNMatrix<9,REAL> Grad_vx_i(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_i(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vx_j(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_j(fdimension,1,0.0);
    
    if (!fSimulationData->IsCurrentStateQ()) {
        
        return;
    }
    
    for (int iu = 0; iu < nphi_u; iu++) {
        
        // Computing Gradient of the test function for each component
        for (int d = 0; d < fdimension; d++) {
            Grad_vx_i(d,0) = grad_phi_u(d,iu);
            Grad_vy_i(d,0) = grad_phi_u(d,iu);
        }
        
        ef(2*iu + first_u, 0)   += weight * ((S(0,0) - falpha * p) * Grad_vx_i(0,0) + S(0,1) * Grad_vx_i(1,0));
        ef(2*iu+1 + first_u, 0)	+= weight * (S(1,0) * Grad_vy_i(0,0) + (S(1,1) - falpha * p) * Grad_vy_i(1,0));
        
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
        
    }

}


void TPZElasticBiot::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
}

void TPZElasticBiot::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    int u_b = 0;
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[u_b].fShapeType;
    REAL c_big = 1.0;
    if(shapetype == datavec[u_b].EVecShape)
    {
        this->ContributeRB_BC(datavec, weight, ek, ef, bc);
        return;
    }
    
    TPZFMatrix<REAL>  &phiu = datavec[u_b].phi;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    
    int phru = phiu.Rows();
    short in,jn;
    REAL v[2];
    v[0] = bc.Val2()(0,0);	//	Ux displacement
    v[1] = bc.Val2()(1,0);	//	Uy displacement
    
    REAL time = this->SimulationData()->t();
    REAL dt  = this->SimulationData()->dt();
    REAL Value = bc.Val2()(0,0);
    if (bc.HasTimedependentBCForcingFunction()) {
        TPZManVector<REAL,3> f(2);
        TPZFMatrix<REAL> gradf;
        bc.TimedependentBCForcingFunction()->Execute(datavec[u_b].x, time, f, gradf);
        v[0] = f[0];	//	Ux displacement or Tx
        v[1] = f[1];	//	Uy displacement or Ty
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
                ef(2*in,0)		+= c_big*gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                ef(2*in+1,0)	+= c_big*gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= c_big*gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(2*in+1,2*jn+1)	+= c_big*gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
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
                ef(2*in,0)		+= c_big*gBigNumber*(u[0] - v[0])*phiu(in,0)*weight;	// X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= c_big*gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                }
            }
            break;
            
        }
            
        case 2 :
        {
            //	Dirichlet condition for y each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in+1,0)	+= c_big*gBigNumber*(u[1] - v[1])*phiu(in,0)*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)	+= c_big*gBigNumber*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
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
            
            break;
        }
            
        case 4 :
        {
            //	Neumann condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(2*in,0)		+= -1.0* v[0]*phiu(in,0)*weight;		//	Tnx
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
                ef(2*in+1,0)	+= -1.0* v[1]*phiu(in,0)*weight;		//	Tny
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
void TPZElasticBiot::ContributeRB(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int u_b = 0;
    
    // RB functions must to include axes transformations .
    
    // Getting RB functions and solution form integration points
    // Get the data at the integrations points
    
    long global_point_index = datavec[u_b].intGlobPtIndex;
    TPZElasticBiotMemory &point_memory = GetMemory()[global_point_index];
    TPZFMatrix<REAL> & phiu = point_memory.phi_u();
    TPZFMatrix<REAL> & dphiu = point_memory.grad_phi_u();
    
    TPZFNMatrix<3,REAL> & int_u = point_memory.u();
    TPZFNMatrix<3,REAL> & int_u_n = point_memory.u_n();
    
    TPZFNMatrix<9,REAL> & int_grad_u = point_memory.grad_u();
    TPZFNMatrix<9,REAL> & int_grad_u_n = point_memory.grad_u_n();
    
    // Getting the space functions
    //    TPZFMatrix<REAL>    &phiu   =   datavec[u_b].phi;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,10> & u = datavec[u_b].sol[0];
    TPZFNMatrix <15,REAL> & du = datavec[u_b].dsol[0];
    
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
    
    REAL p = 0.0;
    
    int nphi_u = phiu.Rows();
    
    int first_u = 0;
    
    REAL dt = fSimulationData->dt();
    REAL div_u = grad_u(0,0) + grad_u(1,1);
    
    TPZFNMatrix<9,REAL> S(3,3,0.0);
    Compute_Sigma(S,grad_u);
    
    TPZFNMatrix<9,REAL> Grad_vx_i(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_i(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vx_j(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_j(fdimension,1,0.0);
    
    if (!fSimulationData->IsCurrentStateQ()) {
        
        return;
    }
    
    
    for (int iu = 0; iu < nphi_u; iu++) {
        
        // Computing Gradient of the test function for each component
        Grad_vx_i(0,0) = dphiu(iu,0); // dvx/dx
        Grad_vx_i(1,0) = dphiu(iu,1); // dvx/dy
        
        Grad_vy_i(0,0) = dphiu(iu,2); // dvy/dx
        Grad_vy_i(1,0) = dphiu(iu,3); // dvy/dy
        
        REAL iterm1 = weight * ((S(0,0) - falpha * p) * Grad_vx_i(0,0) + S(0,1) * Grad_vx_i(1,0));
        REAL iterm2 = weight * (S(1,0) * Grad_vy_i(0,0) + (S(1,1) - falpha * p) * Grad_vy_i(1,0));
        ef(iu + first_u, 0) += iterm1 + iterm2;
        
        for (int ju = 0; ju < nphi_u; ju++) {
            
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
        
    }
    
}

void TPZElasticBiot::ContributeRB_BC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    int u_b = 0;
    
    int ux_id = 0;
    int uy_id = 1;
    REAL c_big = 1.0;
    
    // Getting RB functions and solution form integration points
    // Get the data at the integrations points
    TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond>  & material_bc_mem = dynamic_cast<TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond > & >(bc);
    
    TPZFNMatrix<3,REAL>  int_u, int_u_n;
    TPZFNMatrix<9,REAL>  int_grad_u, int_grad_u_n;
    long global_point_index = datavec[u_b].intGlobPtIndex;
    TPZElasticBiotMemory &point_memory = material_bc_mem.GetMemory()[global_point_index];
    TPZFMatrix<REAL> & phiu = point_memory.phi_u();
    
    int_u = point_memory.u();
    int_u_n = point_memory.u_n();
    
    int_grad_u = point_memory.grad_u();
    int_grad_u_n = point_memory.grad_u_n();
    
    
    //    TPZFMatrix<REAL>  &phiu = datavec[u_b].phi;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    
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
    
    int first_u = 0;
    
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
        bc.TimedependentBCForcingFunction()->Execute(datavec[u_b].x, time, f, gradf);
        v[0] = f[0];	//	Ux displacement or Tx
        v[1] = f[1];	//	Uy displacement or Ty
        v[2] = f[2];	//	Pressure
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
            break;
        }
            
        case 1 :
        {
            //	Neumann condition for y state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                REAL iterm2		= -1.0* v[1]*phiu(in,uy_id)*weight;		//	Tny
                ef(in + first_u,0) += iterm2;
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
void TPZElasticBiot::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    
    int u_b = 0;
    
    // Getting the space functions
    TPZFNMatrix <9,REAL>	&axes_u	=	datavec[u_b].axes;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    TPZFNMatrix <6,REAL> du = datavec[u_b].dsol[0];
    
    
    REAL to_Mpa     = 1.0e-6;
    
    // Computing Gradient of the Solution
    TPZFNMatrix<6,REAL> Grad_u(3,3,0.0),Grad_u_n(3,3,0.0),e_e(3,3,0.0),e_p(3,3,0.0),S;
    
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
    
    //	darcy
    if(var == 1) {
        
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
    
    std::cout  << "not implemented. " << std::endl;
    DebugStop();
    
}


/** @brief Unique identifier for serialization purposes */
int TPZElasticBiot::ClassId() const {
    return -6378637866;
}

/** @brief Save object data to a stream */
void TPZElasticBiot::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
}

/** @brief Read object data from a stream */
void TPZElasticBiot::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}