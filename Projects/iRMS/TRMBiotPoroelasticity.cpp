//
//  TRMBiotPoroelasticity.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMBiotPoroelasticity.h"


/** @brief Default constructor */
TRMBiotPoroelasticity::TRMBiotPoroelasticity() : TPZMatWithMem<TRMMemory, TPZDiscontinuousGalerkin>()
{
    fdimension = 0;
}

/** @brief Constructor based on a material id */
TRMBiotPoroelasticity::TRMBiotPoroelasticity(int matid, int dimension) : TPZMatWithMem<TRMMemory, TPZDiscontinuousGalerkin>(matid)
{
    fdimension = dimension;
}

/** @brief Default desconstructor */
TRMBiotPoroelasticity::~TRMBiotPoroelasticity(){
    
}

/** @brief Copy constructor $ */
TRMBiotPoroelasticity::TRMBiotPoroelasticity(const TRMBiotPoroelasticity& other){
    this->fdimension    = other.fdimension;
    this->fSimulationData    = other.fSimulationData;
}

/** @brief Copy assignemnt operator $ */
TRMBiotPoroelasticity& TRMBiotPoroelasticity::operator = (const TRMBiotPoroelasticity& other){
    
    if (this != & other) // prevent self-assignment
    {
        this->fdimension    = other.fdimension;
        this->fSimulationData    = other.fSimulationData;
    }
    return *this;
}


void TRMBiotPoroelasticity::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
    
}

void TRMBiotPoroelasticity::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TRMBiotPoroelasticity::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

int TRMBiotPoroelasticity::VariableIndex(const std::string &name) {
    if (!strcmp("u", name.c_str()))   return 0;
    if (!strcmp("s_x", name.c_str())) return 1;
    if (!strcmp("s_y", name.c_str())) return 2;
    if (!strcmp("s_z", name.c_str())) return 3;
    if (!strcmp("s_v", name.c_str())) return 4;
    if (!strcmp("id", name.c_str()))  return 5;
    return TPZMatWithMem::VariableIndex(name);
}

int TRMBiotPoroelasticity::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return fdimension; // vector
        case 1:
            return fdimension; // vector
        case 2:
            return fdimension; // vector
        case 3:
            return fdimension; // vector
        case 4:
            return 1; // scalar
        case 5:
            return 1; // integer
    }
    
    DebugStop();
    return TPZMatWithMem::NSolutionVariables(var);
}

void TRMBiotPoroelasticity::Compute_Sigma(REAL & l, REAL & mu, TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_u){
    
    
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

void TRMBiotPoroelasticity::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int u_b = 0;
    
    if (!fSimulationData->IsCurrentStateQ()) {
        
        return;
    }
    
    // Getting the space functions from memory
    long global_point_index = datavec[u_b].intGlobPtIndex;
    TRMMemory & memory = GetMemory()[global_point_index];
    
    TPZFMatrix<REAL> & phi_u        = memory.phi_u();
    TPZFMatrix<REAL> & grad_phi_u   = memory.grad_phi_u();
    
    TPZFMatrix<REAL> & grad_u   = memory.grad_u();
    TPZFMatrix<REAL> & grad_u_n = memory.grad_u_n();
    
    REAL & p_n  = memory.p_n();
    
    int nphi_u = phi_u.Rows();
    int first_u = 0;
    
    TPZFNMatrix<9,REAL> Grad_vx_i(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_i(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vx_j(fdimension,1,0.0);
    TPZFNMatrix<9,REAL> Grad_vy_j(fdimension,1,0.0);
    
    REAL l_dr   = 2.30769e9;
    REAL mu_dr  = 1.53846e9;
    REAL alpha  = 0.8;
    REAL rho    = 2500.0;
    TPZManVector<REAL,2> b(2,0.0);
    b[0] = 0.0*rho;
    b[1] = -1.0*9.81*rho;
    
    TPZFNMatrix<9,REAL> S(3,3),S_n(3,3);
    Compute_Sigma(l_dr, mu_dr, S, grad_u);
    Compute_Sigma(l_dr, mu_dr, S_n, grad_u_n);
    
    REAL source = alpha * p_n;
    
    for (int iu = 0; iu < nphi_u; iu++) {
        
        // Computing Gradient of the test function for each component
        for (int d = 0; d < fdimension; d++) {
            Grad_vx_i(d,0) = grad_phi_u(d,iu);
            Grad_vy_i(d,0) = grad_phi_u(d,iu);
        }
        
        ef(2*iu + first_u, 0)   += weight * ((S_n(0,0) - source) * Grad_vx_i(0,0) + S_n(0,1) * Grad_vx_i(1,0) - b[0] * phi_u(iu,0));
        ef(2*iu+1 + first_u, 0)	+= weight * (S_n(1,0) * Grad_vy_i(0,0) + (S_n(1,1) - source) * Grad_vy_i(1,0) - b[1] * phi_u(iu,0));
        
        for (int ju = 0; ju < nphi_u; ju++) {
            
            // Computing Gradient of the test function for each component
            for (int d = 0; d < fdimension; d++) {
                Grad_vx_j(d,0) = grad_phi_u(d,ju);
                Grad_vy_j(d,0) = grad_phi_u(d,ju);
            }
            
            ek(2*iu + first_u, 2*ju + first_u)      += weight * ( ( (2.0*mu_dr + l_dr) * Grad_vx_j(0,0) ) * Grad_vx_i(0,0) + mu_dr * Grad_vx_j(1,0) * Grad_vx_i(1,0));
            ek(2*iu + first_u, 2*ju+1 + first_u)    += weight * ( (l_dr * Grad_vy_j(1,0) ) * Grad_vx_i(0,0) + mu_dr * Grad_vy_j(0,0) * Grad_vx_i(1,0)  );
            ek(2*iu+1 + first_u, 2*ju + first_u)	+= weight * ( mu_dr * Grad_vx_j(1,0) * Grad_vy_i(0,0) + l_dr * Grad_vx_j(0,0) * Grad_vy_i(1,0));
            ek(2*iu+1 + first_u, 2*ju+1 + first_u)	+= weight * ( (2.0*mu_dr + l_dr) * Grad_vy_j(1,0) * Grad_vy_i(1,0) + mu_dr * Grad_vy_j(0,0) * Grad_vy_i(0,0) );
        }
        
    }
    
}

void TRMBiotPoroelasticity::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ek_fake, ef);
    
}

void TRMBiotPoroelasticity::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){

    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int u_b = 0;

    // Get the data at the integrations points
    TPZMatWithMem<TRMMemory,TPZBndCond>  & material_mem = dynamic_cast<TPZMatWithMem<TRMMemory,TPZBndCond > & >(bc);
    long global_point_index = datavec[u_b].intGlobPtIndex;
    TRMMemory & memory = material_mem.GetMemory()[global_point_index];
    
    TPZFMatrix<REAL>  &phiu = memory.phi_u();
    TPZFNMatrix<3,REAL> u_m = memory.u_n();
    TPZManVector<REAL,2> u(2,0.0);
    u[0] = u_m(0,0);
    u[1] = u_m(0,1);
    
    int phru = phiu.Rows();
    short in,jn;
    REAL v[2];
    v[0] = bc.Val2()(0,0);	//	Ux displacement
    v[1] = bc.Val2()(1,0);	//	Uy displacement
    
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
            //	Dirichlet condition for y each state variable
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
            for(in = 0 ; in < phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(2*in,0)		+= -1.0 * weight * v[0] * phiu(in,0);		//	Tnx
                ef(2*in+1,0)	+= -1.0 * weight * v[1] * phiu(in,0);		//	Tny
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
            
        case 6 :
        {
            
            //	Neumann condition for each state variable
            //	Elasticity Equation
            
            REAL well_pressure = bc.Val2()(0,0);
            if (bc.HasTimedependentBCForcingFunction()) {
                TPZManVector<STATE,2> f(1);
                TPZFMatrix<double> gradf;
                REAL time = fSimulationData->t();
                bc.TimedependentBCForcingFunction()->Execute(datavec[u_b].x, time, f, gradf);
                well_pressure = f[0];
            }
            

            REAL Tn[2];
            TPZManVector<STATE,3> n =  datavec[u_b].normal;
            
            for (int d = 0; d < Dimension(); d++) {
                Tn[d] = well_pressure * n[d];
            }
            
            for(in = 0 ; in < phru; in++)
            {
                //	Normal Tension Components on neumman boundary
                ef(2*in,0)		+= 1.0 * weight * Tn[0] * phiu(in,0);		//	Tnx
                ef(2*in+1,0)	+= 1.0 * weight * Tn[1] * phiu(in,0);		//	Tny
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


void TRMBiotPoroelasticity::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    DebugStop();
}

void TRMBiotPoroelasticity::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    DebugStop();
}

void TRMBiotPoroelasticity::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    DebugStop();
}

void TRMBiotPoroelasticity::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
    DebugStop();
}

void TRMBiotPoroelasticity::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {

    Solout.Resize( this->NSolutionVariables(var));

    int u_b = 0;
    
    // Getting the solutions and derivatives
    TPZManVector<REAL,2> u = datavec[u_b].sol[0];
    
    REAL l_dr   = 2.30769e9;
    REAL mu_dr  = 1.53846e9;
    
    TPZFNMatrix <6,REAL> du     = datavec[u_b].dsol[0];
    TPZFNMatrix <9,REAL> axes_u	= datavec[u_b].axes;
    
    // Computing Gradient of the Solution
    TPZFNMatrix<9,REAL> Grad_u(3,3,0.0),S(3,3,0.0);
    TPZAxesTools<STATE>::Axes2XYZ(du, Grad_u, axes_u);
    
//    Grad_u(0,0) = du(0,0)*axes_u(0,0)+du(1,0)*axes_u(1,0); // dux/dx
//    Grad_u(0,1) = du(0,0)*axes_u(0,1)+du(1,0)*axes_u(1,1); // dux/dy
//    
//    Grad_u(1,0) = du(0,1)*axes_u(0,0)+du(1,1)*axes_u(1,0); // duy/dx
//    Grad_u(1,1) = du(0,1)*axes_u(0,1)+du(1,1)*axes_u(1,1); // duy/dy
    
    Grad_u.Resize(3, 3);
    this->Compute_Sigma(l_dr, mu_dr, S, Grad_u);
    
    //	Displacements
    if(var == 0){
        for (int d = 0; d < Dimension(); d++) {
            Solout[d] = u[d];
        }
        return;
    }
    
    if(var == 1){
        for (int d = 0; d < Dimension(); d++) {
            Solout[d] = S(0,d);
        }
        return;
    }
    
    if(var == 2){
        for (int d = 0; d < Dimension(); d++) {
            Solout[d] = S(1,d);
        }
        return;
    }
    
//    if(var == 3){
//        for (int d = 0; d < Dimension(); d++) {
//            Solout[d] = S(0,d);
//        }
//        return;
//    }
    
    if(var == 4){
        for (int d = 0; d < Dimension(); d++) {
            Solout[0] += S(d,d);
        }
        Solout[0] *= 1.0/3.0;
        return;
    }
    
    if(var == 5){
        Solout[0] = this->Id();
        return;
    }
    
    std::cout  << " not implemented. " << std::endl;
    DebugStop();
    
}


int TRMBiotPoroelasticity::ClassId() const {
    return -63786378;
}

// -------------------------------------------------------------------------------------------

void TRMBiotPoroelasticity::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
}

// -------------------------------------------------------------------------------------------

void TRMBiotPoroelasticity::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}

// Update element memory by copying the n+1 data to the n data
void TRMBiotPoroelasticity::UpdateMemory()
{
    DebugStop();
}