//
//  TRMPhaseInterfaceTransport.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMPhaseInterfaceTransport.h"

/**
 * Empty Constructor
 */
TRMPhaseInterfaceTransport::TRMPhaseInterfaceTransport() : TPZMatWithMem<TRMPhaseInterfaceMemory, TPZDiscontinuousGalerkin>()
{
    
}

/** Creates a material object and inserts it in the vector of
 *  material pointers of the mesh.
 */
TRMPhaseInterfaceTransport::TRMPhaseInterfaceTransport(int matid) : TPZMatWithMem<TRMPhaseInterfaceMemory, TPZDiscontinuousGalerkin>(matid)
{
    
}


/** Creates a material object based on the referred object and
 *  inserts it in the vector of material pointers of the mesh.
 */
TRMPhaseInterfaceTransport::TRMPhaseInterfaceTransport(const TRMPhaseInterfaceTransport &mat) : TPZMatWithMem<TRMPhaseInterfaceMemory, TPZDiscontinuousGalerkin>(mat)
{
    
}

/**
 * Destructor
 */
TRMPhaseInterfaceTransport::~TRMPhaseInterfaceTransport()
{

}

/** Fill material data parameter with necessary requirements for the
 * Contribute method. Here, in base class, all requirements are considered
 * as necessary. Each derived class may optimize performance by selecting
 * only the necessary data.
 * @since April 10, 2007
 */
void TRMPhaseInterfaceTransport::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    DebugStop();
}

void TRMPhaseInterfaceTransport::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

/** print out the data associated with the material */
void TRMPhaseInterfaceTransport::Print(std::ostream &out)
{
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}


/** returns the variable index associated with the name */
int TRMPhaseInterfaceTransport::VariableIndex(const std::string &name)
{
    DebugStop();
    return 0;
}


/** returns the number of variables associated with the variable
 indexed by var.  var is obtained by calling VariableIndex */
int TRMPhaseInterfaceTransport::NSolutionVariables(int var)
{
    DebugStop();
    return 0;
}


/** Computes the divergence over the parametric space */
void TRMPhaseInterfaceTransport::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU)
{
    DebugStop();
}


/** returns the solution associated with the var index based on
 * the finite element approximation */
void TRMPhaseInterfaceTransport::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout)
{
    DebugStop();
}





// Contribute Methods being used

/**
 * It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}


/**
 * It computes a contribution to the load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}


/**
 * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}


/**
 * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    int sb_a    = 0;
    
    TPZFNMatrix<100,STATE> phi_ss_l       = datavecleft[sb_a].phi;
    
    int nphis_a_l     = phi_ss_l.Rows();
    int firsts_a_l    = 0;
    
    // Get the pressure at the integrations points
    TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond>  * material_bc_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond> *>(bc.Material());
//    long global_point_index = data.intGlobPtIndex;
//    TRMPhaseInterfaceMemory &point_memory = material_bc_mem->GetMemory()[global_point_index];
    REAL p_avg_n    = 0.0;//point_memory.p_avg_n_l();
    REAL sa_avg_n    = 0.0;//point_memory.sa_n_l();
    REAL un_l    = 0.0;//point_memory.un();
    
    TPZManVector<STATE,3> n = data.normal;
    REAL p_l                  = p_avg_n;
    REAL s_l                  = sa_avg_n;
    
    
    //  Average values p_a
    STATE p_a_l    = p_l;
    STATE s_a_l    = s_l;
    
    STATE beta = 0.0;
    
    TPZManVector<STATE, 10> fa_l,v_l(nvars+1);
    
    
    REAL Value_m    = 0.0;
    REAL Value_s    = 0.0;
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<STATE,2> f(2);
        TPZFMatrix<double> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavecleft[sb_a].x, time, f, gradf);
        Value_m = f[0];
        Value_s = f[1];
    }
    else{
        Value_m = bc.Val2()(0,0);
    }
    
    switch (bc.Type()) {
            
        case 0 :    // Dirichlet BC  PD outlet
        {
            
            // upwinding
            if (un_l > 0) {
                beta = 1.0;
            }
            
            STATE p_D = Value_m;
            
            v_l[0] = p_D;
            v_l[1] = s_a_l;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0])*phi_ss_l(is,0)*un_l;
                
                for (int js = 0; js < nphis_a_l; js++) {
                    ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ss_l(js,0) * phi_ss_l(is,0)*un_l;
                }
                
                
            }
            
        }
            break;
            
        case 1 :    // Neumann BC  QN outlet
        {
            
            // upwinding
            if (Value_m > 0) {
                beta = 1.0;
            }
            
            STATE un_N = Value_m;
            
            v_l[0] = p_a_l;
            v_l[1] = s_a_l;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ss_l(is,0)*un_N;
                
                for (int js = 0; js < nphis_a_l; js++) {
                    ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ss_l(js,0) * phi_ss_l(is,0)*un_N;
                }
                
                
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD inlet
        {
            
            // upwinding
            if (un_l > 0) {
                beta = 1.0;
            }
            
            STATE p_D = Value_m;
            
            v_l[0] = p_D;
            v_l[1] = Value_s;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ss_l(is,0)*un_l;
                
            }
            
        }
            break;
            
        case 3 :    // Neumann BC  QN inlet
        {
            
            // upwinding
            if (Value_m < 0) {
                beta = 1.0;
            }
            
            STATE un_N = Value_m;
            
            v_l[0] = p_a_l;
            v_l[1] = Value_s;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ss_l(is,0)*un_N;
                
                
            }
            
        }
            break;
            
        case 4 :    // Neumann BC  Impervious bc
        {
            
            // upwinding
            beta = 1.0;
            
            STATE un_N = 0.0;
            
            v_l[0] = p_a_l;
            v_l[1] = Value_s;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ss_l(is,0)*un_N;
                
                
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


/**
 * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}



/**
 * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    int sb_a    = 0;
    
    TPZFNMatrix<100,STATE> phi_ss_l       = datavecleft[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ss_r       = datavecright[sb_a].phi;
    
    int nphis_a_l     = phi_ss_l.Rows();
    int firsts_a_l    = 0;
    
    int nphis_a_r     = phi_ss_r.Rows();
    int firsts_a_r    = firsts_a_l + nphis_a_l;
    
    TPZManVector<STATE,3> n = data.normal;

    REAL sa_l                  = datavecleft[sb_a].sol[0][0];
    REAL sa_r                  = datavecright[sb_a].sol[0][0];
    REAL p_avg_l = 0.0;
    REAL p_avg_r = 0.0;
    REAL un_l = 0.0;
    
    // Interface memory
    // Get the pressure at the integrations points
    long global_point_index = data.intGlobPtIndex;
    TRMPhaseInterfaceMemory &point_memory = GetMemory()[global_point_index];
    un_l = point_memory.un();
    p_avg_l = point_memory.p_avg_n_l();
    p_avg_r = point_memory.p_avg_n_r();

    
    //  Average values p_a
    
    STATE p_a_l    = p_avg_l;
    STATE s_a_l    = sa_l;
    STATE p_a_r    = p_avg_r;
    STATE s_a_r    = sa_r;
    
    STATE beta = 0.0;
    // upwinding
    if (un_l > 0.0) {
        beta = 1.0;
    }
    
    TPZManVector<STATE, 10> fa_l,v_l(nvars+1),fa_r,v_r(nvars+1);
    v_l[0] = p_a_l;
    v_l[1] = s_a_l;
    v_r[0] = p_a_r;
    v_r[1] = s_a_r;
    
    this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
    this->fSimulationData->PetroPhysics()->fa(fa_r, v_r);

    
    for (int is = 0; is < nphis_a_l; is++) {
        
        ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ss_l(is,0)*un_l;
        
        for (int js = 0; js < nphis_a_l; js++) {
            ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ss_l(js,0) * phi_ss_l(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_a_l, js + firsts_a_r) += +1.0*weight * (1.0-beta) * fa_r[2] * phi_ss_r(js,0) * phi_ss_l(is,0)*un_l;
        }
        
    }
    
    for (int is = 0; is < nphis_a_r; is++) {
        
        ef(is + firsts_a_r) += -1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ss_r(is,0)*un_l;
        
        for (int js = 0; js < nphis_a_l; js++) {
            ek(is + firsts_a_r, js + firsts_a_l) += -1.0*weight * beta * fa_l[2] * phi_ss_l(js,0) * phi_ss_r(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_a_r, js + firsts_a_r) += -1.0*weight * (1.0-beta) * fa_r[2] * phi_ss_r(js,0) * phi_ss_r(is,0)*un_l;
        }
        
    }
}


/**
 * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TRMPhaseInterfaceTransport::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef)
{
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    int sb_a    = 0;
    
    TPZFNMatrix<100,STATE> phi_ss_l       = datavecleft[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ss_r       = datavecright[sb_a].phi;
    
    int nphis_a_l     = phi_ss_l.Rows();
    int firsts_a_l    = 0;
    
    int nphis_a_r     = phi_ss_r.Rows();
    int firsts_a_r    = firsts_a_l + nphis_a_l;
    
    TPZManVector<STATE,3> n = data.normal;
    
    REAL sa_l                  = datavecleft[sb_a].sol[0][0];
    REAL sa_r                  = datavecright[sb_a].sol[0][0];
    REAL p_avg_l = 0.0;
    REAL p_avg_r = 0.0;
    REAL un_l = 0.0;
    
    // Interface memory
    // Get the pressure at the integrations points
    long global_point_index = data.intGlobPtIndex;
    TRMPhaseInterfaceMemory &point_memory = GetMemory()[global_point_index];
    un_l = point_memory.un();
    p_avg_l = point_memory.p_avg_n_l();
    p_avg_r = point_memory.p_avg_n_r();
    
    
    //  Average values p_a
    
    STATE p_a_l    = p_avg_l;
    STATE s_a_l    = sa_l;
    STATE p_a_r    = p_avg_r;
    STATE s_a_r    = sa_r;
    
    STATE beta = 0.0;
    // upwinding
    if (un_l > 0.0) {
        beta = 1.0;
    }
    
    TPZManVector<STATE, 10> fa_l,v_l(nvars+1),fa_r,v_r(nvars+1);
    v_l[0] = p_a_l;
    v_l[1] = s_a_l;
    v_r[0] = p_a_r;
    v_r[1] = s_a_r;
    
    this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
    this->fSimulationData->PetroPhysics()->fa(fa_r, v_r);
    
    for (int is = 0; is < nphis_a_l; is++) {
        
        ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ss_l(is,0)*un_l;
        
    }
    
    for (int is = 0; is < nphis_a_r; is++) {
        
        ef(is + firsts_a_r) += -1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ss_r(is,0)*un_l;
        
    }

}


/**
 * Unique identifier for serialization purposes
 */
int TRMPhaseInterfaceTransport::ClassId() const
{
    DebugStop();
}


/**
 * Save the element data to a stream
 */
void TRMPhaseInterfaceTransport::Write(TPZStream &buf, int withclassid)
{
    DebugStop();
}


/**
 * Read the element data from a stream
 */
void TRMPhaseInterfaceTransport::Read(TPZStream &buf, void *context)
{
    DebugStop();
}



/// Copy the n+1 data to the n data
void TRMPhaseInterfaceTransport::UpdateMemory()
{
    DebugStop();
}

