//
//  TRMPhaseTransport.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMPhaseTransport.h"

TRMPhaseTransport::TRMPhaseTransport() : TPZMatWithMem<TRMPhaseMemory, TPZDiscontinuousGalerkin>()
{
    
    /** @brief define the simulation data */
    fSimulationData = NULL;
    
    /** @brief define the transfer matrices */
    fTransfer = NULL;
    
}

TRMPhaseTransport::TRMPhaseTransport(int matid) : TPZMatWithMem<TRMPhaseMemory, TPZDiscontinuousGalerkin>(matid)
{
    
    /** @brief define the simulation data */
    fSimulationData = NULL;
    
    /** @brief define the transfer matrices */
    fTransfer = NULL;
    
}


TRMPhaseTransport::TRMPhaseTransport(const TRMPhaseTransport &mat) : TPZMatWithMem<TRMPhaseMemory, TPZDiscontinuousGalerkin>(mat)
{
    
}

TRMPhaseTransport::~TRMPhaseTransport()
{
    
}

void TRMPhaseTransport::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

void TRMPhaseTransport::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TRMPhaseTransport::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

int TRMPhaseTransport::VariableIndex(const std::string &name) {
    if (!strcmp("p", name.c_str())) return 0;
    if (!strcmp("u", name.c_str())) return 1;
    if (!strcmp("div_u", name.c_str())) return 2;
    if (!strcmp("s_a", name.c_str())) return 3;

    return TPZMatWithMem::VariableIndex(name);
}

int TRMPhaseTransport::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return 3; // Vector
        case 2:
            return 1; // Scalar
        case 3:
            return 1; // Scalar
    }
    return TPZMatWithMem::NSolutionVariables(var);
}

void TRMPhaseTransport::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int sb_a    = 0;
    REAL sa = datavec[sb_a].sol[0][0];

    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0:
        {
            DebugStop();
        }
            break;
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            DebugStop();
        }
            break;
        case 3:
        {
            Solout[0] = sa;
        }
            break;
        default:
        {
            DebugStop();
        }
    }
    
}

// Jacobian contribution
void TRMPhaseTransport::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    int nvars = 4; // {p,sa,sb,t}
    
    int sb_a    = 0;
    
    TPZFNMatrix<100,STATE> phi_ss       = datavec[sb_a].phi;
    int nphis_a     = phi_ss.Rows();
    int firsts_a    = 0;
    
    REAL s                  = datavec[sb_a].sol[0][0];
    
    // Time
    STATE dt = fSimulationData->dt();
    
    // Get the pressure at the integrations points
    long global_point_index = datavec[sb_a].intGlobPtIndex;
    TRMPhaseMemory &point_memory = GetMemory()[global_point_index];
    REAL p_avg_n    = point_memory.p_avg_n();
    REAL p_avg      = point_memory.p_avg();
    
    REAL sa_avg_n    = point_memory.sa_n();
    REAL sa_avg      = point_memory.sa_n();
    
    //  Average values p_a
    
    REAL p_a    = p_avg_n;
    REAL s_a    = sa_avg_n;
    
    //  Computing closure relationship at given average values
    
    TPZManVector<STATE, 10> v(nvars);
    v[0] = p_a;
    v[1] = s_a;
    
    // Fluid parameters
    TPZManVector<STATE, 10> rho_a,rho_b,l;
    fSimulationData->AlphaProp()->Density(rho_a, v);
    fSimulationData->BetaProp()->Density(rho_b, v);
    fSimulationData->PetroPhysics()->l(l, v);
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    TPZManVector<STATE, 10> phi;
    fSimulationData->Map()->Kappa(datavec[sb_a].x, K, Kinv, v);
    fSimulationData->Map()->phi(datavec[sb_a].x, phi, v);
    
    // Defining local variables
    TPZFNMatrix<3,STATE> lambda_K_inv_u(3,1),lambda_dp_K_inv_u(3,1), lambda_ds_K_inv_u(3,1), lambda_K_inv_phi_u_j(3,1);
    TPZManVector<STATE,3> Gravity = fSimulationData->Gravity();
    
    // Integration point contribution

    
    if(! fSimulationData->IsCurrentStateQ()){

        
        for (int is = 0; is < nphis_a; is++)
        {
            
            ef(is + firsts_a) += weight * (-1.0/dt) * s * rho_a[0] * phi[0] * phi_ss(is,0);
            
        }
        
        return;
    }
    
    
    for (int is = 0; is < nphis_a; is++)
    {
        
        ef(is + firsts_a) += weight * (1.0/dt) * s * rho_a[0] * phi[0] * phi_ss(is,0);
        
        for (int js = 0; js < nphis_a; js++)
        {
            ek(is + firsts_a, js + firsts_a) += weight * (1.0/dt) * rho_a[0] * phi[0] * phi_ss(js,0) * phi_ss(is,0);
        }
        
    }
    
    if(fSimulationData->IsThreePhaseQ()){
        DebugStop();
    }
    
}


void TRMPhaseTransport::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    int nvars = 4; // {p,sa,sb,t}
    
    int sb_a    = 0;
    
    TPZFNMatrix<100,STATE> phi_ss       = datavec[sb_a].phi;
    int nphis_a     = phi_ss.Rows();
    int firsts_a    = 0;
    
    REAL s                  = datavec[sb_a].sol[0][0];
    
    // Time
    STATE dt = fSimulationData->dt();
    
    // Get the pressure at the integrations points
    long global_point_index = datavec[sb_a].intGlobPtIndex;
    TRMPhaseMemory &point_memory = GetMemory()[global_point_index];
    REAL p_avg_n    = point_memory.p_avg_n();
    REAL p_avg      = point_memory.p_avg();
    
    REAL sa_avg_n    = point_memory.sa_n();
    REAL sa_avg      = point_memory.sa_n();
    
    //  Average values p_a
    
    REAL p_a    = p_avg_n;
    REAL s_a    = sa_avg_n;
    
    //  Computing closure relationship at given average values
    
    TPZManVector<STATE, 10> v(nvars);
    v[0] = p_a;
    v[1] = s_a;
    
    // Fluid parameters
    TPZManVector<STATE, 10> rho_a,rho_b,l;
    fSimulationData->AlphaProp()->Density(rho_a, v);
    fSimulationData->BetaProp()->Density(rho_b, v);
    fSimulationData->PetroPhysics()->l(l, v);
    
    // Rock parameters
    TPZFNMatrix<9,STATE> K,Kinv;
    TPZManVector<STATE, 10> phi;
    fSimulationData->Map()->Kappa(datavec[sb_a].x, K, Kinv, v);
    fSimulationData->Map()->phi(datavec[sb_a].x, phi, v);
    
    // Defining local variables
    TPZFNMatrix<3,STATE> lambda_K_inv_u(3,1),lambda_dp_K_inv_u(3,1), lambda_ds_K_inv_u(3,1), lambda_K_inv_phi_u_j(3,1);
    TPZManVector<STATE,3> Gravity = fSimulationData->Gravity();
    
    // Integration point contribution
    
    
    if(! fSimulationData->IsCurrentStateQ()){
        
        
        for (int is = 0; is < nphis_a; is++)
        {
            
            ef(is + firsts_a) += weight * (-1.0/dt) * s * rho_a[0] * phi[0] * phi_ss(is,0);
            
        }
        
        return;
    }
    
    
    for (int is = 0; is < nphis_a; is++)
    {
        
        ef(is + firsts_a) += weight * (1.0/dt) * s * rho_a[0] * phi[0] * phi_ss(is,0);
        
    }
    
    if(fSimulationData->IsThreePhaseQ()){
        DebugStop();
    }
    
}


int TRMPhaseTransport::ClassId() const {
    return -637806378;
}

// -------------------------------------------------------------------------------------------

void TRMPhaseTransport::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
}

// -------------------------------------------------------------------------------------------

void TRMPhaseTransport::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}

// Update element memory by copying the n+1 data to the n data
void TRMPhaseTransport::UpdateMemory()
{
    DebugStop();
//    long nel = fMemory.NElements();
//    for (long el=0; el<nel; el++) {
//        fMemory[el].UpdateSolutionMemory();
//    }
}

