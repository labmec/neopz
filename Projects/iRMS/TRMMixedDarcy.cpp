//
//  TRMMixedDarcy.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMMixedDarcy.h"


TRMMixedDarcy::TRMMixedDarcy() : TPZMatWithMem<TRMMemory, TPZDiscontinuousGalerkin>()
{

}

TRMMixedDarcy::TRMMixedDarcy(int matid) : TPZMatWithMem<TRMMemory, TPZDiscontinuousGalerkin>(matid)
{

}


TRMMixedDarcy::TRMMixedDarcy(const TRMMixedDarcy &mat) : TPZMatWithMem<TRMMemory, TPZDiscontinuousGalerkin>(mat)
{

}

TRMMixedDarcy::~TRMMixedDarcy()
{
    
}

void TRMMixedDarcy::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

void TRMMixedDarcy::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TRMMixedDarcy::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

int TRMMixedDarcy::VariableIndex(const std::string &name) {
    if (!strcmp("p", name.c_str())) return 0;
    if (!strcmp("u", name.c_str())) return 1;
    if (!strcmp("div_u", name.c_str())) return 2;
    if (!strcmp("AWeightedPressure", name.c_str())) return 3;
    if (!strcmp("ABulkVelocity", name.c_str())) return 4;
    if (!strcmp("ADivOfBulkVeclocity", name.c_str())) return 5;
    return TPZMatWithMem::VariableIndex(name);
}

int TRMMixedDarcy::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return 3; // Vector
        case 2:
            return 1; // Scalar
        case 3:
            return 1; // Scalar
        case 4:
            return 3; // Vector
        case 5:
            return 1; // Scalar
    }
    return TPZMatWithMem::NSolutionVariables(var);
}

void TRMMixedDarcy::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int Qblock = 0;
    int Pblock = 1;
    
    TPZManVector<REAL,3> Q = datavec[Qblock].sol[0];
    REAL P = datavec[Pblock].sol[0][0];
    
    TPZFMatrix<STATE> dQdx = datavec[Qblock].dsol[0];
    TPZFMatrix<STATE> dPdx = datavec[Pblock].dsol[0];

    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0:
        {
            Solout[0] = P;
        }
            break;
        case 1:
        {
            Solout[0] = Q[0]; // Bulk mass velocity
            Solout[1] = Q[1]; // Bulk mass velocity
            Solout[2] = Q[2]; // Bulk mass velocity
        }
            break;
        case 2:
        {
            Solout[0] = dQdx(0,0) + dQdx(1,1) + dQdx(2,2);
        }
            break;
        case 3:
        {
            REAL x = datavec[Qblock].x[0];
            REAL y = datavec[Qblock].x[1];
            REAL z = datavec[Qblock].x[2];
            Solout[0] = (1. - x)*x + (1. - y)*y + (1. - z)*z;
        }
            break;
        case 4:
        {
            REAL x = datavec[Qblock].x[0];
            REAL y = datavec[Qblock].x[1];
            REAL z = datavec[Qblock].x[2];
            Solout[0] = -1. + 2*x;//2.*(-0.5 + x)*(-1. + y)*y*(-1. + z)*z; // Bulk mass velocity
            Solout[1] = -1. + 2*y;//2.*(-1. + x)*x*(-0.5 + y)*(-1. + z)*z; // Bulk mass velocity
            Solout[2] = -1. + 2*z;//2.*(-1. + x)*x*(-1. + y)*y*(-0.5 + z); // Bulk mass velocity
        }
            break;
        case 5:
        {
//            REAL x = datavec[Qblock].x[0];
//            REAL y = datavec[Qblock].x[1];
//            REAL z = datavec[Qblock].x[2];
            Solout[0] = 6;//2.*(-1. + x)*x*(-1. + y)*y + 2.*(-1. + x)*x*(-1. + z)*z + 2.*(-1. + y)*y*(-1. + z)*z;
        }
            break;
        default:
        {
            TPZMatWithMem::Solution(datavec, var, Solout);
        }
    }
}

// Divergence on master element
void TRMMixedDarcy::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU)
{
    int ublock = 0;
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;   // For H1  test functions Q
    TPZFMatrix<STATE> dphiuH1       = datavec[ublock].dphi; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiuH1axes   = datavec[ublock].dphix; // Derivative For H1  test functions
    TPZFNMatrix<9,STATE> gradu = datavec[ublock].dsol[0];
    TPZFNMatrix<9,STATE> graduMaster;
    gradu.Transpose();
    
    TPZFNMatrix<660> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[ublock].axes);
    
    int nphiuHdiv = datavec[ublock].fVecShapeIndex.NElements();
    
    DivergenceofPhi.Resize(nphiuHdiv,1);
    
    REAL JacobianDet = datavec[ublock].detjac;
    
    TPZFMatrix<STATE> Qaxes = datavec[ublock].axes;
    TPZFMatrix<STATE> QaxesT;
    TPZFMatrix<STATE> Jacobian = datavec[ublock].jacobian;
    TPZFMatrix<STATE> JacobianInverse = datavec[ublock].jacinv;
    
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
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            VectorOnXYZ(0,0) = datavec[ublock].fNormalVec(0,ivectorindex);
            VectorOnXYZ(1,0) = datavec[ublock].fNormalVec(1,ivectorindex);
            VectorOnXYZ(2,0) = datavec[ublock].fNormalVec(2,ivectorindex);
            
            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
            VectorOnMaster *= JacobianDet;
            
            /* Contravariant Piola mapping preserves the divergence */
            // the division by the jacobianDet is to make the integral on the master element???
            DivergenceofPhi(iq,0) = ( dphiuH1(0,ishapeindex)*VectorOnMaster(0,0) +
                                                          dphiuH1(1,ishapeindex)*VectorOnMaster(1,0) +
                                                          dphiuH1(2,ishapeindex)*VectorOnMaster(2,0) );
        }
        
        GradOfXInverse.Multiply(gradu, graduMaster);
        graduMaster *= JacobianDet;
        DivergenceofU = (graduMaster(0,0)+graduMaster(1,1)+graduMaster(2,2));
    }
    else
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            /* Computing the divergence for constant jacobian elements */
            DivergenceofPhi(iq,0) =  datavec[ublock].fNormalVec(0,ivectorindex)*GradphiuH1(0,ishapeindex) +
            datavec[ublock].fNormalVec(1,ivectorindex)*GradphiuH1(1,ishapeindex) +
            datavec[ublock].fNormalVec(2,ivectorindex)*GradphiuH1(2,ishapeindex) ;
        }
    }
    
    return;
    
}

// Jacobian contribution
void TRMMixedDarcy::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2         = datavec[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<STATE> dphiuH1   = datavec[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2   = datavec[Pblock].dphix; // Derivative For L2  test functions
    
    TPZFMatrix<STATE> DivergenceOnMaster;
    STATE divFlux;
    // Compute the divergence on deformed element by piola contravariant transformation
    this->ComputeDivergenceOnMaster(datavec, DivergenceOnMaster,divFlux);
    
    REAL JacobianDet = datavec[ublock].detjac;

    // Blocks dimensions and lengths
    int nphiuHdiv   = datavec[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2     = phiPL2.Rows();                                    // For L2   P
    int iniu    = 0;
    int iniP    = nphiuHdiv     + iniu;
    
    
    // Getting linear combinations from different approximation spaces
    TPZManVector<REAL,3> u      = datavec[ublock].sol[0];
//    REAL                 P      = datavec[Pblock].sol[0][0];
    
    TPZFMatrix<STATE> Graduaxes = datavec[ublock].dsol[0]; // Piola divengence may works, needed set piola computation on the solution elchiv method!!!
    TPZFMatrix<STATE> GradPaxes = datavec[Pblock].dsol[0];
    TPZFNMatrix<660> GradP;
    TPZAxesTools<REAL>::Axes2XYZ(GradPaxes, GradP, datavec[Pblock].axes);
    TPZFMatrix<STATE> KInverse(3,3,0.0);
    
    KInverse(0,0) = 1.0;
    KInverse(1,1) = 1.0;
    KInverse(2,2) = 1.0;
    
    
    
    // Defining local variables
    TPZFMatrix<STATE> oneoverlambda_Kinv_u(3,1);
    TPZFMatrix<STATE> oneoverlambda_Kinv_jphiuHdiv(3,1);
    
    oneoverlambda_Kinv_u(0,0) = (1.0/1.0)* (KInverse(0,0)*u[0] + KInverse(0,1)*u[1] + KInverse(0,2)*u[2]);
    oneoverlambda_Kinv_u(1,0) = (1.0/1.0)* (KInverse(1,0)*u[0] + KInverse(1,1)*u[1] + KInverse(1,2)*u[2]);
    oneoverlambda_Kinv_u(2,0) = (1.0/1.0)* (KInverse(2,0)*u[0] + KInverse(2,1)*u[1] + KInverse(2,2)*u[2]);
    
    
    TPZFMatrix<STATE> iphiuHdiv(3,1);
    TPZFMatrix<STATE> jphiuHdiv(3,1);
    int ishapeindex;
    int ivectorindex;
    int jshapeindex;
    int jvectorindex;
    
    for (int iq = 0; iq < nphiuHdiv; iq++)
    {
        
        /* $ \underset{\Omega_{e}}{\int}\left(K\lambda\right)^{-1}\mathbf{q}\cdot\mathbf{v}\;\partial\Omega_{e}-\underset{\Omega_{e}}{\int}P\; div\left(\mathbf{v}\right)\partial\Omega-\underset{\Omega_{e}}{\int}\nabla\left(\rho_{f}g\; z\right)\cdot\mathbf{v} $ */
        
        ivectorindex    = datavec[ublock].fVecShapeIndex[iq].first;
        ishapeindex     = datavec[ublock].fVecShapeIndex[iq].second;
        
        iphiuHdiv(0,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(0,ivectorindex);
        iphiuHdiv(1,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(1,ivectorindex);
        iphiuHdiv(2,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(2,ivectorindex);
        
        // du/dalphau terms
        for (int jq = 0; jq < nphiuHdiv; jq++)
        {
            jvectorindex = datavec[ublock].fVecShapeIndex[jq].first;
            jshapeindex = datavec[ublock].fVecShapeIndex[jq].second;
            
            jphiuHdiv(0,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(0,jvectorindex);
            jphiuHdiv(1,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(1,jvectorindex);
            jphiuHdiv(2,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(2,jvectorindex);
            
            oneoverlambda_Kinv_jphiuHdiv(0,0) = (1.0/1.0) * (KInverse(0,0)*jphiuHdiv(0,0) + KInverse(0,1)*jphiuHdiv(1,0) + KInverse(0,2)*jphiuHdiv(2,0));
            oneoverlambda_Kinv_jphiuHdiv(1,0) = (1.0/1.0) * (KInverse(1,0)*jphiuHdiv(0,0) + KInverse(1,1)*jphiuHdiv(1,0) + KInverse(1,2)*jphiuHdiv(2,0));
            oneoverlambda_Kinv_jphiuHdiv(2,0) = (1.0/1.0) * (KInverse(2,0)*jphiuHdiv(0,0) + KInverse(2,1)*jphiuHdiv(1,0) + KInverse(2,2)*jphiuHdiv(2,0));
            
            ek(iq + iniu,jq + iniu) += weight * ((oneoverlambda_Kinv_jphiuHdiv(0,0)*iphiuHdiv(0,0) + oneoverlambda_Kinv_jphiuHdiv(1,0)*iphiuHdiv(1,0)) + oneoverlambda_Kinv_jphiuHdiv(2,0)*iphiuHdiv(2,0));
        }
        
        
        // du/dalphau terms
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(iq + iniu, jp + iniP) += -1.0 * weight/JacobianDet * phiPL2(jp,0) * DivergenceOnMaster(iq,0);
        }

    }
    
    
    /* $ - \underset{\Omega}{\int}w\; div\left(\mathbf{q}\right)\partial\Omega $ */
    for (int ip = 0; ip < nphiPL2; ip++)
    {
        // du/dalphau terms
        for (int jq = 0; jq < nphiuHdiv; jq++)
        {
            ek(ip + iniP, jq + iniu) += -1.0 * weight/JacobianDet * DivergenceOnMaster(jq,0) * phiPL2(ip,0);
        }
    }
        
    
    this->Contribute(datavec,weight,ef);
    
}

// Residual contribution
void TRMMixedDarcy::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2         = datavec[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<STATE> dphiuH1   = datavec[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2   = datavec[Pblock].dphix; // Derivative For L2  test functions
    
    TPZFNMatrix<40,STATE> DivergenceOnMaster;
    STATE divflux;
    // Compute the divergence on deformed element by piola contravariant transformation
    this->ComputeDivergenceOnMaster(datavec, DivergenceOnMaster,divflux);
    REAL JacobianDet = datavec[ublock].detjac;
    
    // Blocks dimensions and lengths
    int nphiuHdiv   = datavec[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2     = phiPL2.Rows();                                    // For L2   P
    int iniu    = 0;
    int iniP    = nphiuHdiv     + iniu;    
    
    // Getting linear combinations from different approximation spaces
    TPZManVector<REAL,3> u      = datavec[ublock].sol[0];
    REAL P              = datavec[Pblock].sol[0][0];
    
    // Get the pressure at the integrations points
    long global_point_index = datavec[0].intGlobPtIndex;
    TRMMemory &point_memory = GetMemory()[global_point_index];
//    STATE pressure = point_memory.GetPressure();
//    STATE rhs = point_memory.GetRhs();
//    STATE w = point_memory.GetWeight();
//    STATE det = point_memory.GetDetJac();
    TPZManVector<STATE> flux = point_memory.GetTotal_Flux();

//    std::cout << "flux = " << flux << std::endl;
//    std::cout << "u    = " << u << std::endl;
//    std::cout << "flux x difference = " << flux[0] - u[0] << std::endl;
//    std::cout << "flux y difference = " << flux[1] - u[1] << std::endl;
//    std::cout << "flux z difference = " << flux[2] - u[2] << std::endl;
//    std::cout << "Pressure difference = " << pressure - P << std::endl;
    
    TPZFMatrix<STATE> Graduaxes = datavec[ublock].dsol[0]; // Piola divengence may works, needed set piola computation on the solution elchiv method!!!
    TPZFMatrix<STATE> GradPaxes = datavec[Pblock].dsol[0];
    
    TPZFNMatrix<660> GradP;
    TPZAxesTools<REAL>::Axes2XYZ(GradPaxes, GradP, datavec[Pblock].axes);

    
    // Rock and fluids parameters
    TPZFMatrix<STATE> KInverse(3,3,0.0);
    KInverse(0,0) = 1.0;
    KInverse(1,1) = 1.0;
    KInverse(2,2) = 1.0;
    
    
    // Defining local variables
    TPZFNMatrix<3,STATE> oneoverlambda_Kinv_u(3,1);
    TPZFNMatrix<3,STATE> Gravity(3,1);
    
    oneoverlambda_Kinv_u(0,0) = (1.0/1.0)* (KInverse(0,0)*u[0] + KInverse(0,1)*u[1] + KInverse(0,2)*u[2]);
    oneoverlambda_Kinv_u(1,0) = (1.0/1.0)* (KInverse(1,0)*u[0] + KInverse(1,1)*u[1] + KInverse(1,2)*u[2]);
    oneoverlambda_Kinv_u(2,0) = (1.0/1.0)* (KInverse(2,0)*u[0] + KInverse(2,1)*u[1] + KInverse(2,2)*u[2]);
    
    Gravity(0,0) = -0.0;
    Gravity(1,0) = -0.0;
    Gravity(2,0) = -0.0;
    
    REAL divu = 0.0;
    TPZFNMatrix<3,STATE> iphiuHdiv(3,1);
    int ishapeindex;
    int ivectorindex;
    
    for (int iq = 0; iq < nphiuHdiv; iq++)
    {
        
        /* $ \underset{\Omega_{e}}{\int}\left(K\lambda\right)^{-1}\mathbf{q}\cdot\mathbf{v}\;\partial\Omega_{e}-\underset{\Omega_{e}}{\int}P\; div\left(\mathbf{v}\right)\partial\Omega-\underset{\Omega_{e}}{\int}\nabla\left(\rho_{f}g\; z\right)\cdot\mathbf{v} $ */
        
        ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
        
        iphiuHdiv(0,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(0,ivectorindex);
        iphiuHdiv(1,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(1,ivectorindex);
        iphiuHdiv(2,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(2,ivectorindex);
        
        ef(iq + iniu) += weight * ((oneoverlambda_Kinv_u(0,0)*iphiuHdiv(0,0) + oneoverlambda_Kinv_u(1,0)*iphiuHdiv(1,0) + oneoverlambda_Kinv_u(2,0)*iphiuHdiv(2,0)));
        ef(iq + iniu) += weight/JacobianDet *(-P) * DivergenceOnMaster(iq,0);
        
    }
    
    TPZManVector<STATE,1> fvalue(1,0.0);
    if(fForcingFunction)
    {
        fForcingFunction->Execute(datavec[Pblock].x,fvalue);
    }
    
    
    divu = (Graduaxes(0,0) + Graduaxes(1,1) + Graduaxes(2,2)); // uses this for constant jacobian elements
//    REAL divu2 = point_memory.GetDiv_Flux();
    
//    std::cout << "divu = " << divu << std::endl;
//    std::cout << "divflux = " << divflux << std::endl;
//    std::cout << "divu2 = " << divu2 << std::endl;
//    std::cout << "rhs = " << rhs << std::endl;
    
//    std::cout << "fvalue[0] = " << fvalue[0] << std::endl;
//    std::cout << "rhs = " << rhs << std::endl;
//    std::cout << "diff = " << fvalue[0] - rhs << std::endl;
    
    /* $ - \underset{\Omega}{\int}w\; div\left(\mathbf{q}\right)\partial\Omega $ */
    for (int ip = 0; ip < nphiPL2; ip++)
    {
        ef(ip + iniP) += -1.0 * weight * (divu - fvalue[0]) * phiPL2(ip,0);
    }
    
    return;
    
}

void TRMMixedDarcy::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    std::cout << " This method should be called only for Capillary pressure terms " << std::endl;
    DebugStop();
    
}

void TRMMixedDarcy::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef)
{
    std::cout << " This method should be called only for Capillary pressure terms " << std::endl;
    DebugStop();
}


void TRMMixedDarcy::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    std::cout << " This method should be called only for Capillary pressure terms " << std::endl;
    DebugStop();
}

void TRMMixedDarcy::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    std::cout << " This method should be called only for Capillary pressure terms " << std::endl;
    DebugStop();
}


void TRMMixedDarcy::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
    // At each Integration Point.
    
    int Qblock = 0;
    int Pblock = 1;

    // Getting test and basis functions
    TPZFNMatrix<9,STATE> PhiH1 = datavec[Qblock].phi; // For H1   test functions
    TPZFNMatrix<9,STATE> WL2   = datavec[Pblock].phi; // For HL2  test functions
    
    TPZFMatrix<STATE> dPhiH1 = datavec[Qblock].dphix; // Derivative For H1   test functions
    TPZFMatrix<STATE> dWL2   = datavec[Pblock].dphix; // Derivative For HL2  test functions
    
    // Getting Linear combinations of basis functions
    TPZManVector<STATE> Q = datavec[Qblock].sol[0];
    TPZManVector<STATE> P = datavec[Pblock].sol[0];
    
    
    // Computing normal flux
    STATE Qn = Q[0];
    
    TPZFMatrix<STATE> dQdx = datavec[Qblock].dsol[0];
    TPZFMatrix<STATE> dPdx = datavec[Pblock].dsol[0];
    
    // Number of phis
    int nPhiHdiv = PhiH1.Rows();  // For Hdiv
//    int nPhiL2   = WL2.Rows();                                  // For L2
    TPZVec<STATE> &x = datavec[0].x;
    
    REAL Value = bc.Val2()(0,0);
    if (bc.HasfTimedependentBCForcingFunction()) {
        TPZManVector<STATE,2> force(1);
        REAL time = 0.0;
        TPZFMatrix<double> gradf;
        bc.TimedependentBCForcingFunction()->Execute(x, time, force, gradf);
        Value = force[0];
    }
    else{
        Value = bc.Val2()(0,0);
    }
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD inflow
        {
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                ef(iq) += weight * ( (Value ) * PhiH1(iq,0));
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN inflow
        {
            
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                
                ef(iq) += weight * (gBigNumber * (Qn - Value)) * PhiH1(iq,0);
                
                for (int jq = 0; jq < nPhiHdiv; jq++)
                {
                    
                    
                    ek(iq,jq) += weight * ( gBigNumber * PhiH1(jq,0)) * PhiH1(iq,0);
                }
                
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


int TRMMixedDarcy::ClassId() const {
    return -6378;
}

// -------------------------------------------------------------------------------------------

void TRMMixedDarcy::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
}

// -------------------------------------------------------------------------------------------

void TRMMixedDarcy::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}

// Update element memory by copying the n+1 data to the n data
void TRMMixedDarcy::UpdateMemory()
{
    long nel = fMemory.NElements();
    for (long el=0; el<nel; el++) {
        fMemory[el].UpdateSolutionMemory();
    }
}

