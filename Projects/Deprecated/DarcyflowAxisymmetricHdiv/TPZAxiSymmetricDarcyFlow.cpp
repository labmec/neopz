/*
 *  TPZAxiSymmetricDarcyFlow.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZAxiSymmetricDarcyFlow.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzfmatrix.h"

TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow() : TPZDiscontinuousGalerkin()
{
    
    fepsilon = 10.0e-8;
    fSalpha_max = 0.387426;
    
    fSimulationData=NULL;
    fReservoirdata=NULL;
    fPetrophysicdata=NULL;
    fluid_alpha=NULL;
    fluid_beta=NULL;
    fluid_gamma=NULL;
    
    fnstate_vars = 0;
    fBulkVelocity.Resize(3,0.0);
    fAveragePressure = 0.0;
    fWaterSaturation = 0.0;
    fOilSaturation = 0.0;
    
    fWaterVelocity.Resize(3,0.0);
    fOilVelocity.Resize(3,0.0);
    fGasVelocity.Resize(3,0.0);
    
    fWaterPressure.Resize(4,0.0);
    fOilPressure.Resize(4,0.0);
    fGasPressure.Resize(4,0.0);
    
    fWaterDensity.Resize(4,0.0);
    fOilDensity.Resize(4,0.0);
    fGasDensity.Resize(4,0.0);
    
    fFWater.Resize(4,0.0);
    fFOil.Resize(4,0.0);
    fFGas.Resize(4,0.0);
    
    fWaterMobility.Resize(4,0.0);
    fOilMobility.Resize(4,0.0);
    fGasMobility.Resize(4,0.0);
    
    fTotalMobility.Resize(4,0.0);
    fTotalMobility.Resize(4,0.0);
    
}

TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow(int matid) : TPZDiscontinuousGalerkin(matid)
{
    fepsilon = 10.0e-8;
    fSalpha_max = 0.387426;
    
    fSimulationData=NULL;
    fReservoirdata=NULL;
    fPetrophysicdata=NULL;
    
    fBulkVelocity.Resize(3,0.0);
    fAveragePressure = 0.0;
    fWaterSaturation = 0.0;
    fOilSaturation = 0.0;
    
    fWaterVelocity.Resize(3,0.0);
    fOilVelocity.Resize(3,0.0);
    fGasVelocity.Resize(3,0.0);
    
    fWaterPressure.Resize(4,0.0);
    fOilPressure.Resize(4,0.0);
    fGasPressure.Resize(4,0.0);
    
    fWaterDensity.Resize(4,0.0);
    fOilDensity.Resize(4,0.0);
    fGasDensity.Resize(4,0.0);
    
    fFWater.Resize(4,0.0);
    fFOil.Resize(4,0.0);
    fFGas.Resize(4,0.0);
    
    fWaterMobility.Resize(4,0.0);
    fOilMobility.Resize(4,0.0);
    fGasMobility.Resize(4,0.0);
    
    fTotalMobility.Resize(4,0.0);
    fTotalMobility.Resize(4,0.0);
}


TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow(const TPZAxiSymmetricDarcyFlow &mat) : TPZDiscontinuousGalerkin(mat)
{
    fSimulationData     = mat.fSimulationData;
    fReservoirdata      = mat.fReservoirdata;
    fPetrophysicdata    = mat.fPetrophysicdata;
    
    fBulkVelocity       = mat.fBulkVelocity;
    fAveragePressure    = mat.fAveragePressure;
    fWaterSaturation    = mat.fWaterSaturation;
    fOilSaturation      = mat.fOilSaturation;
    
    fWaterVelocity      = mat.fWaterVelocity;
    fOilVelocity        = mat.fOilVelocity;
    fGasVelocity        = mat.fGasVelocity;
    
    fWaterPressure      = mat.fWaterPressure;
    fOilPressure        = mat.fOilPressure;
    fGasPressure        = mat.fGasPressure;
    
    fWaterDensity       = mat.fWaterDensity;
    fOilDensity         = mat.fOilDensity;
    fGasDensity         = mat.fGasDensity;
    
    fFWater             = mat.fFWater;
    fFOil               = mat.fFOil;
    fFGas               = mat.fFGas;
    
    fWaterMobility      = mat.fWaterMobility;
    fOilMobility        = mat.fOilMobility;
    fGasMobility        = mat.fGasMobility;
    
    fTotalMobility      = mat.fTotalMobility;
    fTotalDensity       = mat.fTotalDensity;
}

TPZAxiSymmetricDarcyFlow::~TPZAxiSymmetricDarcyFlow()
{
    
}

void TPZAxiSymmetricDarcyFlow::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(true);
        datavec[idata].fNeedsSol = true;
    }
}

void TPZAxiSymmetricDarcyFlow::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(true);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TPZAxiSymmetricDarcyFlow::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

int TPZAxiSymmetricDarcyFlow::VariableIndex(const std::string &name) {
    if (!strcmp("P", name.c_str())) return 0;
    if (!strcmp("u", name.c_str())) return 1;
    if (!strcmp("S_alpha", name.c_str())) return 2;
    if (!strcmp("S_beta", name.c_str())) return 3;
    if (!strcmp("S_gamma", name.c_str())) return 4;
    if (!strcmp("Rho_alpha", name.c_str())) return 5;
    if (!strcmp("Rho_beta", name.c_str())) return 6;
    if (!strcmp("Rho_gamma", name.c_str())) return 7;
    if (!strcmp("Porosity", name.c_str())) return 8;
    if (!strcmp("div_u", name.c_str())) return 9;
    if (!strcmp("Exact_S", name.c_str())) return 10;
    if (!strcmp("Exact_GradS", name.c_str())) return 11;
    if (!strcmp("Rhs", name.c_str())) return 12;
    if (!strcmp("Pc_beta_alpha", name.c_str())) return 13;
    if (!strcmp("P_alpha", name.c_str())) return 14;
    if (!strcmp("P_beta", name.c_str())) return 15;
    if (!strcmp("u_alpha", name.c_str())) return 16;
    if (!strcmp("u_beta", name.c_str())) return 17;
    if (!strcmp("u_gamma", name.c_str())) return 18;
    if (!strcmp("u_alpha_sc", name.c_str())) return 19;
    if (!strcmp("u_beta_sc", name.c_str())) return 20;
    if (!strcmp("u_gamma_sc", name.c_str())) return 21;
    if (!strcmp("f_alpha", name.c_str())) return 22;
    if (!strcmp("f_beta", name.c_str())) return 23;
    if (!strcmp("f_gamma", name.c_str())) return 24;
    if (!strcmp("kappa", name.c_str())) return 25;
    std::cout  << " Var index not implemented " << std::endl;
    DebugStop();
    return 0;
}

int TPZAxiSymmetricDarcyFlow::NSolutionVariables(int var) {
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
            return 1; // Scalar
        case 5:
            return 1; // Scalar
        case 6:
            return 1; // Scalar
        case 7:
            return 1; // Scalar
        case 8:
            return 1; // Scalar
        case 9:
            return 1; // Scalar
        case 10:
            return 1; // Scalar
        case 11:
            return 3; // Vector
        case 12:
            return 1; // Scalar
        case 13:
            return 1; // Scalar
        case 14:
            return 1; // Scalar
        case 15:
            return 1; // Scalar
        case 16:
            return 3; // Vector
        case 17:
            return 3; // Vector
        case 18:
            return 3; // Vector
        case 19:
            return 3; // Vector
        case 20:
            return 3; // Vector
        case 21:
            return 3; // Vector
        case 22:
            return 1; // Scalar
        case 23:
            return 1; // Scalar
        case 24:
            return 1; // Scalar
        case 25:
            return 3; // Vector
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return 0;
}

void TPZAxiSymmetricDarcyFlow::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int ublock = 0;
    int Pblock = 1;
    int Sablock = 2;
    int Sbblock = 3;
    
    TPZManVector<REAL,3> u = datavec[ublock].sol[0];
    REAL P = datavec[Pblock].sol[0][0];
    TPZFMatrix<STATE> dudx = datavec[ublock].dsol[0];
    TPZFMatrix<STATE> dPdx = datavec[Pblock].dsol[0];
    TPZManVector<REAL> state_vars(2,0.0);
    state_vars[1] = P;
    
    REAL S_alpha = 0.0;
    REAL S_beta = 0.0;
    REAL S_gamma = 0.0;
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    this->ComputeProperties(datavec, props);
    TPZManVector<REAL> rho_alpha        = props[0];
    TPZManVector<REAL> rho_beta         = props[1];
    TPZManVector<REAL> rho_gamma        = props[2];
    TPZManVector<REAL> f_alpha          = props[3];
    TPZManVector<REAL> f_beta           = props[4];
    TPZManVector<REAL> f_gamma          = props[5];
    TPZManVector<REAL> l                = props[6];
    TPZManVector<REAL> rho              = props[7];
    TPZManVector<REAL> rhof             = props[8];
    TPZManVector<REAL> Pc_beta_alpha    = props[9];

    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = datavec[0].x[0];
    REAL r = Norm(x_spatial);
    REAL s = 1.0;
    
    if (fSimulationData->IsAxisymmetricQ()) {
        s *= 2.0*M_PI*r;
        u[0] *= (1.0/s);
        u[1] *= (1.0/s);
        u[2] *= (1.0/s);
        dudx *= (1.0/s);
    }
    

    
    TPZManVector<STATE> GradS(2,0.0);
    TPZManVector<STATE> Grad_Pc(2,0.0);
    
    if (fSimulationData->IsTwoPhaseQ()) {
        state_vars.Resize(3);
        
        fluid_beta->Density(rho_beta, state_vars);
        S_alpha = datavec[Sablock].sol[0][0];
        state_vars[2] = S_alpha;
        S_beta = 1.0 - S_alpha;
        
        TPZFMatrix<STATE> GradSaxes = datavec[Sablock].dsol[0];
        //  Compute grad(S)
        GradS[0] = GradSaxes(0,0)*datavec[Sablock].axes(0,0)+GradSaxes(1,0)*datavec[Sablock].axes(1,0);
        GradS[1] = GradSaxes(0,0)*datavec[Sablock].axes(0,1)+GradSaxes(1,0)*datavec[Sablock].axes(1,1);
        
        //  Compute grad(Pc)
        Grad_Pc[0] = Pc_beta_alpha[3] * GradS[0];
        Grad_Pc[1] = Pc_beta_alpha[3] * GradS[1];
        
    }
    
    if (fSimulationData->IsThreePhaseQ()) {
        state_vars.Resize(4);
        
        fluid_beta->Density(rho_beta, state_vars);
        fluid_gamma->Density(rho_gamma, state_vars);
        S_alpha = datavec[Sablock].sol[0][0];
        S_beta = datavec[Sbblock].sol[0][0];
        state_vars[2] = S_alpha;
        state_vars[3] = S_beta;
    }
    
    
    S_gamma = 1.0 - S_alpha - S_beta;

    REAL time = fSimulationData->GetTime();
    

    TPZFMatrix<STATE> K = fReservoirdata->Kabsolute();
    TPZFMatrix<STATE> G = fSimulationData->GetGravity();
    
    if (fSimulationData->IsHeterogeneousQ()) {
        int gelid = datavec[0].gelElId;
        K = fReservoirdata->GetKvector()[gelid];
    }
    
    REAL rock_phi;
    REAL drock_phidP;
    this->fReservoirdata->Porosity(P, rock_phi, drock_phidP);
    

    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0:
        {
            Solout[0] = P;
        }
            break;
        case 1:
        {
            Solout[0] = u[0]; // Bulk mass velocity
            Solout[1] = u[1]; // Bulk mass velocity
        }
            break;
        case 2:
        {
            Solout[0] = S_alpha;
        }
            break;
        case 3:
        {
            Solout[0] = S_beta;
        }
            break;
        case 4:
        {
            Solout[0] = S_gamma;
        }
            break;
        case 5:
        {
            Solout[0] = rho_alpha[0];
        }
            break;
        case 6:
        {
            Solout[0] = rho_beta[0];
        }
            break;
        case 7:
        {
            Solout[0] = rho_gamma[0];
        }
            break;
        case 8:
        {
            Solout[0] = rock_phi;
        }
            break;
        case 9:
        {
            Solout[0] = (dudx(0,0) + dudx(1,1) + dudx(2,2));
        }
            break;
        case 10:
        {
            TPZVec<STATE> S(1,0.0);
            TPZFMatrix<STATE> GradS(3,1,0.0);
            fTimedependentFunctionExact->Execute(datavec[ublock].x, fSimulationData->GetTime(), S,GradS);
            Solout[0] = S[0];
            
        }
            break;
        case 11:
        {
            
            TPZVec<STATE> S(1,0.0);
            TPZFMatrix<STATE> GradS(3,1,0.0);
            fTimedependentFunctionExact->Execute(datavec[ublock].x, fSimulationData->GetTime(), S,GradS);
            Solout[0] = GradS(0,0);
            Solout[1] = GradS(1,0);
            Solout[2] = GradS(2,0);
            
//            //            REAL epsilon = 0.01;
//            //            REAL xc = 0.5;
//            //            REAL yc = 0.5;
//            REAL x = datavec[Pblock].x[0];
//            REAL y = datavec[Pblock].x[1];
//            //            REAL Pressure = exp(-((x-xc)*(x-xc)+(y-yc)*(y-yc))/epsilon);
//            
//            REAL Pressure = x*y*(1-x)*(1-y)*(sin(M_PI*x))*(cos(M_PI*y));
//            Solout[0] = Pressure;
            
        }
            break;
        case 12:
        {
            TPZManVector<STATE,1> fvalue(1,0.0);
            TPZFMatrix<REAL> Grad;
            if(fTimeDependentForcingFunction)
            {
                fTimeDependentForcingFunction->Execute(datavec[Pblock].x,fSimulationData->GetTime(),fvalue,Grad);
            }
            Solout[0] = fvalue[0];
        }
            break;
        case 13:
        {
            Solout[0] = Pc_beta_alpha[0];
        }
            break;
        case 14:
        {
            Solout[0] = P - S_beta * (Pc_beta_alpha[0]);
        }
            break;
        case 15:
        {
            Solout[0] = P + S_alpha * (Pc_beta_alpha[0]);
        }
            break;
        case 16:
        {
            Solout[0] = f_alpha[0]*u[0] + l[0]*f_alpha[0]*f_beta[0] * (K(0,0)*Grad_Pc[0] + K(0,1)*Grad_Pc[1]) + l[0]*f_alpha[0]*f_beta[0] * (rho_alpha[0]-rho_beta[0]) * (K(0,0)*G(0,0) + K(0,1)*G(1,0));
            Solout[1] = f_alpha[0]*u[1] + l[0]*f_alpha[0]*f_beta[0] * (K(1,0)*Grad_Pc[0] + K(1,1)*Grad_Pc[1]) + l[0]*f_alpha[0]*f_beta[0] * (rho_alpha[0]-rho_beta[0]) * (K(1,0)*G(0,0) + K(1,1)*G(1,0));
        }
            break;
        case 17:
        {
            Solout[0] = f_beta[0]*u[0] - l[0]*f_alpha[0]*f_beta[0] * (K(0,0)*Grad_Pc[0] + K(0,1)*Grad_Pc[1]) - l[0]*f_alpha[0]*f_beta[0] * (rho_alpha[0]-rho_beta[0]) * (K(0,0)*G(0,0) + K(0,1)*G(1,0));
            Solout[1] = f_beta[0]*u[1] - l[0]*f_alpha[0]*f_beta[0] * (K(1,0)*Grad_Pc[0] + K(1,1)*Grad_Pc[1]) - l[0]*f_alpha[0]*f_beta[0] * (rho_alpha[0]-rho_beta[0]) * (K(1,0)*G(0,0) + K(1,1)*G(1,0));
        }
            break;
        case 18:
        {
            Solout[0] = f_alpha[0]*u[0] + l[0]*f_alpha[0]*f_beta[0] * (K(0,0)*Grad_Pc[0] + K(0,1)*Grad_Pc[1]) + l[0]*f_alpha[0]*f_beta[0] * (rho_alpha[0]-rho_beta[0]) * (K(0,0)*G(0,0) + K(0,1)*G(1,0));
            Solout[1] = f_alpha[0]*u[1] + l[0]*f_alpha[0]*f_beta[0] * (K(1,0)*Grad_Pc[0] + K(1,1)*Grad_Pc[1]) + l[0]*f_alpha[0]*f_beta[0] * (rho_alpha[0]-rho_beta[0]) * (K(1,0)*G(0,0) + K(1,1)*G(1,0));
        }
            break;
        case 19:
        {
            Solout[0] = f_alpha[0]*u[0] + l[0]*f_alpha[0]*f_beta[0] * (K(0,0)*Grad_Pc[0] + K(0,1)*Grad_Pc[1]) + l[0]*f_alpha[0]*f_beta[0] * (rho_alpha[0]-rho_beta[0]) * (K(0,0)*G(0,0) + K(0,1)*G(1,0));
            Solout[1] = f_alpha[0]*u[1] + l[0]*f_alpha[0]*f_beta[0] * (K(1,0)*Grad_Pc[0] + K(1,1)*Grad_Pc[1]) + l[0]*f_alpha[0]*f_beta[0] * (rho_alpha[0]-rho_beta[0]) * (K(1,0)*G(0,0) + K(1,1)*G(1,0));
            
            Solout[0] *= (1.0)/(rho_alpha[0]);
            Solout[1] *= (1.0)/(rho_alpha[0]);
            
            REAL B_alpha = fluid_alpha->GetRho()/rho_alpha[0];
            Solout[0] *= (1.0)/(B_alpha);
            Solout[1] *= (1.0)/(B_alpha);
            
        }
            break;
        case 20:
        {
            Solout[0] = f_beta[0]*u[0] - l[0]*f_alpha[0]*f_beta[0] * (K(0,0)*Grad_Pc[0] + K(0,1)*Grad_Pc[1]) - l[0]*f_alpha[0]*f_beta[0] * (rho_alpha[0]-rho_beta[0]) * (K(0,0)*G(0,0) + K(0,1)*G(1,0));
            Solout[1] = f_beta[0]*u[1] - l[0]*f_alpha[0]*f_beta[0] * (K(1,0)*Grad_Pc[0] + K(1,1)*Grad_Pc[1]) - l[0]*f_alpha[0]*f_beta[0] * (rho_alpha[0]-rho_beta[0]) * (K(1,0)*G(0,0) + K(1,1)*G(1,0));
            
            Solout[0] *= (1.0)/(rho_beta[0]);
            Solout[1] *= (1.0)/(rho_beta[0]);
            
            REAL B_beta = fluid_beta->GetRho()/rho_beta[0];
            Solout[0] *= (1.0)/(B_beta);
            Solout[1] *= (1.0)/(B_beta);
            
        }
            break;
        case 21:
        {
            Solout[0] = f_alpha[0]*u[0] + l[0]*f_alpha[0]*f_beta[0] * (K(0,0)*Grad_Pc[0] + K(0,1)*Grad_Pc[1]) + l[0]*f_alpha[0]*f_beta[0] * (rho_alpha[0]-rho_beta[0]) * (K(0,0)*G(0,0) + K(0,1)*G(1,0));
            Solout[1] = f_alpha[0]*u[1] + l[0]*f_alpha[0]*f_beta[0] * (K(1,0)*Grad_Pc[0] + K(1,1)*Grad_Pc[1]) + l[0]*f_alpha[0]*f_beta[0] * (rho_alpha[0]-rho_beta[0]) * (K(1,0)*G(0,0) + K(1,1)*G(1,0));
            
            Solout[0] *= (1.0)/(rho_beta[0]);
            Solout[1] *= (1.0)/(rho_beta[0]);
            
            REAL B_beta = fluid_beta->GetRho()/rho_beta[0];
            Solout[0] *= (rock_phi)/(B_beta);
            Solout[1] *= (rock_phi)/(B_beta);
            
        }
        case 22:
        {
            Solout[0] = f_alpha[0];
        }
            break;
        case 23:
        {
            Solout[0] = f_beta[0];
        }
            break;
        case 24:
        {
            Solout[0] = f_gamma[0];
        }
        case 25:
        {
            Solout[0] = K(0,0);
            Solout[1] = K(1,1);
        }
            break;
            break;
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
}

// Divergence on deformed element
void TPZAxiSymmetricDarcyFlow::ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi)
{
    int ublock = 0;
    
   
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;   // For H1  test functions Q
    TPZFMatrix<REAL> dphiuH1       = datavec[ublock].dphi; // Derivative For H1  test functions
    TPZFMatrix<REAL> dphiuH1axes   = datavec[ublock].dphix; // Derivative For H1  test functions
    
    TPZFNMatrix<660,REAL> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[ublock].axes);
    
    int nphiuHdiv = datavec[ublock].fVecShapeIndex.NElements();
    
    DivergenceofPhi.Resize(nphiuHdiv,1);
    
    REAL JacobianDet = datavec[ublock].detjac;
    
    TPZFMatrix<REAL> Qaxes = datavec[ublock].axes;
    TPZFMatrix<REAL> QaxesT;
    TPZFMatrix<REAL> Jacobian = datavec[ublock].jacobian;
    TPZFMatrix<REAL> JacobianInverse = datavec[ublock].jacinv;
    
    TPZFMatrix<REAL> GradOfX;
    TPZFMatrix<REAL> GradOfXInverse;
    TPZFMatrix<REAL> VectorOnMaster;
    TPZFMatrix<REAL> VectorOnXYZ(3,1,0.0);
    Qaxes.Transpose(&QaxesT);
    QaxesT.Multiply(Jacobian, GradOfX);
    JacobianInverse.Multiply(Qaxes, GradOfXInverse);
    
    int ivectorindex = 0;
    int ishapeindex = 0;
    
    if (HDivPiola)
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
            DivergenceofPhi(iq,0) =  (1.0/JacobianDet) * (dphiuH1(0,ishapeindex)*VectorOnMaster(0,0) + dphiuH1(1,ishapeindex)*VectorOnMaster(1,0) );
        }
    }
    else
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            /* Computing the divergence for constant jacobian elements */
            DivergenceofPhi(iq,0) =  (datavec[ublock].fNormalVec(0,ivectorindex)*GradphiuH1(0,ishapeindex) + datavec[ublock].fNormalVec(1,ivectorindex)*GradphiuH1(1,ishapeindex));
        }
    }
    
    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = datavec[0].x[0];
    REAL r = Norm(x_spatial);
    REAL s = 1.0;
    
    if (fSimulationData->IsAxisymmetricQ()) {
        s *= 2.0*M_PI*r;
        DivergenceofPhi *= 1.0/s;
    }
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = datavec[0].x[0];
    REAL r = Norm(x_spatial);
    REAL s = 1.0;
    
     if (fSimulationData->IsAxisymmetricQ()) {
         s *= 2.0*M_PI*r;
         
         weight *= s;
         datavec[0].phi *= 1.0/s;
         datavec[0].sol[0][0]   *= 1.0/s;
         datavec[0].sol[0][1]   *= 1.0/s;
         datavec[0].sol[0][2]   *= 1.0/s;
         datavec[0].dsol[0]   *= 1.0/s;

         
     }
    
    this->ContributeDarcy(datavec, weight, ek, ef);
    
    if (fSimulationData->IsTwoPhaseQ()) {
        this->ContributeAlpha(datavec, weight, ek, ef);
    }
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = datavec[0].x[0];
    REAL r = Norm(x_spatial);
    REAL s = 1.0;
    
    if (fSimulationData->IsAxisymmetricQ()) {
        s *= 2.0*M_PI*r;
        
        weight *= s;
        datavec[0].phi *= 1.0/s;
        datavec[0].sol[0][0]   *= 1.0/s;
        datavec[0].sol[0][1]   *= 1.0/s;
        datavec[0].sol[0][2]   *= 1.0/s;
        datavec[0].dsol[0]   *= 1.0/s;
        
        
    }
    
    this->ContributeDarcy(datavec, weight, ef);
    
    if (fSimulationData->IsTwoPhaseQ()) {
        this->ContributeAlpha(datavec, weight, ef);
    }
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    
    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = data.x[0];
    REAL r = Norm(x_spatial);
    REAL s = 1.0;
    
    if (fSimulationData->IsAxisymmetricQ()) {
        s *= 2.0*M_PI*r;
        
        weight *= s;
        datavecleft[0].phi *= 1.0/s;
        datavecleft[0].sol[0][0]   *= 1.0/s;
        datavecleft[0].sol[0][1]   *= 1.0/s;
        datavecleft[0].sol[0][2]   *= 1.0/s;
        datavecleft[0].dsol[0]   *= 1.0/s;
        
        datavecright[0].phi *= 1.0/s;
        datavecright[0].sol[0][0]   *= 1.0/s;
        datavecright[0].sol[0][1]   *= 1.0/s;
        datavecright[0].sol[0][2]   *= 1.0/s;
        datavecright[0].dsol[0]   *= 1.0/s;

        
    }
    
    this->ContributeInterfaceDarcy(data, datavecleft, datavecright, weight, ek, ef);
    
    if (fSimulationData->IsTwoPhaseQ()) {
        this->ContributeInterfaceAlpha(data, datavecleft, datavecright, weight, ek, ef);
    }
    
    
}

void TPZAxiSymmetricDarcyFlow::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef)
{
 
    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = data.x[0];
    REAL r = Norm(x_spatial);
    REAL s = 1.0;
    
    if (fSimulationData->IsAxisymmetricQ()) {
        s *= 2.0*M_PI*r;
        
        weight *= s;
        datavecleft[0].phi *= 1.0/s;
        datavecleft[0].sol[0][0]   *= 1.0/s;
        datavecleft[0].sol[0][1]   *= 1.0/s;
        datavecleft[0].sol[0][2]   *= 1.0/s;
        datavecleft[0].dsol[0]   *= 1.0/s;
        
        datavecright[0].phi *= 1.0/s;
        datavecright[0].sol[0][0]   *= 1.0/s;
        datavecright[0].sol[0][1]   *= 1.0/s;
        datavecright[0].sol[0][2]   *= 1.0/s;
        datavecright[0].dsol[0]   *= 1.0/s;
        
        
    }
    
    this->ContributeInterfaceDarcy(data, datavecleft, datavecright, weight, ef);
    
    
    if (fSimulationData->IsTwoPhaseQ()) {
        this->ContributeInterfaceAlpha(data, datavecleft, datavecright, weight, ef);
    }
    
    return;
}

void TPZAxiSymmetricDarcyFlow::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = data.x[0];
    REAL r = Norm(x_spatial);
    REAL s = 1.0;
    
    if (fSimulationData->IsAxisymmetricQ()) {
        s *= 2.0*M_PI*r;
        
        weight *= s;
        datavecleft[0].phi *= 1.0/s;
        datavecleft[0].sol[0][0]   *= 1.0/s;
        datavecleft[0].sol[0][1]   *= 1.0/s;
        datavecleft[0].sol[0][2]   *= 1.0/s;
        datavecleft[0].dsol[0]   *= 1.0/s;
        
        
    }
    
    if (fSimulationData->IsTwoPhaseQ()) {
        this->ContributeBCInterfaceAlpha(data, datavecleft, weight, ek, ef, bc);
    }
    
    return;
    
    
}


void TPZAxiSymmetricDarcyFlow::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = data.x[0];
    REAL r = Norm(x_spatial);
    REAL s = 1.0;
    
    if (fSimulationData->IsAxisymmetricQ()) {
        s *= 2.0*M_PI*r;
        
        weight *= s;
        datavecleft[0].phi *= 1.0/s;
        datavecleft[0].sol[0][0]   *= 1.0/s;
        datavecleft[0].sol[0][1]   *= 1.0/s;
        datavecleft[0].sol[0][2]   *= 1.0/s;
        datavecleft[0].dsol[0]   *= 1.0/s;
        
        
    }
    
    if (fSimulationData->IsTwoPhaseQ()) {
        this->ContributeBCInterfaceAlpha(data, datavecleft, weight, ef, bc);
    }
    
    return;
    
}


void TPZAxiSymmetricDarcyFlow::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = datavec[0].x[0];
    REAL r = Norm(x_spatial);
    REAL s = 1.0;
    
    if (fSimulationData->IsAxisymmetricQ()) {
        s *= 2.0*M_PI*r;
        
        weight *= s;
        datavec[0].phi *= 1.0/s;
        datavec[0].sol[0][0]   *= 1.0/s;
        
    }
    
    this->ContributeBCDarcy(datavec, weight, ek, ef, bc);
    return;
    
}




void TPZAxiSymmetricDarcyFlow::ContributeDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // Getting data for Mixed-Darcy flow problem
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2         = datavec[Pblock].phi;  // For L2  test functions P
    
    TPZFMatrix<REAL> dphiuH1   = datavec[ublock].dphix; // Derivative For H1  test functions u
    TPZFMatrix<REAL> dphiPL2   = datavec[Pblock].dphix; // Derivative For L2  test functions P
    
    TPZFMatrix<STATE> DivergenceOnDeformed;
    // Compute the divergence on deformed element by piola contravariant transformation
    this->ComputeDivergenceOnDeformed(datavec, DivergenceOnDeformed);
    
    
    // Blocks dimensions and lengths
    int nphiuHdiv   = datavec[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2     = phiPL2.Rows();                                    // For L2   P
    int iniu    = 0;
    int iniP    = nphiuHdiv     + iniu;
    
    // Getting linear combinations from different approximation spaces
    TPZManVector<REAL,3> u      = datavec[ublock].sol[0];
    REAL P              = datavec[Pblock].sol[0][0];
    
    TPZFMatrix<STATE> Graduaxes = datavec[ublock].dsol[0];
    TPZFMatrix<STATE> GradPaxes = datavec[Pblock].dsol[0];
    TPZFNMatrix<660,REAL> GradP;
    TPZAxesTools<REAL>::Axes2XYZ(GradPaxes, GradP, datavec[Pblock].axes);
    
    // Rock and fluids parameters
    TPZFMatrix<STATE> KInverse = fReservoirdata->KabsoluteInv();
    TPZFMatrix<STATE> Oneoverlambda_Kinv_u(2,1);
    TPZFMatrix<REAL> Oneoverlambda_Kinv_Phiu(2,1);
    TPZFMatrix<STATE> Gravity(2,1);
    TPZFMatrix<STATE> gm(2,1);
    TPZFMatrix<STATE> dgmdP(2,1);

    TPZFMatrix<STATE> K;
    if (fSimulationData->IsHeterogeneousQ()) {
        int gelid = datavec[0].gelElId;
        K = fReservoirdata->GetKvector()[gelid];
        KInverse = fReservoirdata->ComputeInvKabsolute(K);
    }
    
    
    Gravity = fSimulationData->GetGravity();
    
    REAL phi, dphidP;
    this->fReservoirdata->Porosity(P, phi, dphidP);
    
    // time
    REAL dt = fSimulationData->GetDeltaT();
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    this->ComputeProperties(datavec, props);
    TPZManVector<REAL> rho_alpha        = props[0];
    TPZManVector<REAL> rho_beta         = props[1];
    TPZManVector<REAL> rho_gamma        = props[2];
    TPZManVector<REAL> f_alpha          = props[3];
    TPZManVector<REAL> f_beta           = props[4];
    TPZManVector<REAL> f_gamma          = props[5];
    TPZManVector<REAL> lambda           = props[6];
    TPZManVector<REAL> rho              = props[7];
    TPZManVector<REAL> rhof             = props[8];
        
    Oneoverlambda_Kinv_u(0,0) = (1.0/lambda[0])* (KInverse(0,0)*u[0] + KInverse(0,1)*u[1]);
    Oneoverlambda_Kinv_u(1,0) = (1.0/lambda[0])* (KInverse(1,0)*u[0] + KInverse(1,1)*u[1]);
    
    gm(0,0) = rhof[0] * Gravity(0,0);
    gm(1,0) = rhof[0] * Gravity(1,0);
    
    dgmdP(0,0) = rhof[2] * Gravity(0,0);
    dgmdP(1,0) = rhof[2] * Gravity(1,0);
    
    REAL divu = 0.0;
    TPZFMatrix<STATE> iphiuHdiv(2,1);
    int ishapeindex;
    int ivectorindex;
    TPZFMatrix<STATE> jphiuHdiv(2,1);
    int jshapeindex;
    int jvectorindex;
    
    for (int iq = 0; iq < nphiuHdiv; iq++)
    {
        
        /* $ \underset{\Omega_{e}}{\int}\left(K\lambda\right)^{-1}\mathbf{q}\cdot\mathbf{v}\;\partial\Omega_{e}-\underset{\Omega_{e}}{\int}P\; div\left(\mathbf{v}\right)\partial\Omega-\underset{\Omega_{e}}{\int}\nabla\left(\rho_{f}g\; z\right)\cdot\mathbf{v} $ */
        
        ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
        
        iphiuHdiv(0,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(0,ivectorindex);
        iphiuHdiv(1,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(1,ivectorindex);
        
        ef(iq + iniu) += weight * ((Oneoverlambda_Kinv_u(0,0)*iphiuHdiv(0,0) + Oneoverlambda_Kinv_u(1,0)*iphiuHdiv(1,0))
                                   - P * DivergenceOnDeformed(iq,0)
                                   - (gm(0,0)*iphiuHdiv(0,0) + gm(1,0)*iphiuHdiv(1,0)) );
        
        // du/dalphau terms
        for (int jq = 0; jq < nphiuHdiv; jq++)
        {
            jvectorindex = datavec[ublock].fVecShapeIndex[jq].first;
            jshapeindex = datavec[ublock].fVecShapeIndex[jq].second;
            
            jphiuHdiv(0,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(0,jvectorindex);
            jphiuHdiv(1,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(1,jvectorindex);
            
            Oneoverlambda_Kinv_Phiu(0,0) = (1.0/lambda[0]) * (KInverse(0,0)*jphiuHdiv(0,0) + KInverse(0,1)*jphiuHdiv(1,0));
            Oneoverlambda_Kinv_Phiu(1,0) = (1.0/lambda[0]) * (KInverse(1,0)*jphiuHdiv(0,0) + KInverse(1,1)*jphiuHdiv(1,0));
            
            ek(iq + iniu,jq + iniu) +=  weight * ((Oneoverlambda_Kinv_Phiu(0,0)*iphiuHdiv(0,0) + Oneoverlambda_Kinv_Phiu(1,0)*iphiuHdiv(1,0)));
        }
        
        // du/dalphap terms
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(iq + iniu,jp + iniP) += weight * ((-lambda[2]/lambda[0])*(Oneoverlambda_Kinv_u(0,0)*iphiuHdiv(0,0) + Oneoverlambda_Kinv_u(1,0)*iphiuHdiv(1,0)) - DivergenceOnDeformed(iq,0) - (dgmdP(0,0)*iphiuHdiv(0,0) + dgmdP(1,0)*iphiuHdiv(1,0)) )* phiPL2(jp,0) ;
        }
        
        
    }
    
    TPZManVector<STATE,1> fvalue(1,0.0);
    TPZFMatrix<REAL> Grad;
    if(fTimeDependentForcingFunction)
    {
        fTimeDependentForcingFunction->Execute(datavec[Pblock].x,fSimulationData->GetTime(),fvalue,Grad);
    }

    divu = (Graduaxes(0,0) + Graduaxes(1,1) + Graduaxes(2,2)); // Note::  use it for HdivPiola = 1 or constant Jacobian mappings
    
    /* $ - \underset{\Omega}{\int}w\; div\left(\mathbf{q}\right)\partial\Omega $ */
    for (int ip = 0; ip < nphiPL2; ip++)
    {
        
        ef(ip + iniP) += weight * (- divu + fvalue[0] - (1.0/dt) * phi * rho[0]) * phiPL2(ip,0);
        
        for (int jq = 0; jq < nphiuHdiv; jq++)
        {
            ek(ip + iniP,jq + iniu) += weight * (- DivergenceOnDeformed(jq,0)) * phiPL2(ip,0);
        }
        
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(ip + iniP,jp + iniP) += weight * (- (1.0/dt) * phi * rho[2] * phiPL2(jp,0) ) * phiPL2(ip,0);
        }
        
    }
    
    
}


void TPZAxiSymmetricDarcyFlow::ContributeDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    // Getting data for Mixed-Darcy flow problem
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2         = datavec[Pblock].phi;  // For L2  test functions P
    
    TPZFMatrix<REAL> dphiuH1   = datavec[ublock].dphix; // Derivative For H1  test functions u
    TPZFMatrix<REAL> dphiPL2   = datavec[Pblock].dphix; // Derivative For L2  test functions P
    
    TPZFMatrix<STATE> DivergenceOnDeformed;
    // Compute the divergence on deformed element by piola contravariant transformation
    this->ComputeDivergenceOnDeformed(datavec, DivergenceOnDeformed);
    
    // Blocks dimensions and lengths
    int nphiuHdiv   = datavec[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2     = phiPL2.Rows();                                    // For L2   P
    int iniu    = 0;
    int iniP    = nphiuHdiv     + iniu;
    
    // Getting linear combinations from different approximation spaces
    TPZManVector<REAL,3> u      = datavec[ublock].sol[0];
    REAL P              = datavec[Pblock].sol[0][0];
    
    TPZFMatrix<STATE> Graduaxes = datavec[ublock].dsol[0];
    TPZFMatrix<STATE> GradPaxes = datavec[Pblock].dsol[0];
    TPZFNMatrix<660,REAL> GradP;
    TPZAxesTools<REAL>::Axes2XYZ(GradPaxes, GradP, datavec[Pblock].axes);
    
    // Rock and fluids parameters
    TPZFMatrix<STATE> KInverse = fReservoirdata->KabsoluteInv();
    TPZFMatrix<STATE> Oneoverlambda_Kinv_u(2,1);
    TPZFMatrix<STATE> Gravity(2,1);
    TPZFMatrix<STATE> gm(2,1);
    
    TPZFMatrix<STATE> K;
    if (fSimulationData->IsHeterogeneousQ()) {
        int gelid = datavec[0].gelElId;
        K = fReservoirdata->GetKvector()[gelid];
        KInverse = fReservoirdata->ComputeInvKabsolute(K);
    }
    
    Gravity = fSimulationData->GetGravity();
    
    REAL phi, dphidP;
    this->fReservoirdata->Porosity(P, phi, dphidP);
    
    // time
    REAL dt = fSimulationData->GetDeltaT();
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    this->ComputeProperties(datavec, props);
    TPZManVector<REAL> rho_alpha        = props[0];
    TPZManVector<REAL> rho_beta         = props[1];
    TPZManVector<REAL> rho_gamma        = props[2];
    TPZManVector<REAL> f_alpha          = props[3];
    TPZManVector<REAL> f_beta           = props[4];
    TPZManVector<REAL> f_gamma          = props[5];
    TPZManVector<REAL> lambda           = props[6];
    TPZManVector<REAL> rho              = props[7];
    TPZManVector<REAL> rhof             = props[8];
    
    
    Oneoverlambda_Kinv_u(0,0) = (1.0/lambda[0])* (KInverse(0,0)*u[0] + KInverse(0,1)*u[1]);
    Oneoverlambda_Kinv_u(1,0) = (1.0/lambda[0])* (KInverse(1,0)*u[0] + KInverse(1,1)*u[1]);
    
    gm(0,0) = rhof[0] * Gravity(0,0);
    gm(1,0) = rhof[0] * Gravity(1,0);
    
    REAL divu = 0.0;
    TPZFMatrix<STATE> iphiuHdiv(2,1);
    int ishapeindex;
    int ivectorindex;
    
    /////////////////////////////////
    // Last State n
    if (fSimulationData->IsnStep()) {
        
        //  n state computations
        for (int ip = 0; ip < nphiPL2; ip++)
        {
            ef(ip + iniP) += 1.0 * weight * (1.0/dt) * phi * rho[0]*  phiPL2(ip,0);
        }
        
        return;
    }
    // Last State n
    /////////////////////////////////
    
    for (int iq = 0; iq < nphiuHdiv; iq++)
    {
        
        /* $ \underset{\Omega_{e}}{\int}\left(K\lambda\right)^{-1}\mathbf{q}\cdot\mathbf{v}\;\partial\Omega_{e}-\underset{\Omega_{e}}{\int}P\; div\left(\mathbf{v}\right)\partial\Omega-\underset{\Omega_{e}}{\int}\nabla\left(\rho_{f}g\; z\right)\cdot\mathbf{v} $ */
        
        ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
        
        iphiuHdiv(0,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(0,ivectorindex);
        iphiuHdiv(1,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(1,ivectorindex);
        
        ef(iq + iniu) += weight * ((Oneoverlambda_Kinv_u(0,0)*iphiuHdiv(0,0) + Oneoverlambda_Kinv_u(1,0)*iphiuHdiv(1,0))
                                   - P * DivergenceOnDeformed(iq,0)
                                   - (gm(0,0)*iphiuHdiv(0,0) + gm(1,0)*iphiuHdiv(1,0)) );
        
    }
    
    TPZManVector<STATE,1> fvalue(1,0.0);
    TPZFMatrix<REAL> Grad;
    if(fTimeDependentForcingFunction)
    {
        fTimeDependentForcingFunction->Execute(datavec[Pblock].x,fSimulationData->GetTime(),fvalue,Grad);
    }
    
    divu = (Graduaxes(0,0) + Graduaxes(1,1) + Graduaxes(2,2)); // Note::  use it for HdivPiola = 1 or constant Jacobian mappings
    
    /* $ - \underset{\Omega}{\int}w\; div\left(\mathbf{q}\right)\partial\Omega $ */
    for (int ip = 0; ip < nphiPL2; ip++)
    {
        ef(ip + iniP) += weight * (- divu + fvalue[0] - (1.0/dt) * phi * rho[0]) * phiPL2(ip,0);
        
    }
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::ContributeInterfaceDarcy(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2L        = datavecleft[Sablock].phi; // For L2  test functions Sw
    TPZFMatrix<STATE> dphiuH1L       = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L       = datavecleft[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Getting test and basis functions for the right side
    TPZFMatrix<REAL> phiuH1R         = datavecright[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2R         = datavecright[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2R        = datavecright[Sablock].phi; // For L2  test functions Sa
    TPZFMatrix<STATE> dphiuH1R       = datavecright[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2R       = datavecright[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSaL2L    = phiSaL2L.Rows();                                   // For L2   Sa
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSaL   = nphiPL2L       + iniPL;
    
    // Blocks dimensions and lengths for the right side
    int nphiuHdivR   = datavecright[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2R     = phiPL2R.Rows();                                    // For L2   P
    int nphiSaL2R    = phiSaL2R.Rows();                                   // For L2   Sa
    int iniuR    = 0;
    int iniPR    = nphiuHdivR     + iniuR;
    int iniSaR   = nphiPL2R       + iniPR;
    int iblock = iniSaL + nphiSaL2L;
    
    TPZManVector<REAL,3> n =  data.normal;
    
    // Getting linear combinations from different approximation spaces for the left side
    REAL P_L             = datavecleft[Pblock].sol[0][0];
    //    REAL Salpha_L        = datavecleft[Sablock].sol[0][0];
    
    // Getting linear combinations from different approximation spaces for the right side
    REAL P_R             = datavecright[Pblock].sol[0][0];
    //    REAL Salpha_R        = datavecright[Sablock].sol[0][0];
    
    TPZManVector<REAL,3> vi(2,0.0);
    int ivectorindex = 0;
    int ishapeindex =0;
    
    for (int iq = 0; iq < nphiuHdivL; iq++)
    {
        ivectorindex = datavecleft[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavecleft[ublock].fVecShapeIndex[iq].second;
        vi[0] = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
        vi[1] = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
        REAL vin = vi[0]*n[0] + vi[1]*n[1];
        
        ef(iq + iniuL) += weight * (P_R - P_L) * vin;
        
        for (int jp = 0; jp < nphiPL2L; jp++)
        {
            ek(iq + iniuL,jp + iniPL) += - 1.0 * weight * phiPL2L(jp,0) * vin;
        }
        
        for (int jp = 0; jp < nphiPL2R; jp++)
        {
            ek(iq + iniuL,iblock + jp + iniPR) += 1.0 * weight * phiPL2R(jp,0) * vin;
        }
        
    }
    
}

void TPZAxiSymmetricDarcyFlow::ContributeInterfaceDarcy(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
    
    //    // Full implicit case: there is no n state computations here
    //    if (fSimulationData->IsnStep()) {
    //        return;
    //    }
    
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2L        = datavecleft[Sablock].phi; // For L2  test functions Sw
    TPZFMatrix<STATE> dphiuH1L       = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L       = datavecleft[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Getting test and basis functions for the right side
    TPZFMatrix<REAL> phiuH1R         = datavecright[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2R         = datavecright[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2R        = datavecright[Sablock].phi; // For L2  test functions Sa
    TPZFMatrix<STATE> dphiuH1R       = datavecright[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2R       = datavecright[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSaL2L    = phiSaL2L.Rows();                                   // For L2   Sa
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSaL   = nphiPL2L       + iniPL;
    
    // Blocks dimensions and lengths for the right side
    int nphiuHdivR   = datavecright[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2R     = phiPL2R.Rows();                                    // For L2   P
    int nphiSaL2R    = phiSaL2R.Rows();                                   // For L2   Sa
    int iniuR    = 0;
    int iniPR    = nphiuHdivR     + iniuR;
//    int iniSaR   = nphiPL2R       + iniPR;
//    int iblock = iniSaL + nphiSaL2L;
    
    
    TPZManVector<REAL,3> n =  data.normal;
    
    // Getting linear combinations from different approximation spaces for the left side
    REAL P_L             = datavecleft[Pblock].sol[0][0];
//    REAL Salpha_L        = datavecleft[Sablock].sol[0][0];
    
    // Getting linear combinations from different approximation spaces for the right side
    REAL P_R             = datavecright[Pblock].sol[0][0];
//    REAL Salpha_R        = datavecright[Sablock].sol[0][0];
    
    TPZManVector<REAL,3> vi(2,0.0);
    int ivectorindex = 0;
    int ishapeindex =0;
    
    for (int iq = 0; iq < nphiuHdivL; iq++)
    {
        ivectorindex = datavecleft[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavecleft[ublock].fVecShapeIndex[iq].second;
        vi[0] = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
        vi[1] = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
        REAL vin = vi[0]*n[0] + vi[1]*n[1];
        
        ef(iq + iniuL) += weight * (P_R - P_L) * vin;
        
    }
    
}


void TPZAxiSymmetricDarcyFlow::ContributeBCDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // At each Integration Point.
    
    int ublock = 0;
    int Pblock = 1;
    
    
    // Getting test and basis functions
    TPZFNMatrix<9,STATE> phiuH1 = datavec[ublock].phi; // For H1   test functions
    TPZFNMatrix<9,STATE> WL2   = datavec[Pblock].phi; // For HL2  test functions
    
    TPZFMatrix<STATE> dPhiH1 = datavec[ublock].dphix; // Derivative For H1   test functions
    TPZFMatrix<STATE> dWL2   = datavec[Pblock].dphix; // Derivative For HL2  test functions
    
    // Getting Linear combinations of basis functions
    TPZManVector<STATE> u = datavec[ublock].sol[0];
    TPZManVector<STATE> P = datavec[Pblock].sol[0];
    
    // Computing normal flux
    REAL un = u[0];
    
    TPZFMatrix<STATE> dQdx = datavec[ublock].dsol[0];
    TPZFMatrix<STATE> dPdx = datavec[Pblock].dsol[0];
    
    // Number of phis
    int nPhiHdiv = phiuH1.Rows();  // For Hdiv
    //    int nPhiL2   = WL2.Rows();                                  // For L2
    
    TPZFMatrix<STATE> iPhiHdiv(2,1);
    TPZFMatrix<STATE> jPhiHdiv(2,1);
    
//    // Computing the radius
//    TPZFMatrix<REAL> x_spatial(3,1,0.0);
//    x_spatial(0,0) = datavec[0].x[0];
//    x_spatial(2,0) = datavec[0].x[2];
//    REAL r = Norm(x_spatial);
//    REAL s = 1.0;
    
//    if (fSimulationData->IsAxisymmetricQ()) {
//        s *= 2.0*M_PI*r;
//    }
//    
//    // Applying the scaling
//    phiuH1      *= 1.0/s;
//    un          *= 1.0/s;
    
    REAL Value;
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD inflow
        {
            
            Value = bc.Val2()(0,0);         //  Pressure
            TPZManVector<STATE,1> P_hydro(1,0.0);
            TPZFMatrix<STATE> GradP_hydro(1,0);            
            
            if (fSimulationData->IsHydrostaticBCQ()) {
                if (fTimedependentBCForcingFunction) {
                    fTimedependentBCForcingFunction->Execute(datavec[ublock].x, fSimulationData->GetTime(), P_hydro,GradP_hydro);
                    Value = P_hydro[0];
                }
                
                if (fSimulationData->IsTwoPhaseQ()) {
                    
                    REAL sw = datavec[2].sol[0][0];
                    
                    //    REAL x = datavec[ublock].x[0];
                    REAL y = datavec[ublock].x[1];
                    
                    //    REAL Kstr           = 1.0e-13;
                    REAL Pstr           = 2.0e7;
                    REAL Lstr           = 100.0;
                    //    REAL Mustr          = 0.001;
                    REAL Rhostr         = 1000.0;
                    REAL rho_alpha = 1000.0/Rhostr;
                    REAL rho_beta = 800.0/Rhostr;
                    REAL rho = sw * rho_alpha + (1-sw)*rho_beta;
                    REAL P_at_datum = bc.Val2()(0,0);//2.0*1.0e7;
                    REAL g = -10.0*((Lstr*Rhostr)/Pstr);
                    Value = (rho * g * y)+P_at_datum;
                }
            }
            
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                ef(iq) += weight * ( Value ) * phiuH1(iq,0);
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN inflow
        {
            
            Value = bc.Val2()(0,0);         //  NormalFlux
            
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                ef(iq) += weight * (gBigNumber * (un - Value)) * phiuH1(iq,0);
                
                for (int jq = 0; jq < nPhiHdiv; jq++)
                {
                    ek(iq,jq) += gBigNumber * weight * (phiuH1(jq,0)) * phiuH1(iq,0);
                }
                
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD outflow
        {
            
            Value = bc.Val2()(0,0);         //  Pressure
            TPZVec<STATE> P_hydro(1,0.0);
            TPZFMatrix<STATE> GradP_hydro(1,0);
            
            if (fSimulationData->IsHydrostaticBCQ()) {
                if (fTimedependentBCForcingFunction) {
                    fTimedependentBCForcingFunction->Execute(datavec[ublock].x, fSimulationData->GetTime(), P_hydro,GradP_hydro);
                    Value = P_hydro[0];
                }
                
                if (fSimulationData->IsTwoPhaseQ()) {
                    
                    REAL sw = datavec[2].sol[0][0];
                    
                    //    REAL x = datavec[ublock].x[0];
                    REAL y = datavec[ublock].x[1];
                    
                    //    REAL Kstr           = 1.0e-13;
                    REAL Pstr           = 2.0e7;
                    //    REAL Tstr           = 355.37;
                    //    REAL Tres           = 355.37;
                    REAL Lstr           = 100.0;
                    //    REAL Mustr          = 0.001;
                    REAL Rhostr         = 1000.0;
                    REAL rho_alpha = 1000.0/Rhostr;
                    REAL rho_beta = 800.0/Rhostr;
                    REAL rho = sw * rho_alpha + (1-sw)*rho_beta;
                    REAL P_at_datum = bc.Val2()(0,0);//2.0*1.0e7;
                    REAL g = -10.0*((Lstr*Rhostr)/Pstr);
                    Value = (rho * g * y)+P_at_datum;
                }
                
            }

            
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                
                ef(iq) += weight * ( Value ) * phiuH1(iq,0);
                
            }
        }
            break;
            
        case 3 :    // Neumann BC  QN outflow
        {
            Value = bc.Val2()(0,0);         //  NormalFlux
            
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                ef(iq) += weight * (gBigNumber * (un - Value) ) * phiuH1(iq,0);
                
                for (int jq = 0; jq < nPhiHdiv; jq++)
                {
                    ek(iq,jq) += gBigNumber * weight * (phiuH1(jq,0)) * phiuH1(iq,0);
                }
                
            }
            
        }
            break;
            
        case 4 :    // Neumann BC  QN impervious
        {
            
            
            Value = bc.Val2()(0,0);         //  NormalFlux
            
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                ef(iq) += weight * (gBigNumber * (un - Value) ) * phiuH1(iq,0);
                
                for (int jq = 0; jq < nPhiHdiv; jq++)
                {
                    ek(iq,jq) += gBigNumber * weight * (phiuH1(jq,0)) * phiuH1(iq,0);
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


void TPZAxiSymmetricDarcyFlow::ContributeAlpha(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // Getting data for Mixed-Darcy flow problem
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;         // Salpha Saturation needs L2 scalar functions     (phiPL2)
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2         = datavec[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> Sa_phiPL2         = datavec[Sablock].phi;  // For L2  test functions S_alpha
    
    TPZFMatrix<REAL> dphiuH1   = datavec[ublock].dphix;         // Derivative For H1  test functions u
    TPZFMatrix<REAL> dphiPL2   = datavec[Pblock].dphix;         // Derivative For L2  test functions P
    TPZFMatrix<REAL> dphiSL2   = datavec[Sablock].dphix;    // Derivative For L2  test functions S
    
    TPZFMatrix<STATE> DivergenceOnDeformed;
    // Compute the divergence on deformed element by piola contravariant transformation
    this->ComputeDivergenceOnDeformed(datavec, DivergenceOnDeformed);
    
    
    // Blocks dimensions and lengths
    int nphiuHdiv   = datavec[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2     = phiPL2.Rows();                                    // For L2   P
    int nphiSaL2     = Sa_phiPL2.Rows();                                    // For L2   S_alpha
    int iniu    = 0;
    int iniP    = nphiuHdiv     + iniu;
    int iniSa    = nphiPL2     + iniP;
    
    // Getting linear combinations from different approximation spaces
    TPZManVector<REAL,3> u      = datavec[ublock].sol[0];
    REAL P              = datavec[Pblock].sol[0][0];
    REAL S_alpha              = datavec[Sablock].sol[0][0];
    //
    //    TPZFMatrix<STATE> Graduaxes = datavec[ublock].dsol[0];
    TPZFMatrix<STATE> GradSaxes = datavec[Sablock].dsol[0];
    //    TPZFNMatrix<660,REAL> GradP;
    //    TPZAxesTools<REAL>::Axes2XYZ(GradPaxes, GradP, datavec[Pblock].axes);
    
    
    // Rock and fluids parameters
    TPZFMatrix<STATE> KInverse = fReservoirdata->KabsoluteInv();
    TPZFMatrix<STATE> K = fReservoirdata->Kabsolute();
    TPZFMatrix<STATE> Oneoverlambda_Kinv_u(2,1);
    TPZFMatrix<REAL> Oneoverlambda_Kinv_Phiu(2,1);
    TPZFMatrix<STATE> Gravity(2,1);
    TPZFMatrix<STATE> gm(2,1);
    TPZFMatrix<STATE> dgmdS(2,1);
    
    if (fSimulationData->IsHeterogeneousQ()) {
        int gelid = datavec[0].gelElId;
        K = fReservoirdata->GetKvector()[gelid];
        KInverse = fReservoirdata->ComputeInvKabsolute(K);
    }
    
    Gravity = fSimulationData->GetGravity();
    
    REAL phi, dphidP;
    this->fReservoirdata->Porosity(P, phi, dphidP);
    
    // time
    REAL dt = fSimulationData->GetDeltaT();
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    this->ComputeProperties(datavec, props);
    TPZManVector<REAL> rho_alpha        = props[0];
    TPZManVector<REAL> rho_beta         = props[1];
    TPZManVector<REAL> rho_gamma        = props[2];
    TPZManVector<REAL> f_alpha          = props[3];
    TPZManVector<REAL> f_beta           = props[4];
    TPZManVector<REAL> f_gamma          = props[5];
    TPZManVector<REAL> lambda           = props[6];
    TPZManVector<REAL> rho              = props[7];
    TPZManVector<REAL> rhof             = props[8];
    TPZManVector<REAL> Pc_beta_alpha    = props[9];
    
    Oneoverlambda_Kinv_u(0,0) = (1.0/lambda[0])* (KInverse(0,0)*u[0] + KInverse(0,1)*u[1]);
    Oneoverlambda_Kinv_u(1,0) = (1.0/lambda[0])* (KInverse(1,0)*u[0] + KInverse(1,1)*u[1]);
    
    gm(0,0) = rhof[0] * Gravity(0,0);
    gm(1,0) = rhof[0] * Gravity(1,0);
    
    dgmdS(0,0) = rhof[3] * Gravity(0,0);
    dgmdS(1,0) = rhof[3] * Gravity(1,0);
    
    TPZManVector<STATE> GradS(2,0.0);
    TPZManVector<STATE> Grad_Pc(2,0.0);
    
    //  Compute grad(S)
    GradS[0] = GradSaxes(0,0)*datavec[Sablock].axes(0,0)+GradSaxes(1,0)*datavec[Sablock].axes(1,0);
    GradS[1] = GradSaxes(0,0)*datavec[Sablock].axes(0,1)+GradSaxes(1,0)*datavec[Sablock].axes(1,1);
    
    //  Compute grad(Pc)
    Grad_Pc[0] = Pc_beta_alpha[3] * GradS[0];
    Grad_Pc[1] = Pc_beta_alpha[3] * GradS[1];
    
    TPZManVector<REAL,3> vi(2,0.0);
    TPZManVector<REAL,3> vj(2,0.0);
    int ishapeindex, jshapeindex;
    int ivectorindex, jvectorindex;
    
    
    for (int iq = 0; iq < nphiuHdiv; iq++)
    {
        
        /* $ \underset{\Omega_{e}}{\int}\left(K\lambda\right)^{-1}\mathbf{q}\cdot\mathbf{v}\;\partial\Omega_{e}-\underset{\Omega_{e}}{\int}P\; div\left(\mathbf{v}\right)\partial\Omega-\underset{\Omega_{e}}{\int}\nabla\left(\rho_{f}g\; z\right)\cdot\mathbf{v} $ */
        
        ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
        
        vi[0] = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(0,ivectorindex);
        vi[1] = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(1,ivectorindex);
        
        ef(iq + iniu) += weight * (-1.0 * (S_alpha * Pc_beta_alpha[0]) * DivergenceOnDeformed(iq,0) - f_alpha[0] * (Grad_Pc[0]*vi[0] + Grad_Pc[1]*vi[1]));
        
        
        // du/dalphaSa terms
        for (int jsa = 0; jsa < nphiSaL2; jsa++)
        {
            
            //  Compute grad(phi)
            TPZManVector<STATE> Gradphis(2,0);
            Gradphis[0] = dphiSL2(0,jsa)*datavec[Sablock].axes(0,0)+dphiSL2(1,jsa)*datavec[Sablock].axes(1,0);
            Gradphis[1] = dphiSL2(0,jsa)*datavec[Sablock].axes(0,1)+dphiSL2(1,jsa)*datavec[Sablock].axes(1,1);
            
            ek(iq + iniu,jsa + iniSa) += weight * ( (-lambda[3]/lambda[0])*(Oneoverlambda_Kinv_u(0,0)*vi[0] + Oneoverlambda_Kinv_u(1,0)*vi[1]) - (dgmdS(0,0)*vi[0] + dgmdS(1,0)*vi[1])
                                                   - ( S_alpha * Pc_beta_alpha[3] + 1.0 * Pc_beta_alpha[0] ) * DivergenceOnDeformed(iq,0)
                                                   - ( f_alpha[3]*(Grad_Pc[0]*vi[0] + Grad_Pc[1]*vi[1]) + (f_alpha[0]*Pc_beta_alpha[3])*(Gradphis[0]*vi[0] + Gradphis[1]*vi[1])) )* Sa_phiPL2(jsa,0) ;
            
//            ek(iq + iniu,jsa + iniSa) += weight * ( (-lambda[3]/lambda[0])*(Oneoverlambda_Kinv_u(0,0)*vi[0] + Oneoverlambda_Kinv_u(1,0)*vi[1]) - (dgmdS(0,0)*vi[0] + dgmdS(1,0)*vi[1])
//                                                   - ( S_alpha * Pc_beta_alpha[3] + 1.0 * Pc_beta_alpha[0] ) * DivergenceOnDeformed(iq,0)
//                                                   - ( f_alpha[3]*(Grad_Pc[0]*vi[0] + Grad_Pc[1]*vi[1]) ) )* Sa_phiPL2(jsa,0) ;
            
        }
        
    }
    
    /* $ - \underset{\Omega}{\int}w\; div\left(\mathbf{q}\right)\partial\Omega $ */
    for (int ip = 0; ip < nphiPL2; ip++)
    {
        for (int jsa = 0; jsa < nphiSaL2; jsa++)
        {
            ek(ip + iniP, jsa + iniSa) += weight * (- (1.0/dt) * phi * rho[3] * Sa_phiPL2(jsa,0) ) * phiPL2(ip,0);
        }
    }
    
    
//    // Gravitational Segregation
//    TPZVec< TPZManVector<REAL,3> > qg;
//    this->GravitationalSegregation(datavec,qg);
//    
//    // Capillary Segregation
//    TPZVec< TPZManVector<REAL,3> > qc;
//    TPZManVector<REAL,3> grad_s(3,0.0);
//    this->CapillarySegregation(datavec,qc,grad_s);
    
    for (int isw = 0; isw < nphiSaL2; isw++)
    {
        
//        //  Compute grad(phi)
//        TPZManVector<STATE> Gradphis(2,0);
//        Gradphis[0] = dphiSL2(0,isw)*datavec[Sablock].axes(0,0)+dphiSL2(1,isw)*datavec[Sablock].axes(1,0);
//        Gradphis[1] = dphiSL2(0,isw)*datavec[Sablock].axes(0,1)+dphiSL2(1,isw)*datavec[Sablock].axes(1,1);
//        
//        ef(isw  + iniSa ) += weight * ( (1.0/dt) * phi * S_alpha * rho_alpha[0] * Sa_phiPL2(isw,0)
//                                       - f_alpha[0]*(u[0]*Gradphis[0] + u[1]*Gradphis[1])
//                                       - (qg[0][0]*Gradphis[0] + qg[0][1]*Gradphis[1])
//                                       - (qc[0][0]*Gradphis[0] + qc[0][1]*Gradphis[1]));
        
        ef(isw  + iniSa ) += weight * ( (1.0/dt) * phi * S_alpha * rho_alpha[0] * Sa_phiPL2(isw,0) );

        
//        // du/dalphau terms
//        for (int jq = 0; jq < nphiuHdiv; jq++)
//        {
//            jvectorindex = datavec[ublock].fVecShapeIndex[jq].first;
//            jshapeindex = datavec[ublock].fVecShapeIndex[jq].second;
//            
//            vi[0] = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(0,jvectorindex);
//            vi[1] = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(1,jvectorindex);
//            
//            ek(isw  + iniSa, jq  + iniu ) = -1.0 * weight * (f_alpha[0]*(vi[0]*Gradphis[0] + vi[1]*Gradphis[1]));
//            
//        }
        
        for (int jp = 0; jp < nphiPL2; jp++)
        {
//            ek(isw  + iniSa, jp  + iniP ) += 1.0 * weight * ( (1.0/dt) * phi * S_alpha * rho_alpha[2] * Sa_phiPL2(isw,0)
//                                                             - f_alpha[2]*(u[0]*Gradphis[0] + u[1]*Gradphis[1])
//                                                             - (qg[2][0]*Gradphis[0] + qg[2][1]*Gradphis[1])
//                                                             - (qc[2][0]*Gradphis[0] + qc[2][1]*Gradphis[1]) ) * phiPL2(jp,0);
            
            ek(isw  + iniSa, jp  + iniP ) += 1.0 * weight * ( (1.0/dt) * phi * S_alpha * rho_alpha[2] * Sa_phiPL2(isw,0) ) * phiPL2(jp,0);
            
            
        }
        
        for (int jsw = 0; jsw < nphiSaL2; jsw++)
        {
            
//            //  Compute grad(phi)
//            TPZManVector<STATE> Gradphis_j(2,0);
//            Gradphis_j[0] = dphiSL2(0,jsw)*datavec[Sablock].axes(0,0)+dphiSL2(1,jsw)*datavec[Sablock].axes(1,0);
//            Gradphis_j[1] = dphiSL2(0,jsw)*datavec[Sablock].axes(0,1)+dphiSL2(1,jsw)*datavec[Sablock].axes(1,1);
//            
//            ek(isw  + iniSa, jsw  + iniSa ) += 1.0 * weight * ((1.0/dt) * phi * rho_alpha[0] * Sa_phiPL2(isw,0)
//                                                               - f_alpha[3]*(u[0]*Gradphis[0] + u[1]*Gradphis[1])
//                                                               - (qg[3][0]*Gradphis[0] + qg[3][1]*Gradphis[1])
//                                                               - (qc[3][0]*Gradphis[0] + qc[3][1]*Gradphis[1])
//                                                               - ( (qc[4][0]) * (K(0,0)*Gradphis_j[0] + K(0,1)*Gradphis_j[1]) *Gradphis[0]
//                                                                  + (qc[4][1]) * (K(1,0)*Gradphis_j[0] + K(1,1)*Gradphis_j[1]) *Gradphis[1] ) ) * Sa_phiPL2(jsw,0) ;
            
            ek(isw  + iniSa, jsw  + iniSa ) += 1.0 * weight * ((1.0/dt) * phi * rho_alpha[0] * Sa_phiPL2(isw,0)) * Sa_phiPL2(jsw,0) ;
            
        }
    }
    
}


void TPZAxiSymmetricDarcyFlow::ContributeAlpha(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    
    // Getting data for Mixed-Darcy flow problem
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;         // Salpha Saturation needs L2 scalar functions     (phiPL2)
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2         = datavec[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> Sa_phiPL2         = datavec[Sablock].phi;  // For L2  test functions S_alpha
    
    TPZFMatrix<REAL> dphiuH1   = datavec[ublock].dphix;         // Derivative For H1  test functions u
    TPZFMatrix<REAL> dphiPL2   = datavec[Pblock].dphix;         // Derivative For L2  test functions P
    TPZFMatrix<REAL> dphiSL2   = datavec[Sablock].dphix;    // Derivative For L2  test functions S
    
    TPZFMatrix<STATE> DivergenceOnDeformed;
    // Compute the divergence on deformed element by piola contravariant transformation
    this->ComputeDivergenceOnDeformed(datavec, DivergenceOnDeformed);
    
    
    // Blocks dimensions and lengths
    int nphiuHdiv   = datavec[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2     = phiPL2.Rows();                                    // For L2   P
    int nphiSaL2     = Sa_phiPL2.Rows();                                    // For L2   S_alpha
    int iniu    = 0;
    int iniP    = nphiuHdiv     + iniu;
    int iniSa    = nphiPL2     + iniP;
    
    // Getting linear combinations from different approximation spaces
    TPZManVector<REAL,3> u      = datavec[ublock].sol[0];
    REAL P              = datavec[Pblock].sol[0][0];
    REAL S_alpha              = datavec[Sablock].sol[0][0];

    TPZFMatrix<STATE> GradSaxes = datavec[Sablock].dsol[0];
////    TPZFMatrix<STATE> GradPaxes = datavec[Pblock].dsol[0];
//    TPZFNMatrix<660,REAL> GradS;
//    TPZAxesTools<REAL>::Axes2XYZ(GradSaxes, GradS, datavec[Pblock].axes);
    
    
    
    // Rock and fluids parameters
    TPZFMatrix<STATE> KInverse = fReservoirdata->KabsoluteInv();
    TPZFMatrix<STATE> Oneoverlambda_Kinv_u(2,1);
    TPZFMatrix<REAL> Oneoverlambda_Kinv_Phiu(2,1);
    TPZFMatrix<STATE> Gravity(2,1);
    TPZFMatrix<STATE> gm(2,1);
    TPZFMatrix<STATE> dgmdP(2,1);
    
    TPZFMatrix<STATE> K;
    if (fSimulationData->IsHeterogeneousQ()) {
        int gelid = datavec[0].gelElId;
        K = fReservoirdata->GetKvector()[gelid];
        KInverse = fReservoirdata->ComputeInvKabsolute(K);
    }
    
    Gravity = fSimulationData->GetGravity();
    
    REAL phi, dphidP;
    this->fReservoirdata->Porosity(P, phi, dphidP);
    
    // time
    REAL dt = fSimulationData->GetDeltaT();
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    this->ComputeProperties(datavec, props);
    TPZManVector<REAL> rho_alpha        = props[0];
    TPZManVector<REAL> rho_beta         = props[1];
    TPZManVector<REAL> rho_gamma        = props[2];
    TPZManVector<REAL> f_alpha          = props[3];
    TPZManVector<REAL> f_beta           = props[4];
    TPZManVector<REAL> f_gamma          = props[5];
    TPZManVector<REAL> lambda           = props[6];
    TPZManVector<REAL> rho              = props[7];
    TPZManVector<REAL> rhof             = props[8];
    TPZManVector<REAL> Pc_beta_alpha    = props[9];
    
    /////////////////////////////////
    // Last State n
    if (fSimulationData->IsnStep()) {
        
        if(fSimulationData->IsImpesQ()){
            
            for (int isw = 0; isw < nphiSaL2; isw++)
            {
                //  Compute grad(phi)
                TPZManVector<STATE> Gradphis(2,0);
                Gradphis[0] = dphiSL2(0,isw)*datavec[Sablock].axes(0,0)+dphiSL2(1,isw)*datavec[Sablock].axes(1,0);
                Gradphis[1] = dphiSL2(0,isw)*datavec[Sablock].axes(0,1)+dphiSL2(1,isw)*datavec[Sablock].axes(1,1);
                
                ef(isw  + iniSa ) +=  weight * ( (-1.0/dt) * phi * S_alpha * rho_alpha[0]  * Sa_phiPL2(isw,0) - f_alpha[0]*(u[0]*Gradphis[0] + u[1]*Gradphis[1]));
            }
        }
        else{
            for (int isw = 0; isw < nphiSaL2; isw++)
            {
                ef(isw  + iniSa ) += - weight * (1.0/dt) * phi * S_alpha * rho_alpha[0]  * Sa_phiPL2(isw,0);
            }
            
        }
        
        
        return;
    }
    // Last State n
    /////////////////////////////////


    TPZManVector<STATE> GradS(2,0.0);
    TPZManVector<STATE> Grad_Pc(2,0.0);
    
    //  Compute grad(S)
    GradS[0] = GradSaxes(0,0)*datavec[Sablock].axes(0,0)+GradSaxes(1,0)*datavec[Sablock].axes(1,0);
    GradS[1] = GradSaxes(0,0)*datavec[Sablock].axes(0,1)+GradSaxes(1,0)*datavec[Sablock].axes(1,1);
    
    //  Compute grad(Pc)
    Grad_Pc[0] = Pc_beta_alpha[3] * GradS[0];
    Grad_Pc[1] = Pc_beta_alpha[3] * GradS[1];
    
    TPZManVector<REAL,3> vi(2,0);
    int ishapeindex;
    int ivectorindex;
    
    for (int iq = 0; iq < nphiuHdiv; iq++)
    {
        ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
        
        vi[0] = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(0,ivectorindex);
        vi[1] = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(1,ivectorindex);
            
        ef(iq + iniu) += weight * (-1.0 * (S_alpha * Pc_beta_alpha[0]) * DivergenceOnDeformed(iq,0) - f_alpha[0] * (Grad_Pc[0]*vi[0] + Grad_Pc[1]*vi[1]));
        
    }
    
//    // Gravitational Segregation
//    TPZVec< TPZManVector<REAL,3> > qg;
//    this->GravitationalSegregation(datavec,qg);
//
//    // Capillary Segregation
//    TPZVec< TPZManVector<REAL,3> > qc;
//    TPZManVector<REAL,3> grad_s(3,0.0);
//    this->CapillarySegregation(datavec,qc,grad_s);
    
    for (int isw = 0; isw < nphiSaL2; isw++)
    {
//        //  Compute grad(phi)
//        TPZManVector<STATE> Gradphis(2,0);
//        Gradphis[0] = dphiSL2(0,isw)*datavec[Sablock].axes(0,0)+dphiSL2(1,isw)*datavec[Sablock].axes(1,0);
//        Gradphis[1] = dphiSL2(0,isw)*datavec[Sablock].axes(0,1)+dphiSL2(1,isw)*datavec[Sablock].axes(1,1);
//        
//        ef(isw  + iniSa ) += weight * ( (1.0/dt) * phi * S_alpha * rho_alpha[0] * Sa_phiPL2(isw,0)
//                                       - f_alpha[0]*(u[0]*Gradphis[0] + u[1]*Gradphis[1])
//                                       - (qg[0][0]*Gradphis[0] + qg[0][1]*Gradphis[1])
//                                       - (qc[0][0]*Gradphis[0] + qc[0][1]*Gradphis[1]));
        
        ef(isw  + iniSa ) += weight * ( (1.0/dt) * phi * S_alpha * rho_alpha[0] * Sa_phiPL2(isw,0));
        
    }
    
    
}


void TPZAxiSymmetricDarcyFlow::ContributeBCInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }

    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2L        = datavecleft[Sablock].phi; // For L2  test functions Sw
    TPZFMatrix<STATE> dphiuH1L   = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L   = datavecleft[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSaL2L    = phiSaL2L.Rows();                                   // For L2   Sa
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSaL   = nphiPL2L       + iniPL;
    
    
    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
//    REAL PL              = datavecleft[Pblock].sol[0][0];
    REAL Salpha_L             = datavecleft[Sablock].sol[0][0];
    
    TPZFMatrix<STATE> iphiuHdivL(2,1);
    int ishapeindex;
    int ivectorindex;
    
    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = datavecleft[0].x[0];
    REAL r = Norm(x_spatial);
    REAL s = 1.0;
    
//    if (fSimulationData->IsAxisymmetricQ()) {
//        s *= 2.0*M_PI*r;
//    }
//
//    // Applying the scaling
//    phiuH1L      *= 1.0/s;
//    //    Graduaxes   *= 1.0/s;
//    //    DivergenceOnDeformed *= 1.0/s;
//    uL[0] *= 1.0/s;
//    uL[1] *= 1.0/s;
//    uL[2] *= 1.0/s;
    
    
    
    TPZManVector<REAL,3> n =  data.normal;
    REAL uLn = uL[0]*n[0] + uL[1]*n[1] + uL[2]*n[2];
    
    
    // Capillary Pressure Pc_beta_beta
    TPZVec< TPZManVector<REAL>  > propsL;
    TPZVec< TPZManVector<REAL>  > propsR;
    
    ComputeProperties(datavecleft, propsL);
    TPZManVector<REAL> Pc_beta_alphaL          = propsL[9];
    
    
    
    TPZManVector<REAL,3> vi(2,0.0);
    
    for (int iq = 0; iq < nphiuHdivL; iq++)
    {
        ivectorindex = datavecleft[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavecleft[ublock].fVecShapeIndex[iq].second;
        vi[0] = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
        vi[1] = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
        REAL vin = vi[0]*n[0] + vi[1]*n[1];
        
        ef(iq + iniuL) += -1.0*weight * (Salpha_L * Pc_beta_alphaL[0] ) * vin;
        
        for (int jsw = 0; jsw < nphiSaL2L; jsw++)
        {
            ek(iq + iniuL, jsw + iniSaL) += -1.0 * weight * (Pc_beta_alphaL[0]+Salpha_L * Pc_beta_alphaL[3]) * phiSaL2L(jsw,0) * vin;
        }
    }
    
//    // Capillary Segregation
//    TPZVec< TPZManVector<REAL,3> > qc_L;
//    this->CapillarySegregation(datavecleft,qc_L);
//    
//    // Comuting the contribution with the minimum gravitational flux
//    
//    for (int isw = 0; isw < nphiSaL2L; isw++)
//    {
//        ef(isw + iniSaL) += 1.0 * weight * (qc_L[0][0]*n[0] + qc_L[0][1]*n[1]) * phiSaL2L(isw,0);
//    }
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD inflow
        {
            
            REAL PD         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon) {
                //                std::cout << "Outflow condition detected in inflow condition qn = " << uLn << std::endl;
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                    for (int jq = 0; jq < nphiuHdivL; jq++)
                    {
                        ivectorindex = datavecleft[ublock].fVecShapeIndex[jq].first;
                        ishapeindex = datavecleft[ublock].fVecShapeIndex[jq].second;
                        iphiuHdivL(0,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
                        iphiuHdivL(1,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
                        REAL vn = iphiuHdivL(0,0)*n[0] + iphiuHdivL(1,0)*n[1];
                        
                        ek(isw + iniSaL, jq + iniuL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * vn;
                    }
                    
                    for (int jp = 0; jp < nphiPL2L; jp++)
                    {
                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * uLn;
                    }
                    
                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
                    {
                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * uLn;
                    }
                }
            }
            else
            {
                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
                state_vars[0] = uLn;                            //  qn
                state_vars[1] = PD;                             //  P
                state_vars[2] = S_alpha;                        //  S_alpha
                state_vars[3] = S_beta;                         //  S_beta
                
                this->ComputeProperties(state_vars, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
//                    for (int jp = 0; jp < nphiPL2L; jp++)
//                    {
//                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * uLn;
//                    }
//
//                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
//                    {
//                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * uLn;
//                    }
                }
                
            }
            
        }
            break;
            
        case 1 :    // Neumann BC  QN inflow
        {
            REAL qn         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon ) {
                std::cout << "N:: Outflow condition detected in inflow condition qn = " << uLn << "; qN = " << qn << std::endl;
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                    for (int jq = 0; jq < nphiuHdivL; jq++)
                    {
                        ivectorindex = datavecleft[ublock].fVecShapeIndex[jq].first;
                        ishapeindex = datavecleft[ublock].fVecShapeIndex[jq].second;
                        iphiuHdivL(0,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
                        iphiuHdivL(1,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
                        REAL vn = iphiuHdivL(0,0)*n[0] + iphiuHdivL(1,0)*n[1];
                        
                        ek(isw + iniSaL, jq + iniuL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * vn;
                    }
                    
                    for (int jp = 0; jp < nphiPL2L; jp++)
                    {
                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * uLn;
                    }
                    
                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
                    {
                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * uLn;
                    }
                }
            }
            else
            {
                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
                state_vars[0] = qn;                             //  qn
                state_vars[1] = datavecleft[Pblock].sol[0][0];  //  P
                state_vars[2] = S_alpha;                        //  S_alpha
                state_vars[3] = S_beta;                         //  S_beta
                
                this->ComputeProperties(state_vars, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * qn;
                    
//                    for (int jq = 0; jq < nphiuHdivL; jq++)
//                    {
//                        ivectorindex = datavecleft[ublock].fVecShapeIndex[jq].first;
//                        ishapeindex = datavecleft[ublock].fVecShapeIndex[jq].second;
//                        iphiuHdivL(0,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
//                        iphiuHdivL(1,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
//                        REAL vn = iphiuHdivL(0,0)*n[0] + iphiuHdivL(1,0)*n[1];
//                        
//                        ek(isw + iniSaL, jq + iniuL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * vn;
//                    }
//                    
//                    for (int jp = 0; jp < nphiPL2L; jp++)
//                    {
//                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * qn;
//                    }
//
//                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
//                    {
//                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * qn;
//                    }
                    
                }
                
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD outflow
        {
//            REAL PD         = bc.Val2()(0,0);
//            REAL S_alpha    = bc.Val2()(1,0);
//            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon) {
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
//                    for (int jq = 0; jq < nphiuHdivL; jq++)
//                    {
//                        ivectorindex = datavecleft[ublock].fVecShapeIndex[jq].first;
//                        ishapeindex = datavecleft[ublock].fVecShapeIndex[jq].second;
//                        iphiuHdivL(0,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
//                        iphiuHdivL(1,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
//                        REAL vn = iphiuHdivL(0,0)*n[0] + iphiuHdivL(1,0)*n[1];
//                        ek(isw + iniSaL, jq + iniuL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * vn;
//                    }
//                    
//                    for (int jp = 0; jp < nphiPL2L; jp++)
//                    {
//                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * uLn;
//                    }
//                    
//                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
//                    {
//                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * uLn;
//                    }
                }
            }
            else
            {
                std::cout << "D:: Inflow  condition detected in outflow condition qn = " <<  uLn << std::endl;
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                    for (int jq = 0; jq < nphiuHdivL; jq++)
                    {
                        ivectorindex = datavecleft[ublock].fVecShapeIndex[jq].first;
                        ishapeindex = datavecleft[ublock].fVecShapeIndex[jq].second;
                        iphiuHdivL(0,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
                        iphiuHdivL(1,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
                        REAL vn = iphiuHdivL(0,0)*n[0] + iphiuHdivL(1,0)*n[1];
                        
                        ek(isw + iniSaL, jq + iniuL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * vn;
                    }
                    
                    for (int jp = 0; jp < nphiPL2L; jp++)
                    {
                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * uLn;
                    }
                    
                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
                    {
                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * uLn;
                    }
                }
                
            }
            
        }
            
        break;
            
        case 3 :    // Neumann BC  QN outflow
        {
            REAL qn         = bc.Val2()(0,0);
//            REAL S_alpha    = bc.Val2()(1,0);
//            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon ) {
                
//                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
//                state_vars[0] = qn;                             //  qn
//                state_vars[1] = datavecleft[Pblock].sol[0][0];  //  P
//                state_vars[2] = S_alpha;                        //  S_alpha
//                state_vars[3] = S_beta;                         //  S_beta
                
//                this->ComputeProperties(state_vars, props);
//                TPZManVector<REAL> f_alpha          = props[3];
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * qn/s;
                    
                    for (int jq = 0; jq < nphiuHdivL; jq++)
                    {
                        ivectorindex = datavecleft[ublock].fVecShapeIndex[jq].first;
                        ishapeindex = datavecleft[ublock].fVecShapeIndex[jq].second;
                        iphiuHdivL(0,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
                        iphiuHdivL(1,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
                        REAL vn = iphiuHdivL(0,0)*n[0] + iphiuHdivL(1,0)*n[1];
                        
                        ek(isw + iniSaL, jq + iniuL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * vn;
                    }
                    
                    for (int jp = 0; jp < nphiPL2L; jp++)
                    {
                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * qn;
                    }
                    
                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
                    {
                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * qn;
                    }
                }
            }
            else
            {
                std::cout << "Inflow  condition detected in outflow condition qn = " << uLn << std::endl;
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                    for (int jq = 0; jq < nphiuHdivL; jq++)
                    {
                        ivectorindex = datavecleft[ublock].fVecShapeIndex[jq].first;
                        ishapeindex = datavecleft[ublock].fVecShapeIndex[jq].second;
                        iphiuHdivL(0,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
                        iphiuHdivL(1,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
                        REAL vn = iphiuHdivL(0,0)*n[0] + iphiuHdivL(1,0)*n[1];
                        
                        ek(isw + iniSaL, jq + iniuL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * vn;
                    }
                    
                    for (int jp = 0; jp < nphiPL2L; jp++)
                    {
                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * uLn;
                    }
                    
                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
                    {
                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * uLn;
                    }
                }
                
            }
            
        }
            break;
            
        case 4 :    // Neumann BC QN impervious
        {
            REAL qn         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
            state_vars[0] = qn;                             //  qn
            state_vars[1] = datavecleft[Pblock].sol[0][0];  //  P
            state_vars[2] = S_alpha;                        //  S_alpha
            state_vars[3] = S_beta;                         //  S_beta
            
            this->ComputeProperties(state_vars, props);
            TPZManVector<REAL> f_alpha          = props[3];
            
            for (int isw = 0; isw < nphiSaL2L; isw++)
            {
                ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * qn;
                
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


void TPZAxiSymmetricDarcyFlow::ContributeBCInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
//    // Full implicit case: there is no n state computations here
//    if (fSimulationData->IsnStep()) {
//        return;
//    }

    // Getting data from different approximation spaces

    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)

    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2L        = datavecleft[Sablock].phi; // For L2  test functions Sw
    TPZFMatrix<STATE> dphiuH1L   = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L   = datavecleft[Pblock].dphix; // Derivative For L2  test functions


    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSaL2L    = phiSaL2L.Rows();                                   // For L2   Sa
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSaL   = nphiPL2L       + iniPL;


    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
    REAL P_L              = datavecleft[Pblock].sol[0][0];
    REAL Salpha_L             = datavecleft[Sablock].sol[0][0];

    // Computing the radius
    TPZFMatrix<REAL> x_spatial(3,1,0.0);
    x_spatial(0,0) = datavecleft[0].x[0];
    REAL r = Norm(x_spatial);
    REAL s = 1.0;
    
    if (fSimulationData->IsAxisymmetricQ()) {
        s *= 2.0*M_PI*r;
    }
    
//    // Applying the scaling
//    phiuH1L      *= 1.0/s;
//    //    Graduaxes   *= 1.0/s;
//    //    DivergenceOnDeformed *= 1.0/s;
//    uL[0] *= 1.0/s;
//    uL[1] *= 1.0/s;
//    uL[2] *= 1.0/s;

    TPZManVector<REAL,3> n =  data.normal;
    REAL uLn = uL[0]*n[0] + uL[1]*n[1] + uL[2]*n[2];

    // Capillary Pressure Pc_beta_beta
    TPZVec< TPZManVector<REAL>  > propsL;
    TPZVec< TPZManVector<REAL>  > propsR;
    
    ComputeProperties(datavecleft, propsL);
    TPZManVector<REAL> Pc_beta_alphaL          = propsL[9];
    
    
    TPZManVector<REAL,3> vi(2,0.0);
    int ivectorindex = 0;
    int ishapeindex =0;
    
    for (int iq = 0; iq < nphiuHdivL; iq++)
    {
        ivectorindex = datavecleft[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavecleft[ublock].fVecShapeIndex[iq].second;
        vi[0] = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
        vi[1] = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
        REAL vin = vi[0]*n[0] + vi[1]*n[1];
        
        ef(iq + iniuL) += -1.0*weight * (Salpha_L * Pc_beta_alphaL[0]) * vin;
        
    }
    
//    // Capillary Segregation
//    TPZVec< TPZManVector<REAL,3> > qc_L;
//    this->CapillarySegregation(datavecleft,qc_L);
//    
//    // Comuting the contribution with the minimum gravitational flux
//    
//    for (int isw = 0; isw < nphiSaL2L; isw++)
//    {
//        ef(isw + iniSaL) += 1.0 * weight * (qc_L[0][0]*n[0] + qc_L[0][1]*n[1]) * phiSaL2L(isw,0);
//    }

    

    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD inflow
        {
            
            REAL PD         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon) {
                //                std::cout << "Outflow condition detected in inflow condition qn = " << uLn << std::endl;
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                }
            }
            else
            {
                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
                state_vars[0] = uLn;                            //  qn
                state_vars[1] = PD;                             //  P
                state_vars[2] = S_alpha;                        //  S_alpha
                state_vars[3] = S_beta;                         //  S_beta
                
                this->ComputeProperties(state_vars, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                }
                
            }
            
        }
            break;
            
        case 1 :    // Neumann BC  QN inflow
        {
            REAL qn         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon ) {
                //                std::cout << "Outflow condition detected in inflow condition qn = " << uLn << std::endl;
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                }
            }
            else
            {
                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
                state_vars[0] = qn;                             //  qn
                state_vars[1] = datavecleft[Pblock].sol[0][0];  //  P
                state_vars[2] = S_alpha;                        //  S_alpha
                state_vars[3] = S_beta;                         //  S_beta
                
                this->ComputeProperties(state_vars, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * qn/s;
                    
                }
                
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD outflow
        {
            REAL PD         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon) {
                
                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
                state_vars[0] = uLn;                            //  qn
                state_vars[1] = PD;                             //  P
                state_vars[2] = S_alpha;                        //  S_alpha
                state_vars[3] = S_beta;                         //  S_beta
                
                //                this->ComputeProperties(state_vars, props);
                //                TPZManVector<REAL> f_alpha          = props[3];
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                }
            }
            else
            {
                std::cout << "Inflow  condition detected in outflow condition qn = " <<  uLn << std::endl;
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                }
                
            }
            
        }
            
            break;
            
        case 3 :    // Neumann BC  QN outflow
        {
            REAL qn         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon ) {
                
                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
                state_vars[0] = qn;                             //  qn
                state_vars[1] = datavecleft[Pblock].sol[0][0];  //  P
                state_vars[2] = S_alpha;                        //  S_alpha
                state_vars[3] = S_beta;                         //  S_beta
                
                this->ComputeProperties(state_vars, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                //                this->ComputeProperties(datavecleft, props);
                //                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * qn/s;
                    
                }
            }
            else
            {
                std::cout << "Inflow  condition detected in outflow condition qn = " << uLn << std::endl;
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                }
                
            }
            
        }
            break;
            
        case 4 :    // Neumann BC QN impervious
        {
            REAL qn         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
            state_vars[0] = qn;                             //  qn
            state_vars[1] = datavecleft[Pblock].sol[0][0];  //  P
            state_vars[2] = S_alpha;                        //  S_alpha
            state_vars[3] = S_beta;                         //  S_beta
            
            this->ComputeProperties(state_vars, props);
            TPZManVector<REAL> f_alpha          = props[3];
            
            for (int isw = 0; isw < nphiSaL2L; isw++)
            {
                ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * qn/s;
                
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
 * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZAxiSymmetricDarcyFlow::ContributeInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }

    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2L        = datavecleft[Sablock].phi; // For L2  test functions Sw
    TPZFMatrix<STATE> dphiuH1L   = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L   = datavecleft[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Getting test and basis functions for the right side
    TPZFMatrix<REAL> phiuH1R         = datavecright[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2R         = datavecright[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2R        = datavecright[Sablock].phi; // For L2  test functions Sa
    TPZFMatrix<STATE> dphiuH1R   = datavecright[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2R   = datavecright[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSaL2L    = phiSaL2L.Rows();                                   // For L2   Sa
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSaL   = nphiPL2L       + iniPL;
    
    // Blocks dimensions and lengths for the right side
    int nphiuHdivR   = datavecright[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2R     = phiPL2R.Rows();                                    // For L2   P
    int nphiSaL2R    = phiSaL2R.Rows();                                   // For L2   Sa
    int iniuR    = 0;
    int iniPR    = nphiuHdivR     + iniuR;
    int iniSaR   = nphiPL2R       + iniPR;
    
    int iblock = iniSaL + nphiSaL2L;
    int iblockt = 0;
    
    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
//    REAL P_L             = datavecleft[Pblock].sol[0][0];
    REAL Salpha_L        = datavecleft[Sablock].sol[0][0];
    
    // Getting linear combinations from different approximation spaces for the right side
    TPZManVector<REAL,3> uR      = datavecright[ublock].sol[0];
//    REAL P_R             = datavecright[Pblock].sol[0][0];
    REAL Salpha_R        = datavecright[Sablock].sol[0][0];
    
    
    TPZManVector<REAL,3> n =  data.normal;
    REAL uLn = uL[0]*n[0] + uL[1]*n[1] + uL[2]*n[2];
    
    
    // Getting test and basis functions for the left side
    TPZVec<TPZMaterialData> datavec;
    TPZFMatrix<REAL> phiuH1;
    TPZFMatrix<REAL> phiPL2;
    TPZFMatrix<REAL> phiSaL2;
    TPZFMatrix<STATE> dphiuH1;
    TPZFMatrix<STATE> dphiPL2;
    int nphiuHdiv   = 0;
    int nphiPL2     = 0;
    int nphiSaL2    = 0;
    int iniu    = 0;
    int iniP    = 0;
    int iniSa   = 0;
    
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    
    if ((uLn >= 0.0) && fabs(uLn) > fepsilon) {
        this->ComputeProperties(datavecleft, props);
        datavec = datavecleft;
        phiuH1      = phiuH1L;
        phiPL2      = phiPL2L;
        phiSaL2     = phiSaL2L;
        dphiuH1     = dphiuH1L;
        dphiPL2     = dphiPL2L;
        nphiuHdiv   = nphiuHdivL;
        nphiPL2     = nphiPL2L;
        nphiSaL2    = nphiSaL2L;
        iniu        = iniuL;
        iniP        = iniPL;
        iniSa       = iniSaL;
        iblockt      = 0;
    }
    else
    {
        this->ComputeProperties(datavecright, props);
        datavec = datavecright;
        phiuH1      = phiuH1R;
        phiPL2      = phiPL2R;
        phiSaL2     = phiSaL2R;
        dphiuH1     = dphiuH1R;
        dphiPL2     = dphiPL2R;
        nphiuHdiv   = nphiuHdivR;
        nphiPL2     = nphiPL2R;
        nphiSaL2    = nphiSaL2R;
        iniu        = iniuR;
        iniP        = iniPR;
        iniSa       = iniSaR;
        iblockt     = iblock;
    }
    
    TPZManVector<REAL> f_alpha          = props[3];
//    TPZFMatrix<STATE> jphiuHdiv(2,1);
//    int jshapeindex;
//    int jvectorindex;
    
    for (int isw = 0; isw < nphiSaL2L; isw++)
    {
        ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
        
//        for (int jq = 0; jq < nphiuHdiv; jq++)
//        {
//            jvectorindex = datavec[ublock].fVecShapeIndex[jq].first;
//            jshapeindex = datavec[ublock].fVecShapeIndex[jq].second;
//            jphiuHdiv(0,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(0,jvectorindex);
//            jphiuHdiv(1,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(1,jvectorindex);
//            REAL vn = jphiuHdiv(0,0)*n[0] + jphiuHdiv(1,0)*n[1];
//            
//            ek(isw + iniSaL,iblockt + jq + iniu) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * vn;
//            
//        }
        
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(isw + iniSaL,iblockt + jp + iniP) += 1.0 * weight * f_alpha[2] * phiPL2(jp,0)  * phiSaL2L(isw,0) * uLn;
        }
        
        for (int jsw = 0; jsw < nphiSaL2; jsw++)
        {
            ek(isw + iniSaL,iblockt + jsw + iniSa) += 1.0 * weight * f_alpha[3] * phiSaL2(jsw,0)  * phiSaL2L(isw,0) * uLn;
        }
    }
    
    for (int isw = 0; isw < nphiSaL2R; isw++)
    {
        ef(iblock + isw + iniSaR) += -1.0 * weight * f_alpha[0] * phiSaL2R(isw,0) * uLn;
        
//        for (int jq = 0; jq < nphiuHdiv; jq++)
//        {
//            jvectorindex = datavec[ublock].fVecShapeIndex[jq].first;
//            jshapeindex = datavec[ublock].fVecShapeIndex[jq].second;
//            jphiuHdiv(0,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(0,jvectorindex);
//            jphiuHdiv(1,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(1,jvectorindex);
//            REAL vn = jphiuHdiv(0,0)*n[0] + jphiuHdiv(1,0)*n[1];
//            
//            ek(iblock + isw + iniSaR,iblockt + jq + iniu) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * vn;
//            
//        }
        
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(iblock + isw + iniSaR, iblockt + jp + iniP ) += -1.0 * weight * f_alpha[2] * phiPL2(jp,0) * phiSaL2R(isw,0) * uLn;
        }
        
        for (int jsw = 0; jsw < nphiSaL2; jsw++)
        {
            ek(iblock + isw + iniSaR, iblockt + jsw + iniSa) += -1.0 * weight * f_alpha[3] * phiSaL2(jsw,0) * phiSaL2R(isw,0) * uLn;
        }
        
    }
    
    
    // Gravitational Segregation
    
    TPZVec<TPZManVector<REAL> > GravityFluxes;
    TPZManVector<REAL> fstar;
    this->GravitationalSegregation(data, datavecleft, datavecright, GravityFluxes,fstar);
    
    // Computing the minimum gravitational flux at the edge
    
    TPZManVector<REAL> qgLdotn = GravityFluxes[0];
    TPZManVector<REAL> qgRdotn = GravityFluxes[1];
    TPZManVector<REAL> gqdotn;

    
    // Taking the minimum module
    if (fstar[1] >= fstar[0]) {
        gqdotn = qgLdotn;
        phiuH1      = phiuH1L;
        phiPL2      = phiPL2L;
        phiSaL2     = phiSaL2L;
        dphiuH1     = dphiuH1L;
        dphiPL2     = dphiPL2L;
        nphiuHdiv   = nphiuHdivL;
        nphiPL2     = nphiPL2L;
        nphiSaL2    = nphiSaL2L;
        iniu        = iniuL;
        iniP        = iniPL;
        iniSa       = iniSaL;
        iblockt      = 0;
    }
    else
    {
        gqdotn = qgRdotn;
        phiuH1      = phiuH1R;
        phiPL2      = phiPL2R;
        phiSaL2     = phiSaL2R;
        dphiuH1     = dphiuH1R;
        dphiPL2     = dphiPL2R;
        nphiuHdiv   = nphiuHdivR;
        nphiPL2     = nphiPL2R;
        nphiSaL2    = nphiSaL2R;
        iniu        = iniuR;
        iniP        = iniPR;
        iniSa       = iniSaR;
        iblockt     = iblock;
    }
    
    // Comuting the contribution with the minimum gravitational flux

    for (int isw = 0; isw < nphiSaL2L; isw++)
    {
        ef(isw + iniSaL) += 1.0 * weight * gqdotn[0] * phiSaL2L(isw,0);
        
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(isw + iniSaL, jp + iniP + iblockt) += 1.0 * weight * gqdotn[2] * phiPL2(jp,0) * phiSaL2L(isw,0);
        }
        
        for (int jsw = 0; jsw < nphiSaL2; jsw++)
        {
            ek(isw + iniSaL, jsw + iniSa + iblockt) += 1.0 * weight * gqdotn[3] * phiSaL2(jsw,0) * phiSaL2L(isw,0);
        }
    }
    
    for (int isw = 0; isw < nphiSaL2R; isw++)
    {
        ef(iblock + isw + iniSaR) += -1.0 * weight * gqdotn[0] * phiSaL2R(isw,0);
        
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(iblock + isw + iniSaR, jp + iniP + iblockt) += -1.0 * weight * gqdotn[2] * phiPL2(jp,0) * phiSaL2L(isw,0);
        }
        
        for (int jsw = 0; jsw < nphiSaL2; jsw++)
        {
            ek(iblock + isw + iniSaR, jsw + iniSa + iblockt) += -1.0 * weight * gqdotn[3] * phiSaL2(jsw,0) * phiSaL2L(isw,0);
        }
    }
    
    
    // Capillary Pressure Pc_beta_beta
    TPZVec< TPZManVector<REAL>  > propsL;
    TPZVec< TPZManVector<REAL>  > propsR;
    
    ComputeProperties(datavecleft, propsL);
    TPZManVector<REAL> Pc_beta_alphaL          = propsL[9];
    
    ComputeProperties(datavecright, propsR);
    TPZManVector<REAL> Pc_beta_alphaR          = propsR[9];
    
    
    TPZManVector<REAL,3> vi(2,0.0);
    int ivectorindex = 0;
    int ishapeindex =0;
    
    for (int iq = 0; iq < nphiuHdivL; iq++)
    {
        ivectorindex = datavecleft[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavecleft[ublock].fVecShapeIndex[iq].second;
        vi[0] = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
        vi[1] = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
        REAL vin = vi[0]*n[0] + vi[1]*n[1];
        
        ef(iq + iniuL) += weight * (Salpha_R * Pc_beta_alphaR[0] - Salpha_L * Pc_beta_alphaL[0]) * vin;
        
        for (int jsw = 0; jsw < nphiSaL2L; jsw++)
        {
            ek(iq + iniuL, jsw + iniSaL) += -1.0 * weight * (Pc_beta_alphaL[0]+Salpha_L * Pc_beta_alphaL[3]) * phiSaL2L(jsw,0) * vin;
        }
        
        for (int jsw = 0; jsw < nphiSaL2R; jsw++)
        {
            ek(iq + iniuL, jsw + iniSaR + iblock) += 1.0 * weight * (Pc_beta_alphaR[0]+Salpha_R * Pc_beta_alphaR[3]) * phiSaL2L(jsw,0) * vin;
        }
        
    }
   
    // Aproximation for Grad of s
    TPZManVector<REAL,3> x_L = datavecleft[Sablock].XCenter;
    TPZManVector<REAL,3> x_R = datavecright[Sablock].XCenter;
    TPZManVector<REAL,3> x_dir(3,0.0);

    REAL delta_length = sqrt((x_L[0]-x_R[0])*(x_L[0]-x_R[0]) + (x_L[1]-x_R[1])*(x_L[1]-x_R[1]) + (x_L[2]-x_R[2])*(x_L[2]-x_R[2]));
    
//    if (delta_length == 0) {
//        std::cout << " delta_length = " << delta_length << std::endl;
//        DebugStop();
//    }
    
    
    REAL delta_s = Salpha_L - Salpha_R;
    REAL dsdl = delta_s / delta_length;
    
    x_dir[0] = (1.0/delta_length)*(x_L[0]-x_R[0]);
    x_dir[1] = (1.0/delta_length)*(x_L[1]-x_R[1]);
    x_dir[2] = (1.0/delta_length)*(x_L[2]-x_R[2]);
    
    // Capillary Segregation
    TPZVec< TPZManVector<REAL,3> > qc_L;
    TPZManVector<REAL,3> grads_L(3,0.0);
    //Computing a rough approximation for Grad of S
    grads_L[0] = dsdl * x_dir[0];
    grads_L[1] = dsdl * x_dir[1];
    grads_L[2] = dsdl * x_dir[2];
    this->CapillarySegregation(datavecleft,qc_L,grads_L);
    
    // Capillary Segregation
    TPZVec< TPZManVector<REAL,3> > qc_R;
    TPZManVector<REAL,3> grads_R(3,0.0);
    //Computing a rough approximation for Grad of S
    grads_R[0] = dsdl * x_dir[0];
    grads_R[1] = dsdl * x_dir[1];
    grads_R[2] = dsdl * x_dir[2];
    this->CapillarySegregation(datavecright,qc_R,grads_R);
    
    /* { flux } */
    REAL flux_c_L = qc_L[0][0]*n[0] + qc_L[0][1]*n[1];
    REAL flux_c_R = qc_R[0][0]*n[0] + qc_R[0][1]*n[1];
    REAL flux_avg = 0.5 * (flux_c_L + flux_c_R);
    
    // Comuting the contribution with the minimum gravitational flux
    
    for (int isw = 0; isw < nphiSaL2L; isw++)
    {
        ef(isw + iniSaL) += 1.0 * weight * (flux_avg) * phiSaL2L(isw,0);
    }
    for (int isw = 0; isw < nphiSaL2R; isw++)
    {
        ef(iblock + isw + iniSaR) += -1.0 * weight * (flux_avg) * phiSaL2R(isw,0);
    }
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::ContributeInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
    
    
//    // Full implicit case: there is no n state computations here
//    if (fSimulationData->IsnStep()) {
//        return;
//    }
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2L        = datavecleft[Sablock].phi; // For L2  test functions Sw
    TPZFMatrix<STATE> dphiuH1L   = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L   = datavecleft[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Getting test and basis functions for the right side
    TPZFMatrix<REAL> phiuH1R         = datavecright[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2R         = datavecright[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2R        = datavecright[Sablock].phi; // For L2  test functions Sa
    TPZFMatrix<STATE> dphiuH1R   = datavecright[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2R   = datavecright[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSaL2L    = phiSaL2L.Rows();                                   // For L2   Sa
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSaL   = nphiPL2L       + iniPL;
    
    // Blocks dimensions and lengths for the right side
    int nphiuHdivR   = datavecright[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2R     = phiPL2R.Rows();                                    // For L2   P
    int nphiSaL2R    = phiSaL2R.Rows();                                   // For L2   Sa
    int iniuR    = 0;
    int iniPR    = nphiuHdivR     + iniuR;
    int iniSaR   = nphiPL2R       + iniPR;
    
    int iblock = iniSaL + nphiSaL2L;
    
    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
//    REAL P_L             = datavecleft[Pblock].sol[0][0];
    REAL Salpha_L        = datavecleft[Sablock].sol[0][0];
    
    // Getting linear combinations from different approximation spaces for the right side
    TPZManVector<REAL,3> uR      = datavecright[ublock].sol[0];
//    REAL P_R             = datavecright[Pblock].sol[0][0];
    REAL Salpha_R        = datavecright[Sablock].sol[0][0];
    
//    // Computing the radius
//    TPZFMatrix<REAL> x_spatial(3,1,0.0);
//    x_spatial(0,0) = data.x[0];
//    REAL r = Norm(x_spatial);
//    REAL s = 1.0;
//    
//    if (fSimulationData->IsAxisymmetricQ()) {
//        s *= 2.0*M_PI*r;
//    }
//    
//    // Applying the scaling
//    phiuH1L      *= 1.0/s;
//    phiuH1R      *= 1.0/s;
//    
//    uL[0] *= 1.0/s;
//    uL[1] *= 1.0/s;
//    uL[2] *= 1.0/s;
//    
//    uR[0] *= 1.0/s;
//    uR[1] *= 1.0/s;
//    uR[2] *= 1.0/s;
    
    
    TPZManVector<REAL,3> n =  data.normal;
    REAL uLn = uL[0]*n[0] + uL[1]*n[1] + uL[2]*n[2];
    
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    
    if ((uLn >= 0.0) && fabs(uLn) > fepsilon)
    {
        this->ComputeProperties(datavecleft, props);
    }
    else
    {
        this->ComputeProperties(datavecright, props);
    }
    
    TPZManVector<REAL> f_alpha          = props[3];
    
    for (int isw = 0; isw < nphiSaL2L; isw++)
    {
        ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
    }
    
    for (int isw = 0; isw < nphiSaL2R; isw++)
    {
        ef(iblock + isw + iniSaR) += -1.0 * weight * f_alpha[0] * phiSaL2R(isw,0) * uLn;
    }
    
    
    
    
    // Gravitational Segregation
    
    TPZVec<TPZManVector<REAL> > GravityFluxes;
    TPZManVector<REAL> fstar;
    this->GravitationalSegregation(data, datavecleft, datavecright, GravityFluxes,fstar);
    
    // Computing the minimum gravitational flux at the edge
    
    TPZManVector<REAL> qgLdotn = GravityFluxes[0];
    TPZManVector<REAL> qgRdotn = GravityFluxes[1];
    TPZManVector<REAL> gqdotn;
    
    // Taking the minimum module
    if (fstar[1] >= fstar[0]) {
        gqdotn = qgLdotn;
    }
    else
    {
        gqdotn = qgRdotn;
    }
    
    // Comuting the contribution with the minimum gravitational flux
    
    for (int isw = 0; isw < nphiSaL2L; isw++)
    {
        ef(isw + iniSaL) += 1.0 * weight * gqdotn[0] * phiSaL2L(isw,0);
    }
    for (int isw = 0; isw < nphiSaL2R; isw++)
    {
        ef(iblock + isw + iniSaR) += -1.0 * weight * gqdotn[0] * phiSaL2R(isw,0);
    }
    
    
    // Capillary Pressure Pc_beta_alpha
    TPZVec< TPZManVector<REAL>  > propsL;
    TPZVec< TPZManVector<REAL>  > propsR;
    
    ComputeProperties(datavecleft, propsL);
    TPZManVector<REAL> Pc_beta_alphaL          = propsL[9];
    
    ComputeProperties(datavecright, propsR);
    TPZManVector<REAL> Pc_beta_alphaR          = propsR[9];

    
    TPZManVector<REAL,3> vi(2,0.0);
    int ivectorindex = 0;
    int ishapeindex =0;
    
    for (int iq = 0; iq < nphiuHdivL; iq++)
    {
        ivectorindex = datavecleft[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavecleft[ublock].fVecShapeIndex[iq].second;
        vi[0] = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
        vi[1] = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
        REAL vin = vi[0]*n[0] + vi[1]*n[1];
        
        ef(iq + iniuL) += weight * (Salpha_R * Pc_beta_alphaR[0] - Salpha_L * Pc_beta_alphaL[0]) * vin;
        
    }
    
    
    // Aproximation for Grad of s
    TPZManVector<REAL,3> x_L = datavecleft[Sablock].XCenter;
    TPZManVector<REAL,3> x_R = datavecright[Sablock].XCenter;
    TPZManVector<REAL,3> x_dir(3,0.0);
    
    REAL delta_length = sqrt((x_L[0]-x_R[0])*(x_L[0]-x_R[0]) + (x_L[1]-x_R[1])*(x_L[1]-x_R[1]) + (x_L[2]-x_R[2])*(x_L[2]-x_R[2]));
    REAL delta_s = Salpha_L - Salpha_R;
    REAL dsdl = delta_s / delta_length;
    
//    if (delta_length == 0) {
//        std::cout << " delta_length = " << delta_length << std::endl;
//        DebugStop();
//    }
    
    x_dir[0] = (1.0/delta_length)*(x_L[0]-x_R[0]);
    x_dir[1] = (1.0/delta_length)*(x_L[1]-x_R[1]);
    x_dir[2] = (1.0/delta_length)*(x_L[2]-x_R[2]);
    
    // Capillary Segregation
    TPZVec< TPZManVector<REAL,3> > qc_L;
    TPZManVector<REAL,3> grads_L(3,0.0);
    //Computing a rough approximation for Grad of S
    grads_L[0] = dsdl * x_dir[0];
    grads_L[1] = dsdl * x_dir[1];
    grads_L[2] = dsdl * x_dir[2];
    this->CapillarySegregation(datavecleft,qc_L,grads_L);
    
    // Capillary Segregation
    TPZVec< TPZManVector<REAL,3> > qc_R;
    TPZManVector<REAL,3> grads_R(3,0.0);
    //Computing a rough approximation for Grad of S
    grads_R[0] = dsdl * x_dir[0];
    grads_R[1] = dsdl * x_dir[1];
    grads_R[2] = dsdl * x_dir[2];
    this->CapillarySegregation(datavecright,qc_R,grads_R);
    
    /* { flux } */
    REAL flux_c_L = qc_L[0][0]*n[0] + qc_L[0][1]*n[1];
    REAL flux_c_R = qc_R[0][0]*n[0] + qc_R[0][1]*n[1];
    REAL flux_avg = 0.5 * (flux_c_L + flux_c_R);
    
    // Comuting the contribution with the minimum gravitational flux
    
    for (int isw = 0; isw < nphiSaL2L; isw++)
    {
        ef(isw + iniSaL) += 1.0 * weight * (flux_avg) * phiSaL2L(isw,0);
    }
    for (int isw = 0; isw < nphiSaL2R; isw++)
    {
        ef(iblock + isw + iniSaR) += -1.0 * weight * (flux_avg) * phiSaL2R(isw,0);
    }
    
}

void TPZAxiSymmetricDarcyFlow::GravitationalSegregation( TPZVec<TPZMaterialData> &datavec, TPZVec<TPZManVector<REAL,3> > & qg){
    
//    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    
    // Getting linear combinations from different approximation spaces for the left side
    REAL P             = datavec[Pblock].sol[0][0];
    REAL Salpha        = datavec[Sablock].sol[0][0];
    
    // Getting linear combinations from different approximation spaces for the right side
    
    // Gravitation Segregation
    int n_data = fnstate_vars + 2;
    TPZManVector<REAL> state_vars(n_data,0.0);
    TPZFMatrix<STATE> Gravity(2,1);
    TPZFMatrix<STATE> KGravity(2,1);
    qg.Resize(n_data,state_vars);
    
    Gravity = fSimulationData->GetGravity();
    
    // Computing Gravitational segregational function for the element

    TPZVec< TPZManVector<REAL>  > props;
    this->ComputeProperties(datavec, props);
    TPZFMatrix<STATE> K = fReservoirdata->Kabsolute();
    
    if (fSimulationData->IsHeterogeneousQ()) {
        int gelid = datavec[0].gelElId;
        K = fReservoirdata->GetKvector()[gelid];
    }
    
    TPZManVector<REAL> rho_alpha      = props[0];
    TPZManVector<REAL> rho_beta       = props[1];
    TPZManVector<REAL> lambda         = props[6];
    TPZManVector<REAL> rho_diff       = props[0];
    
    rho_diff[0] = (rho_alpha[0] - rho_beta[0]);
    rho_diff[1] = (rho_alpha[1] - rho_beta[1]);
    rho_diff[2] = (rho_alpha[2] - rho_beta[2]);
    rho_diff[3] = (rho_alpha[3] - rho_beta[3]);
    
    KGravity(0,0) = K(0,0)*Gravity(0,0) + K(0,1)*Gravity(1,0);
    KGravity(1,0) = K(1,0)*Gravity(0,0) + K(1,1)*Gravity(1,0);

    TPZManVector<REAL> f_value(n_data,0.0);
    TPZManVector<REAL> fstr(n_data,0.0);

    if (fSimulationData->IsLinearSegregationQ()) {
        fLinear(P, Salpha, f_value);
    }
    else{
        f(P, Salpha, f_value);
    }
    
    fstr[0] = ((f_value[0] * lambda[0]) * rho_diff[0]);
    fstr[1] = ((f_value[0] * lambda[0]) * rho_diff[1]) + ((f_value[0] * lambda[1] + f_value[1] * lambda[0]) * rho_diff[0]);
    fstr[2] = ((f_value[0] * lambda[0]) * rho_diff[2]) + ((f_value[0] * lambda[2] + f_value[2] * lambda[0]) * rho_diff[0]);
    fstr[3] = ((f_value[0] * lambda[0]) * rho_diff[3]) + ((f_value[0] * lambda[3] + f_value[3] * lambda[0]) * rho_diff[0]);
    

    qg[0][0] = fstr[0] * (KGravity(0,0));
    qg[0][1] = fstr[0] * (KGravity(1,0));
    
    qg[1][0] = fstr[1] * (KGravity(0,0));
    qg[1][1] = fstr[1] * (KGravity(1,0));
    
    qg[2][0] = fstr[2] * (KGravity(0,0));
    qg[2][1] = fstr[2] * (KGravity(1,0));
    
    qg[3][0] = fstr[3] * (KGravity(0,0));
    qg[3][1] = fstr[3] * (KGravity(1,0));
    
    return;
    
}


void TPZAxiSymmetricDarcyFlow::GravitationalSegregation(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft,TPZVec<TPZMaterialData> &datavecright, TPZVec<TPZManVector<REAL> > & GravityFluxes, TPZManVector<REAL> & fstar){
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    TPZManVector<REAL,3> n =  data.normal;
    
    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
    REAL P_L             = datavecleft[Pblock].sol[0][0];
    REAL Salpha_L        = datavecleft[Sablock].sol[0][0];
    
    // Getting linear combinations from different approximation spaces for the right side
    TPZManVector<REAL,3> uR      = datavecright[ublock].sol[0];
    REAL P_R             = datavecright[Pblock].sol[0][0];
    REAL Salpha_R        = datavecright[Sablock].sol[0][0];
    GravityFluxes.Resize(2);
    fstar.Resize(2);
    
    // Gravitation Segregation
    int n_data = fnstate_vars + 2;
    TPZManVector<REAL> state_vars(n_data,0.0);
    TPZFMatrix<STATE> Gravity(2,1);
    TPZFMatrix<STATE> KGravityL(2,1);
    TPZFMatrix<STATE> KGravityR(2,1);
    
    Gravity = fSimulationData->GetGravity();
    
    TPZFMatrix<STATE> KL = fReservoirdata->Kabsolute();
    
    if (fSimulationData->IsHeterogeneousQ()) {
        int gelid = datavecleft[0].gelElId;
        KL = fReservoirdata->GetKvector()[gelid];
    }
    
    REAL ndotG = Gravity(0,0)*n[0] + Gravity(1,0)*n[1];
    
    
    // Computing Gravitational segregational function on the left side
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props_Left;
    this->ComputeProperties(datavecleft, props_Left);
    
    TPZManVector<REAL> rho_alpha_L      = props_Left[0];
    TPZManVector<REAL> rho_beta_L       = props_Left[1];
    TPZManVector<REAL> rho_diff_L       = props_Left[0];
    
    rho_diff_L[0] = (rho_alpha_L[0] - rho_beta_L[0]);
    rho_diff_L[1] = (rho_alpha_L[1] - rho_beta_L[1]);
    rho_diff_L[2] = (rho_alpha_L[2] - rho_beta_L[2]);
    rho_diff_L[3] = (rho_alpha_L[3] - rho_beta_L[3]);
    
    KGravityL(0,0) = KL(0,0)*Gravity(0,0) + KL(0,1)*Gravity(1,0);
    KGravityL(1,0) = KL(1,0)*Gravity(0,0) + KL(1,1)*Gravity(1,0);
    
    TPZManVector<REAL> fstrL(n_data,0.0);
    TPZManVector<REAL> fL(n_data,0.0);
    
    if ((rho_alpha_L[0] - rho_beta_L[0]) <= 0.0 ) {
        ndotG *= -1.0;
    }
    
    
    if (ndotG <= 0.0) {
        
        // Receiving water ndotG <= 0.0
        if (fSimulationData->IsLinearSegregationQ()) {
            fRecLinear(P_L, Salpha_L, fL);
        }
        else{
            fRecLinear(P_L, Salpha_L, fL);
        }

    }
    else
    {
        // Expelling water ndotG >= 0.0
        if (fSimulationData->IsLinearSegregationQ()) {
            fExpLinear(P_L, Salpha_L, fL);
        }
        else{
            fExp(P_L, Salpha_L, fL);
        }
        
    }
    
    fstrL[0] = (fL[0] * rho_diff_L[0]);
    fstrL[1] = (fL[0] * rho_diff_L[1] + fL[1] * rho_diff_L[0]);
    fstrL[2] = (fL[0] * rho_diff_L[2] + fL[2] * rho_diff_L[0]);
    fstrL[3] = (fL[0] * rho_diff_L[3] + fL[3] * rho_diff_L[0]);
    
    
    // Computing Gravitational segregational function on the right side
    
    TPZFMatrix<STATE> KR = fReservoirdata->Kabsolute();
    
    if (fSimulationData->IsHeterogeneousQ()) {
        int gelid = datavecright[0].gelElId;
        KR = fReservoirdata->GetKvector()[gelid];
    }
    
    TPZVec< TPZManVector<REAL>  > props_Right;
    this->ComputeProperties(datavecright, props_Right);
    
    TPZManVector<REAL> rho_alpha_R        = props_Right[0];
    TPZManVector<REAL> rho_beta_R         = props_Right[1];
    TPZManVector<REAL> rho_diff_R    = props_Right[0];
    TPZManVector<REAL> lambda_R      = props_Right[6];
    
    rho_diff_R[0] = (rho_alpha_R[0] - rho_beta_R[0]);
    rho_diff_R[1] = (rho_alpha_R[1] - rho_beta_R[1]);
    rho_diff_R[2] = (rho_alpha_R[2] - rho_beta_R[2]);
    rho_diff_R[3] = (rho_alpha_R[3] - rho_beta_R[3]);
    
    KGravityR(0,0) = KR(0,0)*Gravity(0,0) + KR(0,1)*Gravity(1,0);
    KGravityR(1,0) = KR(1,0)*Gravity(0,0) + KR(1,1)*Gravity(1,0);
    
    TPZManVector<REAL> fstrR(n_data,0.0);
    TPZManVector<REAL> fR(n_data,0.0);
    
    
    if (ndotG <= 0.0) {
        
        // Expelling water ndotG >= 0.0 Linear Version
        
        if (fSimulationData->IsLinearSegregationQ()) {
            fExpLinear(P_R, Salpha_R, fR);
        }
        else{
            fExp(P_R, Salpha_R, fR);
        }
    }
    else
    {
        // Receive water ndotG <= 0.0 Linear Version
        
        if (fSimulationData->IsLinearSegregationQ()) {
            fRecLinear(P_R, Salpha_R, fR);
        }
        else{
            fRec(P_R, Salpha_R, fR);
        }
    }
    
    fstrR[0] = (fR[0] * rho_diff_R[0]);
    fstrR[1] = (fR[0] * rho_diff_R[1] + fR[1] * rho_diff_R[0]);
    fstrR[2] = (fR[0] * rho_diff_R[2] + fR[2] * rho_diff_R[0]);
    fstrR[3] = (fR[0] * rho_diff_R[3] + fR[3] * rho_diff_R[0]);
    
    
    TPZManVector<REAL> qgLdotn(n_data,0.0);
    TPZManVector<REAL> qgRdotn(n_data,0.0);
    
    fstar[0] = fstrL[0];
    fstar[1] = fstrR[0];
    
    if (fstar[1] >= fstar[0]) {
        
        qgLdotn[0] =  fstrL[0] * (KGravityL(0,0)*n[0] + KGravityL(1,0)*n[1]);
        qgLdotn[1] =  fstrL[1] * (KGravityL(0,0)*n[0] + KGravityL(1,0)*n[1]);
        qgLdotn[2] =  fstrL[2] * (KGravityL(0,0)*n[0] + KGravityL(1,0)*n[1]);
        qgLdotn[3] =  fstrL[3] * (KGravityL(0,0)*n[0] + KGravityL(1,0)*n[1]);
        
        qgRdotn[0] =  fstrL[0] * (KGravityR(0,0)*n[0] + KGravityR(1,0)*n[1]);
        qgRdotn[1] =  fstrL[1] * (KGravityR(0,0)*n[0] + KGravityR(1,0)*n[1]);
        qgRdotn[2] =  fstrL[2] * (KGravityR(0,0)*n[0] + KGravityR(1,0)*n[1]);
        qgRdotn[3] =  fstrL[3] * (KGravityR(0,0)*n[0] + KGravityR(1,0)*n[1]);
        
    }
    else
    {
        
        qgLdotn[0] =  fstrR[0] * (KGravityL(0,0)*n[0] + KGravityL(1,0)*n[1]);
        qgLdotn[1] =  fstrR[1] * (KGravityL(0,0)*n[0] + KGravityL(1,0)*n[1]);
        qgLdotn[2] =  fstrR[2] * (KGravityL(0,0)*n[0] + KGravityL(1,0)*n[1]);
        qgLdotn[3] =  fstrR[3] * (KGravityL(0,0)*n[0] + KGravityL(1,0)*n[1]);
        
        qgRdotn[0] =  fstrR[0] * (KGravityR(0,0)*n[0] + KGravityR(1,0)*n[1]);
        qgRdotn[1] =  fstrR[1] * (KGravityR(0,0)*n[0] + KGravityR(1,0)*n[1]);
        qgRdotn[2] =  fstrR[2] * (KGravityR(0,0)*n[0] + KGravityR(1,0)*n[1]);
        qgRdotn[3] =  fstrR[3] * (KGravityR(0,0)*n[0] + KGravityR(1,0)*n[1]);
        
    }
    
    GravityFluxes[0] = qgLdotn;
    GravityFluxes[1] = qgRdotn;
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::CapillarySegregation(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZManVector<REAL,3> > & qc, TPZManVector<REAL,3> & grads){
    
    //    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    
    // Getting linear combinations from different approximation spaces for the left side
    REAL P             = datavec[Pblock].sol[0][0];
    REAL Salpha        = datavec[Sablock].sol[0][0];
    
    TPZFMatrix<STATE> GradSaxes = datavec[Sablock].dsol[0];
    
    // Getting linear combinations from different approximation spaces for the right side
    
    // Gravitation Segregation
    int n_data = fnstate_vars + 2;
    TPZManVector<REAL> state_vars(n_data,0.0);
    qc.Resize(n_data,state_vars);
    // Computing Gravitational segregational function for the element
    
    TPZVec< TPZManVector<REAL>  > props;
    this->ComputeProperties(datavec, props);
    TPZFMatrix<STATE> K = fReservoirdata->Kabsolute();
    
    if (fSimulationData->IsHeterogeneousQ()) {
        int gelid = datavec[0].gelElId;
        K = fReservoirdata->GetKvector()[gelid];
    }
    
    TPZManVector<REAL> rho_alpha      = props[0];
    TPZManVector<REAL> rho_beta       = props[1];
    TPZManVector<REAL> lambda         = props[6];
    TPZManVector<REAL> Pc_beta_alpha  = props[9];
    

    TPZManVector<STATE> GradS(2,0.0);
    TPZManVector<STATE> Grad_Pc(2,0.0);
    TPZFMatrix<STATE>   K_Grad_Pc(2,1);

    //  Compute grad(S)
    if (fSimulationData->GetGR()) {
        GradS[0] = GradSaxes(0,0)*datavec[Sablock].axes(0,0)+GradSaxes(1,0)*datavec[Sablock].axes(1,0);
        GradS[1] = GradSaxes(0,0)*datavec[Sablock].axes(0,1)+GradSaxes(1,0)*datavec[Sablock].axes(1,1);
    }
    else{
//        GradS[0] = fSimulationData->Length_Scale()*grads[0];
//        GradS[1] = fSimulationData->Length_Scale()*grads[1];
        
        GradS[0] = grads[0];
        GradS[1] = grads[1];
    }
    
    
    //  Compute grad(Pc)
    Grad_Pc[0] = Pc_beta_alpha[3] * GradS[0];
    Grad_Pc[1] = Pc_beta_alpha[3] * GradS[1];
    
    K_Grad_Pc(0,0) = K(0,0)*Grad_Pc[0] + K(0,1)*Grad_Pc[1];
    K_Grad_Pc(1,0) = K(1,0)*Grad_Pc[0] + K(1,1)*Grad_Pc[1];
    
    TPZManVector<REAL> f_value(n_data,0.0);
    TPZManVector<REAL> fstr(n_data,0.0);
    
    if (fSimulationData->IsLinearSegregationQ()) {
        fLinear(P, Salpha, f_value);
    }
    else{
        f(P, Salpha, f_value);
    }
    
    fstr[0] = (f_value[0] * lambda[0]);
    fstr[1] = (f_value[0] * lambda[1] + f_value[1] * lambda[0]);
    fstr[2] = (f_value[0] * lambda[2] + f_value[2] * lambda[0]);
    fstr[3] = (f_value[0] * lambda[3] + f_value[3] * lambda[0]);
    

    qc[0][0] = fstr[0] * (K_Grad_Pc(0,0));
    qc[0][1] = fstr[0] * (K_Grad_Pc(1,0));
    
    qc[1][0] = fstr[1] * (K_Grad_Pc(0,0));
    qc[1][1] = fstr[1] * (K_Grad_Pc(1,0));
    
    qc[2][0] = fstr[2] * (K_Grad_Pc(0,0));
    qc[2][1] = fstr[2] * (K_Grad_Pc(1,0));
    
    qc[3][0] = fstr[3] * (K_Grad_Pc(0,0));
    qc[3][1] = fstr[3] * (K_Grad_Pc(1,0));
    
    qc[4][0] = fstr[0] * Pc_beta_alpha[0];
    qc[4][1] = fstr[0] * Pc_beta_alpha[0];
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::CapillarySegregation(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft,TPZVec<TPZMaterialData> &datavecright, TPZVec<TPZManVector<REAL> > & CapillaryFluxes, TPZManVector<REAL> & fstar){
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    TPZManVector<REAL,3> n =  data.normal;
    
    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
    REAL P_L             = datavecleft[Pblock].sol[0][0];
    REAL Salpha_L        = datavecleft[Sablock].sol[0][0];
    
    TPZFMatrix<STATE> GradSaxes_L = datavecleft[Sablock].dsol[0];
    
    // Getting linear combinations from different approximation spaces for the right side
    TPZManVector<REAL,3> uR      = datavecright[ublock].sol[0];
    REAL P_R             = datavecright[Pblock].sol[0][0];
    REAL Salpha_R        = datavecright[Sablock].sol[0][0];
    
    TPZFMatrix<STATE> GradSaxes_R = datavecright[Sablock].dsol[0];
    TPZFMatrix<STATE> KL = fReservoirdata->Kabsolute();
    TPZFMatrix<STATE> KR = fReservoirdata->Kabsolute();
    
    
    
    CapillaryFluxes.Resize(2);
    fstar.Resize(2);
    
    // Gravitation Segregation
    int n_data = fnstate_vars + 2;
    TPZManVector<REAL> state_vars(n_data,0.0);

    if (fSimulationData->IsHeterogeneousQ()) {
        int gelid = datavecleft[0].gelElId;
        KL = fReservoirdata->GetKvector()[gelid];
    }

    // Computing Gravitational segregational function on the left side
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props_Left;
    this->ComputeProperties(datavecleft, props_Left);
    TPZManVector<REAL> lambda_L             = props_Left[6];
    TPZManVector<REAL> Pc_beta_alpha_L      = props_Left[9];
    
    TPZManVector<STATE> GradS_L(2,0.0);
    TPZManVector<STATE> Grad_Pc_L(2,0.0);
    TPZFMatrix<STATE>   K_Grad_Pc_L(2,1);
    
    //  Compute grad(S_L)
    GradS_L[0] = GradSaxes_L(0,0)*datavecleft[Sablock].axes(0,0)+GradSaxes_L(1,0)*datavecleft[Sablock].axes(1,0);
    GradS_L[1] = GradSaxes_L(0,0)*datavecleft[Sablock].axes(0,1)+GradSaxes_L(1,0)*datavecleft[Sablock].axes(1,1);
    
    //  Compute grad(Pc_L)
    Grad_Pc_L[0] = Pc_beta_alpha_L[3] * GradS_L[0];
    Grad_Pc_L[1] = Pc_beta_alpha_L[3] * GradS_L[1];
    
    TPZManVector<REAL> fstrL(n_data,0.0);
    TPZManVector<REAL> fL(n_data,0.0);
    
    TPZManVector<REAL,3> KGradPc_L(2,0.0);
    KGradPc_L[0] = KL(0,0)*Grad_Pc_L[0] + KL(0,1)*Grad_Pc_L[1];
    KGradPc_L[1] = KL(1,0)*Grad_Pc_L[0] + KL(1,1)*Grad_Pc_L[1];
    
    REAL GradPc_Ldotn = Grad_Pc_L[0]*n[0] + Grad_Pc_L[1]*n[1];
    
    if (GradPc_Ldotn <= 0.0) {
        
        // Receiving water ndotG <= 0.0
        if (fSimulationData->IsLinearSegregationQ()) {
            fRecLinear(P_L, Salpha_L, fL);
        }
        else{
            fRecLinear(P_L, Salpha_L, fL);
        }
        
    }
    else
    {
        // Expelling water ndotG >= 0.0
        if (fSimulationData->IsLinearSegregationQ()) {
            fExpLinear(P_L, Salpha_L, fL);
        }
        else{
            fExp(P_L, Salpha_L, fL);
        }
        
    }
    
    fstrL[0] = (fL[0] * lambda_L[0]);
    fstrL[1] = (fL[0] * lambda_L[1] + fL[1] * lambda_L[0]);
    fstrL[2] = (fL[0] * lambda_L[2] + fL[2] * lambda_L[0]);
    fstrL[3] = (fL[0] * lambda_L[3] + fL[3] * lambda_L[0]);
    
    
    // Computing Gravitational segregational function on the right side
    
    if (fSimulationData->IsHeterogeneousQ()) {
        int gelid = datavecright[0].gelElId;
        KR = fReservoirdata->GetKvector()[gelid];
    }
    
    TPZVec< TPZManVector<REAL>  > props_Right;
    this->ComputeProperties(datavecright, props_Right);
    
    TPZManVector<REAL> lambda_R             = props_Right[6];
    TPZManVector<REAL> Pc_beta_alpha_R      = props_Right[9];
    
    TPZManVector<STATE> GradS_R(2,0.0);
    TPZManVector<STATE> Grad_Pc_R(2,0.0);
    TPZFMatrix<STATE>   K_Grad_Pc_R(2,1);
    
    //  Compute grad(S_R)
    GradS_R[0] = GradSaxes_R(0,0)*datavecright[Sablock].axes(0,0)+GradSaxes_R(1,0)*datavecright[Sablock].axes(1,0);
    GradS_R[1] = GradSaxes_R(0,0)*datavecright[Sablock].axes(0,1)+GradSaxes_R(1,0)*datavecright[Sablock].axes(1,1);
    
    //  Compute grad(Pc_L)
    Grad_Pc_R[0] = Pc_beta_alpha_R[3] * GradS_R[0];
    Grad_Pc_R[1] = Pc_beta_alpha_R[3] * GradS_R[1];
    
    TPZManVector<REAL> fstrR(n_data,0.0);
    TPZManVector<REAL> fR(n_data,0.0);
    
    TPZManVector<REAL,3> KGradPc_R(2,0.0);
    KGradPc_R[0] = KR(0,0)*Grad_Pc_R[0] + KR(0,1)*Grad_Pc_R[1];
    KGradPc_R[1] = KR(1,0)*Grad_Pc_R[0] + KR(1,1)*Grad_Pc_R[1];
    
    REAL GradPc_Rdotn = Grad_Pc_R[0]*n[0] + Grad_Pc_R[1]*n[1];
    
    if (GradPc_Rdotn <= 0.0) {
        
        // Expelling water ndotG >= 0.0 Linear Version
        
        if (fSimulationData->IsLinearSegregationQ()) {
            fExpLinear(P_R, Salpha_R, fR);
        }
        else{
            fExp(P_R, Salpha_R, fR);
        }
    }
    else
    {
        // Receive water ndotG <= 0.0 Linear Version
        
        if (fSimulationData->IsLinearSegregationQ()) {
            fRecLinear(P_R, Salpha_R, fR);
        }
        else{
            fRec(P_R, Salpha_R, fR);
        }
    }
    
    fstrR[0] = (fR[0] * lambda_R[0]);
    fstrR[1] = (fR[0] * lambda_R[1] + fR[1] * lambda_R[0]);
    fstrR[2] = (fR[0] * lambda_R[2] + fR[2] * lambda_R[0]);
    fstrR[3] = (fR[0] * lambda_R[3] + fR[3] * lambda_R[0]);
    
    
    TPZManVector<REAL> qcLdotn(n_data,0.0);
    TPZManVector<REAL> qcRdotn(n_data,0.0);
    
    fstar[0] = fstrL[0];
    fstar[1] = fstrR[0];
    
    if (fstar[1] >= fstar[0]) {
        
        qcLdotn[0] =  fstrL[0] * (KGradPc_L[0]*n[0] + KGradPc_L[1]*n[1]);
        qcLdotn[1] =  fstrL[1] * (KGradPc_L[0]*n[0] + KGradPc_L[1]*n[1]);
        qcLdotn[2] =  fstrL[2] * (KGradPc_L[0]*n[0] + KGradPc_L[1]*n[1]);
        qcLdotn[3] =  fstrL[3] * (KGradPc_L[0]*n[0] + KGradPc_L[1]*n[1]);
        
        qcRdotn[0] =  fstrL[0] * (KGradPc_L[0]*n[0] + KGradPc_L[1]*n[1]);
        qcRdotn[1] =  fstrL[1] * (KGradPc_L[0]*n[0] + KGradPc_L[1]*n[1]);
        qcRdotn[2] =  fstrL[2] * (KGradPc_L[0]*n[0] + KGradPc_L[1]*n[1]);
        qcRdotn[3] =  fstrL[3] * (KGradPc_L[0]*n[0] + KGradPc_L[1]*n[1]);
        
    }
    else
    {
        
        qcLdotn[0] =  fstrR[0] * (KGradPc_R[0]*n[0] + KGradPc_R[1]*n[1]);
        qcLdotn[1] =  fstrR[1] * (KGradPc_R[0]*n[0] + KGradPc_R[1]*n[1]);
        qcLdotn[2] =  fstrR[2] * (KGradPc_R[0]*n[0] + KGradPc_R[1]*n[1]);
        qcLdotn[3] =  fstrR[3] * (KGradPc_R[0]*n[0] + KGradPc_R[1]*n[1]);
        
        qcRdotn[0] =  fstrR[0] * (KGradPc_R[0]*n[0] + KGradPc_R[1]*n[1]);
        qcRdotn[1] =  fstrR[1] * (KGradPc_R[0]*n[0] + KGradPc_R[1]*n[1]);
        qcRdotn[2] =  fstrR[2] * (KGradPc_R[0]*n[0] + KGradPc_R[1]*n[1]);
        qcRdotn[3] =  fstrR[3] * (KGradPc_R[0]*n[0] + KGradPc_R[1]*n[1]);
        
    }
    
    CapillaryFluxes[0] = qcLdotn;
    CapillaryFluxes[1] = qcRdotn;
    
//    std::cout << "Grad_Pc_L = " << Grad_Pc_L << std::endl;
//    std::cout << "Grad_Pc_R = " << Grad_Pc_R << std::endl;
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::fRec(REAL P, REAL Salpha, TPZManVector<REAL> & ExpL){
    
    int n_data = fnstate_vars + 2;
    ExpL.Resize(n_data, 0.0);
    
    TPZManVector<REAL> state_vars(n_data,0.0);
    TPZVec< TPZManVector<REAL>  > props;
    TPZManVector<REAL> f_alpha;
    TPZManVector<REAL> f_beta;
    TPZManVector<REAL> lambda;
    
    // Receive water ndotG <= 0.0 Linear Version
    if (Salpha < fSalpha_max) {
        
        state_vars[0] = 0.0;            //  qn
        state_vars[1] = P;            //  P
        state_vars[2] = fSalpha_max;        //  S_alpha
        state_vars[3] = fSalpha_max;        //  S_beta
        
        this->ComputeProperties(state_vars, props);
        f_alpha         = props[3];
        f_beta          = props[4];
        lambda          = props[6];
        f_alpha[3]       = 0.0;
        f_beta[3]        = 0.0;
        lambda[3]        = 0.0;

        ExpL[0] =  (f_alpha[0] * f_beta[0]) * lambda[0];
        ExpL[1] =  (f_alpha[0] * f_beta[0]) * lambda[1] + (f_alpha[0] * f_beta[1] + f_alpha[1] * f_beta[0]) * lambda[0];
        ExpL[2] =  (f_alpha[0] * f_beta[0]) * lambda[2] + (f_alpha[0] * f_beta[2] + f_alpha[2] * f_beta[0]) * lambda[0];
        ExpL[3] =  (f_alpha[0] * f_beta[0]) * lambda[3] + (f_alpha[0] * f_beta[3] + f_alpha[3] * f_beta[0]) * lambda[0];
        
    }
    else
    {
        
        state_vars[0] = 0.0;            //  qn
        state_vars[1] = P;            //  P
        state_vars[2] = Salpha;        //  S_alpha
        state_vars[3] = Salpha;        //  S_beta
        
        this->ComputeProperties(state_vars, props);
        f_alpha         = props[3];
        f_beta          = props[4];
        lambda          = props[6];

        
        ExpL[0] =  (f_alpha[0] * f_beta[0]) * lambda[0];
        ExpL[1] =  (f_alpha[0] * f_beta[0]) * lambda[1] + (f_alpha[0] * f_beta[1] + f_alpha[1] * f_beta[0]) * lambda[0];
        ExpL[2] =  (f_alpha[0] * f_beta[0]) * lambda[2] + (f_alpha[0] * f_beta[2] + f_alpha[2] * f_beta[0]) * lambda[0];
        ExpL[3] =  (f_alpha[0] * f_beta[0]) * lambda[3] + (f_alpha[0] * f_beta[3] + f_alpha[3] * f_beta[0]) * lambda[0];
        
    }
    
}

void TPZAxiSymmetricDarcyFlow::fExp(REAL P, REAL Salpha, TPZManVector<REAL> & RecL){
    
    int n_data = fnstate_vars + 2;
    RecL.Resize(n_data, 0.0);
    
    TPZManVector<REAL> state_vars(n_data,0.0);
    TPZVec< TPZManVector<REAL>  > props;
    TPZManVector<REAL> f_alpha;
    TPZManVector<REAL> f_beta;
    TPZManVector<REAL> lambda;
    
    
    if (Salpha > fSalpha_max) {
        
        state_vars[0] = 0.0;            //  qn
        state_vars[1] = P;            //  P
        state_vars[2] = fSalpha_max;        //  S_alpha
        state_vars[3] = fSalpha_max;        //  S_beta
        
        this->ComputeProperties(state_vars, props);
        f_alpha         = props[3];
        f_beta          = props[4];
        lambda          = props[6];
        f_alpha[3]       = 0.0;
        f_beta[3]        = 0.0;
        lambda[3]        = 0.0;
        
        RecL[0] =  (f_alpha[0] * f_beta[0]) * lambda[0];
        RecL[1] =  (f_alpha[0] * f_beta[0]) * lambda[1] + (f_alpha[0] * f_beta[1] + f_alpha[1] * f_beta[0]) * lambda[0];
        RecL[2] =  (f_alpha[0] * f_beta[0]) * lambda[2] + (f_alpha[0] * f_beta[2] + f_alpha[2] * f_beta[0]) * lambda[0];
        RecL[3] =  (f_alpha[0] * f_beta[0]) * lambda[3] + (f_alpha[0] * f_beta[3] + f_alpha[3] * f_beta[0]) * lambda[0];
        
    }
    else
    {
        
        state_vars[0] = 0.0;            //  qn
        state_vars[1] = P;            //  P
        state_vars[2] = Salpha;        //  S_alpha
        state_vars[3] = Salpha;        //  S_beta
        
        this->ComputeProperties(state_vars, props);
        f_alpha         = props[3];
        f_beta          = props[4];
        lambda          = props[6];
        
        RecL[0] =  (f_alpha[0] * f_beta[0]) * lambda[0];
        RecL[1] =  (f_alpha[0] * f_beta[0]) * lambda[1] + (f_alpha[0] * f_beta[1] + f_alpha[1] * f_beta[0]) * lambda[0];
        RecL[2] =  (f_alpha[0] * f_beta[0]) * lambda[2] + (f_alpha[0] * f_beta[2] + f_alpha[2] * f_beta[0]) * lambda[0];
        RecL[3] =  (f_alpha[0] * f_beta[0]) * lambda[3] + (f_alpha[0] * f_beta[3] + f_alpha[3] * f_beta[0]) * lambda[0];
        
    }
}

void TPZAxiSymmetricDarcyFlow::fRecLinear(REAL P, REAL Salpha, TPZManVector<REAL> & ExpL){
    
    int n_data = fnstate_vars + 2;
    ExpL.Resize(n_data, 0.0);
    REAL fs_Smax = 0.0;
    
    TPZManVector<REAL> state_vars(n_data,0.0);
    TPZVec< TPZManVector<REAL>  > props;
    TPZManVector<REAL> f_alpha;
    TPZManVector<REAL> f_beta;
    TPZManVector<REAL> lambda;
    
    // Receive water ndotG <= 0.0 Linear Version
    if (Salpha < fSalpha_max) {
        
        state_vars[0] = 0.0;            //  qn
        state_vars[1] = P;            //  P
        state_vars[2] = fSalpha_max;        //  S_alpha
        state_vars[3] = fSalpha_max;        //  S_beta
        
        this->ComputeProperties(state_vars, props);
        f_alpha         = props[3];
        f_beta          = props[4];
        lambda          = props[6];
        f_alpha[3]       = 0.0;
        f_beta[3]        = 0.0;
        lambda[3]        = 0.0;
        
        fs_Smax = f_alpha[0] * f_beta[0] * lambda[0];
        
        ExpL[0] = fs_Smax;
        ExpL[1] = 0.0;
        ExpL[2] = ((f_alpha[0] * f_beta[0]) * lambda[2] + (f_alpha[0] * f_beta[2] + f_alpha[2] * f_beta[0] ) * lambda[0]);
        ExpL[3] = 0.0;
        
    }
    else
    {
        
        state_vars[0] = 0.0;            //  qn
        state_vars[1] = P;            //  P
        state_vars[2] = fSalpha_max;        //  S_alpha
        state_vars[3] = fSalpha_max;        //  S_beta
        
        this->ComputeProperties(state_vars, props);
        f_alpha         = props[3];
        f_beta          = props[4];
        lambda          = props[6];
        f_alpha[3]       = 0.0;
        f_beta[3]        = 0.0;
        lambda[3]        = 0.0;
        
        fs_Smax = f_alpha[0] * f_beta[0] * lambda[0];
        REAL denom = 1.0-fSalpha_max-fReservoirdata->S_nwett_r();
        ExpL[0] = (1.0-Salpha-fReservoirdata->S_nwett_r())*(fs_Smax/(denom));
        ExpL[1] = 0.0;
        ExpL[2] = (1.0-Salpha-fReservoirdata->S_nwett_r())*(1.0/(denom))*((f_alpha[0] * f_beta[0]) * lambda[2] + (f_alpha[0] * f_beta[2] + f_alpha[2] * f_beta[0] ) * lambda[0]);
        ExpL[3] = (- 1.0)*(fs_Smax/(denom));
        
    }
    
}

void TPZAxiSymmetricDarcyFlow::fExpLinear(REAL P, REAL Salpha, TPZManVector<REAL> & RecL){
    
    int n_data = fnstate_vars + 2;
    RecL.Resize(n_data, 0.0);
    REAL fs_Smax = 0.0;
    
    TPZManVector<REAL> state_vars(n_data,0.0);
    TPZVec< TPZManVector<REAL>  > props;
    TPZManVector<REAL> f_alpha;
    TPZManVector<REAL> f_beta;
    TPZManVector<REAL> lambda;
    

    if (Salpha > fSalpha_max) {
        
        state_vars[0] = 0.0;            //  qn
        state_vars[1] = P;            //  P
        state_vars[2] = fSalpha_max;        //  S_alpha
        state_vars[3] = fSalpha_max;        //  S_beta
        
        this->ComputeProperties(state_vars, props);
        f_alpha         = props[3];
        f_beta          = props[4];
        lambda          = props[6];
        f_alpha[3]       = 0.0;
        f_beta[3]        = 0.0;
        lambda[3]        = 0.0;
        
        fs_Smax = f_alpha[0] * f_beta[0] * lambda[0];
        
        RecL[0] = fs_Smax;
        RecL[1] = 0.0;
        RecL[2] = ((f_alpha[0] * f_beta[0]) * lambda[2] + (f_alpha[0] * f_beta[2] + f_alpha[2] * f_beta[0] ) * lambda[0]);
        RecL[3] = 0.0;
        
    }
    else
    {
        
        state_vars[0] = 0.0;            //  qn
        state_vars[1] = P;            //  P
        state_vars[2] = fSalpha_max;        //  S_alpha
        state_vars[3] = fSalpha_max;        //  S_beta
        
        this->ComputeProperties(state_vars, props);
        f_alpha         = props[3];
        f_beta          = props[4];
        lambda          = props[6];
        f_alpha[3]       = 0.0;
        f_beta[3]        = 0.0;
        lambda[3]        = 0.0;
        
        fs_Smax = f_alpha[0] * f_beta[0] * lambda[0];
        REAL denom = fSalpha_max - fReservoirdata->S_wett_r();
        RecL[0] = (Salpha-fReservoirdata->S_wett_r())*(fs_Smax/(denom));
        RecL[1] = 0.0;
        RecL[2] = (Salpha-fReservoirdata->S_wett_r())*(1.0/(denom))*((f_alpha[0] * f_beta[0]) * lambda[2] + (f_alpha[0] * f_beta[2] + f_alpha[2] * f_beta[0] ) * lambda[0]);
        RecL[3] = (1.0)*(fs_Smax/(denom));
        
    }
    
}

void TPZAxiSymmetricDarcyFlow::f(REAL P, REAL Salpha, TPZManVector<REAL> & f_value){
    
    int n_data = fnstate_vars + 2;
    
    TPZManVector<REAL> state_vars(n_data,0.0);
    TPZVec< TPZManVector<REAL>  > props;
    TPZManVector<REAL> f_alpha;
    TPZManVector<REAL> f_beta;
    
    f_value.Resize(n_data, 0.0);
    
    state_vars[0] = 0.0;            //  qn
    state_vars[1] = P;            //  P
    state_vars[2] = Salpha;        //  S_alpha
    state_vars[3] = 1.0 - Salpha;        //  S_beta
    
    this->ComputeProperties(state_vars, props);
    f_alpha         = props[3];
    f_beta          = props[4];
    
    f_value[0] = f_alpha[0] * f_beta[0];
    f_value[1] = (f_alpha[0] * f_beta[1] + f_alpha[1] * f_beta[0]);
    f_value[2] = (f_alpha[0] * f_beta[2] + f_alpha[2] * f_beta[0]);
    f_value[3] = (f_alpha[0] * f_beta[3] + f_alpha[3] * f_beta[0]);
    
    
    
}

void TPZAxiSymmetricDarcyFlow::fLinear(REAL P, REAL Salpha, TPZManVector<REAL> & f_value){
    
    int n_data = fnstate_vars + 2;
    REAL fs_Smax = 0.0;
    
    TPZManVector<REAL> state_vars(n_data,0.0);
    TPZVec< TPZManVector<REAL>  > props;
    TPZManVector<REAL> f_alpha;
    TPZManVector<REAL> f_beta;
    
    f_value.Resize(n_data, 0.0);
    
    state_vars[0] = 0.0;            //  qn
    state_vars[1] = P;            //  P
    state_vars[2] = fSalpha_max;        //  S_alpha
    state_vars[3] = 1.0 - fSalpha_max;        //  S_beta
    
    this->ComputeProperties(state_vars, props);
    f_alpha         = props[3];
    f_beta          = props[4];
    
    fs_Smax = f_alpha[0] * f_beta[0];
    
    if (Salpha < fSalpha_max) {
        f_value[0] = (fSalpha_max/fSalpha_max)*Salpha;
        f_value[1] = 0.0;
        f_value[2] = 0.0;
        f_value[3] = (fSalpha_max/fSalpha_max);
    }
    else{
        f_value[0] = (1.0 - Salpha)*(fs_Smax/(1.0-fSalpha_max));
        f_value[1] = 0.0;
        f_value[2] = 0.0;
        f_value[3] = (-1.0)*(fs_Smax/(1.0-fSalpha_max));
    }
    

    
}

void TPZAxiSymmetricDarcyFlow::UpdateStateVariables(TPZManVector<REAL,3> u, REAL P, REAL Sw, REAL So)
{
    fBulkVelocity       = u;
    fAveragePressure    = P;
    fWaterSaturation    = Sw;
    fOilSaturation      = So;
}

// systematic methods computed on the current u, p, Sw and So values.



void TPZAxiSymmetricDarcyFlow::PhasePressures()
{
    // Petrophysics data
    
    REAL Sw, So, Sg, P;
    P = fAveragePressure;
    Sw = fWaterSaturation;
    So = fOilSaturation;
    Sg = 1.0 - Sw - So;
    
    REAL pcow = 0, pcgo = 0, pcgw = 0;
    REAL dPcowdSw = 0, dPcgodSo = 0;
    
    //    this->fPetrophysicdata->Pcow(Sw, pcow, dPcowdSw);
    //    this->fPetrophysicdata->Pcgo(So, pcgo, dPcgodSo);
    //    this->fPetrophysicdata->Pcgw(Sw, pcgw, dPcgwdSw); // or Pcgo(So) + Pcow(Sw)
    
    // Recovering phase pressures and their derivatives
    
    fWaterPressure[0] = P - So * pcow - (1.0 - So - Sw) * pcgw;
    fWaterPressure[1] = 1.0;
    fWaterPressure[2] = pcgo - So * dPcowdSw + pcow - Sg * dPcowdSw;
    fWaterPressure[3] = pcgo - Sg * dPcgodSo;
    
    fOilPressure[0] = P + Sw * pcow - (1.0 - So - Sw) * pcgo;
    fOilPressure[1] = 1.0;
    fOilPressure[2] = Sw * dPcowdSw + pcow + pcgo;
    fOilPressure[3] = pcgo - Sg * dPcgodSo;
    
    // Computing fluid mobilities and derivatives for two-phase
    
    fGasPressure[0] = 0.0;
    fGasPressure[1] = 0.0;
    fGasPressure[2] = 0.0;
    fGasPressure[3] = 0.0;
    
    
    //    // Computing fluid mobilities and derivatives for three-phase
    //
    //    fGasPressure[0] = P + So * pcgo + So * pcgw;
    //    fGasPressure[1] = 1.0;
    //    fGasPressure[2] = pcgo + Sw * dPcowdSw + pcow;
    //    fGasPressure[3] = pcgo + So * dPcgodSo + Sw * dPcgodSo;
    
}

void TPZAxiSymmetricDarcyFlow::PhaseDensities()
{
    
    this->PhasePressures();
    
    //    fFluidmodeldata->WaterDensity(fWaterPressure[0], fWaterDensity[0], fWaterDensity[1]);
    //    fFluidmodeldata->OilDensity(fOilPressure[0], fOilDensity[0], fOilDensity[1]);
    //    fFluidmodeldata->GasDensity(fGasPressure[0], fGasDensity[0], fGasDensity[1]);
    
    // Chain Rule
    fWaterDensity[1]    *= fWaterPressure[1];
    fOilDensity[1]      *= fOilPressure[1];
    fGasDensity[1]      *= fGasPressure[1];
}

void TPZAxiSymmetricDarcyFlow::PhaseMobilities()
{
    REAL Sw, So, Sg;
    Sw = fWaterSaturation;
    So = fOilSaturation;
    Sg = 1.0 - Sw - So;
    
    REAL waterdensity = 0, oildensity= 0, gasdensity= 0;
    REAL dwaterdensitydP= 0, doildensitydP= 0, dgasdensitydP= 0;
    
    REAL waterviscosity= 0, oilviscosity= 0;
    REAL dwaterviscositydPw= 0, doilviscositydPo= 0;
    
    REAL krw= 0, kro= 0;
    REAL dkrwdSw= 0, dkrodSo= 0;
    
    
    this->PhaseDensities();
    waterdensity        = fWaterDensity[0];
    dwaterdensitydP     = fWaterDensity[1];
    oildensity          = fOilDensity[0];
    doildensitydP       = fOilDensity[1];
    gasdensity          = fGasDensity[0];
    dgasdensitydP       = fGasDensity[1];
    
    //    this->fFluidmodeldata->WaterViscosity(fWaterPressure[0], waterviscosity, dwaterviscositydPw);
    //    this->fFluidmodeldata->OilViscosity(fOilPressure[0], oilviscosity, doilviscositydPo);
    //    this->fFluidmodeldata->GasViscosity(fGasPressure[0], gasviscosity, dgasviscositydPg);
    //
    //    this->fPetrophysicdata->Krw(Sw, krw, dkrwdSw);
    //    this->fPetrophysicdata->Kro(1.0-Sw, kro, dkrodSo); // here appears the two-phase dependence
    //    this->fPetrophysicdata->Krg(1.0-(1.0-Sw)-Sw, krg, dkrgdSg);   // here appears the two-phase dependence
    
    // Computing fluid mobilities and derivatives for two-phase
    
    fWaterMobility[0] = waterdensity * krw / waterviscosity;
    fWaterMobility[1] = krw*(dwaterdensitydP/waterviscosity) - krw*((waterdensity * dwaterviscositydPw)/(waterviscosity*waterviscosity));
    fWaterMobility[2] = 1.0*dkrwdSw*(waterdensity/waterviscosity);
    fWaterMobility[3] = (-0.0)*dkrwdSw*(waterdensity/waterviscosity);   // here appears the two-phase dependence
    
    fOilMobility[0] = oildensity * kro / oilviscosity;
    fOilMobility[1] = kro*(doildensitydP/oilviscosity) - kro*((oildensity * doilviscositydPo)/(oilviscosity*oilviscosity));
    fOilMobility[2] = (-1.0)*dkrodSo*(oildensity/oilviscosity);         // here appears the two-phase dependence
    fOilMobility[3] = 0.0*dkrodSo*(oildensity/oilviscosity);
    
    fGasMobility[0] = 0.0;
    fGasMobility[1] = 0.0;
    fGasMobility[2] = 0.0;
    fGasMobility[3] = 0.0;
    
    //    // Computing fluid mobilities and derivatives for three-phase
    //
    //    fWaterMobility[0] = waterdensity * krw / waterviscosity;
    //    fWaterMobility[1] = krw*(dwaterdensitydP/waterviscosity) - krw*((waterdensity * dwaterviscositydPw)/(waterviscosity*waterviscosity));
    //    fWaterMobility[2] = dkrwdSw*(waterdensity/waterviscosity);
    //    fWaterMobility[3] = 0.0;
    //
    //    fOilMobility[0] = oildensity * kro / oilviscosity;
    //    fOilMobility[1] = kro*(doildensitydP/oilviscosity) - kro*((oildensity * doilviscositydPo)/(oilviscosity*oilviscosity));
    //    fOilMobility[2] = 0.0;
    //    fOilMobility[3] = dkrodSo*(oildensity/oilviscosity);
    //
    //    fGasMobility[0] = gasdensity * krg / gasviscosity;
    //    fGasMobility[1] = krg*(dgasdensitydP/gasviscosity) - krw*((gasdensity * dgasviscositydPg)/(gasviscosity*gasviscosity));
    //    fGasMobility[2] = (dkrgdSg*(-1.0)*(gasdensity/gasviscosity));
    //    fGasMobility[3] = (dkrgdSg*(-1.0)*(gasdensity/gasviscosity));
    
}

void TPZAxiSymmetricDarcyFlow::TotalDensity()
{
    REAL Sw, So, Sg;
    Sw = fWaterSaturation;
    So = fOilSaturation;
    Sg = 1.0 - Sw - So;
    
    this->PhaseDensities();
    
    fTotalDensity[0] = fWaterDensity[0] * Sw + fOilDensity[0] * So + fGasDensity[0] * (1.0 - Sw - So);
    fTotalDensity[1] = fWaterDensity[1] * Sw + fOilDensity[1] * So + fGasDensity[1] * (1.0 - Sw - So);
    fTotalDensity[2] = fWaterDensity[0] - fGasDensity[0];
    fTotalDensity[3] = fOilDensity[0]   - fGasDensity[0];
    
}

void TPZAxiSymmetricDarcyFlow::TotalMobility()
{
    this->PhaseMobilities();
    
    fTotalMobility[0] = fWaterMobility[0] + fOilMobility[0] + fGasMobility[0];
    fTotalMobility[1] = fWaterMobility[1] + fOilMobility[1] + fGasMobility[1];
    fTotalMobility[2] = fWaterMobility[2] + fOilMobility[2] + fGasMobility[2];
    fTotalMobility[3] = fWaterMobility[3] + fOilMobility[3] + fGasMobility[3];
    
    
}

void TPZAxiSymmetricDarcyFlow::PhaseFractionalFlows()
{
    this->TotalMobility();
    
    fFWater[0] = (fWaterMobility[0] / fTotalMobility[0]);
    fFWater[1] = (fWaterMobility[1] / fTotalMobility[0]) - ((fWaterMobility[0] * fTotalMobility[1])/(fTotalMobility[0]*fTotalMobility[0])) ;
    fFWater[2] = (fWaterMobility[2] / fTotalMobility[0]) - ((fWaterMobility[0] * fTotalMobility[2])/(fTotalMobility[0]*fTotalMobility[0])) ;
    fFWater[3] = (fWaterMobility[3] / fTotalMobility[0]) - ((fWaterMobility[0] * fTotalMobility[3])/(fTotalMobility[0]*fTotalMobility[0])) ;
    
    fFOil[0] = (fOilMobility[0] / fTotalMobility[0]);
    fFOil[1] = (fOilMobility[1] / fTotalMobility[0]) - ((fOilMobility[0] * fTotalMobility[1])/(fTotalMobility[0]*fTotalMobility[0])) ;
    fFOil[2] = (fOilMobility[2] / fTotalMobility[0]) - ((fOilMobility[0] * fTotalMobility[2])/(fTotalMobility[0]*fTotalMobility[0])) ;
    fFOil[3] = (fOilMobility[3] / fTotalMobility[0]) - ((fOilMobility[0] * fTotalMobility[3])/(fTotalMobility[0]*fTotalMobility[0])) ;
    
    fFGas[0] = (fGasMobility[0] / fTotalMobility[0]);
    fFGas[1] = (fGasMobility[1] / fTotalMobility[0]) - ((fGasMobility[0] * fTotalMobility[1])/(fTotalMobility[0]*fTotalMobility[0])) ;
    fFGas[2] = (fGasMobility[2] / fTotalMobility[0]) - ((fGasMobility[0] * fTotalMobility[2])/(fTotalMobility[0]*fTotalMobility[0])) ;
    fFGas[3] = (fGasMobility[3] / fTotalMobility[0]) - ((fGasMobility[0] * fTotalMobility[3])/(fTotalMobility[0]*fTotalMobility[0])) ;
    
}



int TPZAxiSymmetricDarcyFlow::ClassId() const{
    return -6378;
}

// -------------------------------------------------------------------------------------------

void TPZAxiSymmetricDarcyFlow::Write(TPZStream &buf, int withclassid) const{
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
}

// -------------------------------------------------------------------------------------------

void TPZAxiSymmetricDarcyFlow::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}



// System Properties

void TPZAxiSymmetricDarcyFlow::ComputeProperties(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZManVector<REAL> > & props){
    
    int nphases = fnstate_vars - 1;
    int n = fnstate_vars + 2;
    int p       = 1;
    int s_alpha = 2;
    int props_num = 10;
    //    int s_beta  = 3;
    
    TPZManVector<REAL,4> state_vars(fnstate_vars+1,0.0);
    props.Resize(props_num);
    
    TPZManVector<REAL> rho_alpha(n,0.0);
    TPZManVector<REAL> mu_alpha(n,0.0);
    TPZManVector<REAL> kr_alpha(n,0.0);
    TPZManVector<REAL> l_alpha(n,0.0);
    TPZManVector<REAL> f_alpha(n,0.0);
    TPZManVector<REAL> Pc_alpha(n,0.0);
    TPZManVector<REAL> vars_alpha(n-1,0.0);
    
    TPZManVector<REAL> rho_beta(n,0.0);
    TPZManVector<REAL> mu_beta(n,0.0);
    TPZManVector<REAL> kr_beta(n,0.0);
    TPZManVector<REAL> l_beta(n,0.0);
    TPZManVector<REAL> f_beta(n,0.0);
    TPZManVector<REAL> Pc_beta(n,0.0);
    TPZManVector<REAL> vars_beta(n-1,0.0);
    
    TPZManVector<REAL> rho_gamma(n,0.0);
    TPZManVector<REAL> mu_gamma(n,0.0);
    TPZManVector<REAL> kr_gamma(n,0.0);
    TPZManVector<REAL> l_gamma(n,0.0);
    TPZManVector<REAL> f_gamma(n,0.0);
    TPZManVector<REAL> Pc_gamma(n,0.0);
    TPZManVector<REAL> vars_gamma(n-1,0.0);
    
    // Strong formulation constitutive functions
    TPZManVector<REAL> rho(n,0.0);
    TPZManVector<REAL> rhof(n,0.0);
    TPZManVector<REAL> l(n,0.0);
    
    
    switch (nphases) {
        case 1:
        {
            // Getting state variables
            REAL P = datavec[p].sol[0][0];
            
            // Compute P_alpha
            state_vars[1] = P;
            REAL P_alpha = P;
            state_vars[1] = P_alpha;
            
            fluid_alpha->Density(rho_alpha, state_vars);
            fluid_alpha->Viscosity(mu_alpha, state_vars);
            
            // Rho computation
            rho = rho_alpha;
            
            // Rhof computation
            rhof = rho_alpha;
            
            // lambda computation
            l[0] = rho_alpha[0]/mu_alpha[0];
            l[2] = rho_alpha[2]/mu_alpha[0] - (rho_alpha[0]*mu_alpha[2])/(mu_alpha[0]*mu_alpha[0]);
            
        }
            break;
        case 2:
        {
            
            // Getting state variables
            REAL P          = datavec[p].sol[0][0];
            REAL S_alpha    = datavec[s_alpha].sol[0][0];
            REAL S_beta     = 1.0 - S_alpha;
            
            state_vars[1] = P;
            state_vars[2] = S_alpha;
            
            vars_alpha[1] = P;
            vars_alpha[2] = S_alpha;

            vars_beta[1] = P;
            vars_beta[2] = S_beta;
            
            // Compute P_alpha  and  P_beta
            fluid_beta->Pc(Pc_beta, state_vars);
            vars_alpha[1]   = P - S_alpha   * Pc_beta[1];
            vars_beta[1]    = P + S_beta    * Pc_beta[1];
            
            
            
            fluid_alpha->Density(rho_alpha, vars_alpha);
            fluid_alpha->Viscosity(mu_alpha, vars_alpha);
            fluid_alpha->Kr(kr_alpha, vars_alpha);
            
            
            fluid_beta->Density(rho_beta, vars_beta);
            fluid_beta->Viscosity(mu_beta, vars_beta);
            fluid_beta->Kr(kr_beta, vars_beta);
            
            // Rho computation
            rho[0] = rho_alpha[0] * S_alpha + rho_beta[0] * S_beta;
            rho[1] = rho_alpha[1] * S_alpha + rho_beta[1] * S_beta;
            rho[2] = rho_alpha[2] * S_alpha + rho_beta[2] * S_beta;
            rho[3] = rho_alpha[3] * S_alpha + rho_beta[3] * S_beta;
            
            // Two-Phase restriction
            kr_beta[3] *= -1.0;
            
            // Fluid Mobilities
            
            l_alpha[0] = kr_alpha[0]*(rho_alpha[0]/mu_alpha[0]);
            l_alpha[1] = 0.0;
            l_alpha[2] = kr_alpha[0]*(rho_alpha[2]/mu_alpha[0]-(rho_alpha[0]*mu_alpha[2])/(mu_alpha[0]*mu_alpha[0]));
            l_alpha[3] = kr_alpha[3]*(rho_alpha[0]/mu_alpha[0]);
            
            l_beta[0] = kr_beta[0]*(rho_beta[0]/mu_beta[0]);
            l_beta[1] = 0.0;
            l_beta[2] = kr_beta[0]*(rho_beta[2]/mu_beta[0]-(rho_beta[0]*mu_beta[2])/(mu_beta[0]*mu_beta[0]));
            l_beta[3] = kr_beta[3]*(rho_beta[0]/mu_beta[0]);
            
            // lambda computation
            l[0] = l_alpha[0] + l_beta[0];
            l[1] = l_alpha[1] + l_beta[1];
            l[2] = l_alpha[2] + l_beta[2];
            l[3] = l_alpha[3] + l_beta[3];
            
            // Fractional flows
            f_alpha[0] = l_alpha[0]/l[0];
            f_alpha[1] = l_alpha[1]/l[0]-(l_alpha[0]*l[1])/(l[0]*l[0]);
            f_alpha[2] = l_alpha[2]/l[0]-(l_alpha[0]*l[2])/(l[0]*l[0]);
            f_alpha[3] = l_alpha[3]/l[0]-(l_alpha[0]*l[3])/(l[0]*l[0]);
            
            f_beta[0] = l_beta[0]/l[0];
            f_beta[1] = l_beta[1]/l[0]-(l_beta[0]*l[1])/(l[0]*l[0]);
            f_beta[2] = l_beta[2]/l[0]-(l_beta[0]*l[2])/(l[0]*l[0]);
            f_beta[3] = l_beta[3]/l[0]-(l_beta[0]*l[3])/(l[0]*l[0]);
            
            // Rhof computation
            rhof[0] = f_alpha[0] * rho_alpha[0] +  f_beta[0] * rho_beta[0];
            rhof[1] = f_alpha[0] * rho_alpha[1] + f_alpha[1] * rho_alpha[0]  + f_beta[0] * rho_beta[1] + f_beta[1] * rho_beta[0];
            rhof[2] = f_alpha[0] * rho_alpha[2] + f_alpha[2] * rho_alpha[0]  + f_beta[0] * rho_beta[2] + f_beta[2] * rho_beta[0];
            rhof[3] = f_alpha[0] * rho_alpha[3] + f_alpha[3] * rho_alpha[0]  + f_beta[0] * rho_beta[3] + f_beta[3] * rho_beta[0];
            
        }
            break;
        case 3:
        {
            std::cout<< "3-phases system not implemented"<< std::endl;
            DebugStop();
            
        }
            break;
        default:
        {
            std::cout<< "4-phases system not implemented"<< std::endl;
            DebugStop();
        }
            break;
    }
    
    // Returning needed relationships
    
    props[0] = rho_alpha;
    props[1] = rho_beta;
    props[2] = rho_gamma;
    props[3] = f_alpha;
    props[4] = f_beta;
    props[5] = f_gamma;
    props[6] = l;
    props[7] = rho;
    props[8] = rhof;
    props[9] = Pc_beta;
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::ComputeProperties(TPZManVector<REAL> &state_vars, TPZVec<TPZManVector<REAL> > & props){
    
    int nphases = fnstate_vars - 1;
    int n = fnstate_vars + 2;
    int props_num = 10;
//    int p       = 1;
//    int s_alpha = 2;
    //    int s_beta  = 3;
    
    props.Resize(props_num);
    
    TPZManVector<REAL> rho_alpha(n,0.0);
    TPZManVector<REAL> mu_alpha(n,0.0);
    TPZManVector<REAL> kr_alpha(n,0.0);
    TPZManVector<REAL> l_alpha(n,0.0);
    TPZManVector<REAL> f_alpha(n,0.0);
    TPZManVector<REAL> Pc_alpha(n,0.0);
    TPZManVector<REAL> vars_alpha = state_vars;
    
    TPZManVector<REAL> rho_beta(n,0.0);
    TPZManVector<REAL> mu_beta(n,0.0);
    TPZManVector<REAL> kr_beta(n,0.0);
    TPZManVector<REAL> l_beta(n,0.0);
    TPZManVector<REAL> f_beta(n,0.0);
    TPZManVector<REAL> Pc_beta(n,0.0);
    TPZManVector<REAL> vars_beta = state_vars;
    
    TPZManVector<REAL> rho_gamma(n,0.0);
    TPZManVector<REAL> mu_gamma(n,0.0);
    TPZManVector<REAL> kr_gamma(n,0.0);
    TPZManVector<REAL> l_gamma(n,0.0);
    TPZManVector<REAL> f_gamma(n,0.0);
    TPZManVector<REAL> Pc_gamma(n,0.0);
    TPZManVector<REAL> vars_gamma = state_vars;
    
    // Strong formulation constitutive functions
    TPZManVector<REAL> rho(n,0.0);
    TPZManVector<REAL> rhof(n,0.0);
    TPZManVector<REAL> l(n,0.0);
    
    
    switch (nphases) {
        case 1:
        {
            // Getting state variables
            REAL P = state_vars[1];
            
            // Compute P_alpha
            REAL P_alpha = P;
            state_vars[1] = P_alpha;
            
            fluid_alpha->Density(rho_alpha, state_vars);
            fluid_alpha->Viscosity(mu_alpha, state_vars);
            
            // Rho computation
            rho = rho_alpha;
            
            // Rhof computation
            rhof = rho_alpha;
            
            // lambda computation
            l[0] = rho_alpha[0]/mu_alpha[0];
            l[2] = rho_alpha[2]/mu_alpha[0] - (rho_alpha[0]*mu_alpha[2])/(mu_alpha[0]*mu_alpha[0]);
            
        }
            break;
        case 2:
        {
            
            // Getting state variables
            REAL P          = state_vars[1];
            REAL S_alpha    = state_vars[2];
            REAL S_beta     = 1.0 - S_alpha;
            
            // Compute P_alpha  and  P_beta
            fluid_beta->Pc(Pc_beta, state_vars);
            vars_alpha[1]   = P - S_alpha   * Pc_beta[1];
            vars_beta[1]    = P + S_beta    * Pc_beta[1];
            
            vars_alpha[2] = S_alpha;
            vars_beta[2] = S_beta;
            
            fluid_alpha->Density(rho_alpha, vars_alpha);
            fluid_alpha->Viscosity(mu_alpha, vars_alpha);
            fluid_alpha->Kr(kr_alpha, vars_alpha);
            
            
            fluid_beta->Density(rho_beta, vars_beta);
            fluid_beta->Viscosity(mu_beta, vars_beta);
            fluid_beta->Kr(kr_beta, vars_beta);
            
            // Rho computation
            rho[0] = rho_alpha[0] * S_alpha + rho_beta[0] * S_beta;
            rho[1] = rho_alpha[1] * S_alpha + rho_beta[1] * S_beta;
            rho[2] = rho_alpha[2] * S_alpha + rho_beta[2] * S_beta;
            rho[3] = rho_alpha[3] * S_alpha + rho_beta[3] * S_beta;
            
            // Two-Phase restriction
            kr_beta[3] *= -1.0;
            
            // Fluid Mobilities
            
            l_alpha[0] = kr_alpha[0]*(rho_alpha[0]/mu_alpha[0]);
            l_alpha[1] = 0.0;
            l_alpha[2] = kr_alpha[0]*(rho_alpha[2]/mu_alpha[0]-(rho_alpha[0]*mu_alpha[2])/(mu_alpha[0]*mu_alpha[0]));
            l_alpha[3] = kr_alpha[3]*(rho_alpha[0]/mu_alpha[0]);
            
            l_beta[0] = kr_beta[0]*(rho_beta[0]/mu_beta[0]);
            l_beta[1] = 0.0;
            l_beta[2] = kr_beta[0]*(rho_beta[2]/mu_beta[0]-(rho_beta[0]*mu_beta[2])/(mu_beta[0]*mu_beta[0]));
            l_beta[3] = kr_beta[3]*(rho_beta[0]/mu_beta[0]);
            
            // lambda computation
            l[0] = l_alpha[0] + l_beta[0];
            l[1] = l_alpha[1] + l_beta[1];
            l[2] = l_alpha[2] + l_beta[2];
            l[3] = l_alpha[3] + l_beta[3];
            
            // Fractional flows
            f_alpha[0] = l_alpha[0]/l[0];
            f_alpha[1] = l_alpha[1]/l[0]-(l_alpha[0]*l[1])/(l[0]*l[0]);
            f_alpha[2] = l_alpha[2]/l[0]-(l_alpha[0]*l[2])/(l[0]*l[0]);
            f_alpha[3] = l_alpha[3]/l[0]-(l_alpha[0]*l[3])/(l[0]*l[0]);
            
            f_beta[0] = l_beta[0]/l[0];
            f_beta[1] = l_beta[1]/l[0]-(l_beta[0]*l[1])/(l[0]*l[0]);
            f_beta[2] = l_beta[2]/l[0]-(l_beta[0]*l[2])/(l[0]*l[0]);
            f_beta[3] = l_beta[3]/l[0]-(l_beta[0]*l[3])/(l[0]*l[0]);
            
            // Rhof computation
            rhof[0] = f_alpha[0] * rho_alpha[0] +  f_beta[0] * rho_beta[0];
            rhof[1] = f_alpha[0] * rho_alpha[1] + f_alpha[1] * rho_alpha[0]  + f_beta[0] * rho_beta[1] + f_beta[1] * rho_beta[0];
            rhof[2] = f_alpha[0] * rho_alpha[2] + f_alpha[2] * rho_alpha[0]  + f_beta[0] * rho_beta[2] + f_beta[2] * rho_beta[0];
            rhof[3] = f_alpha[0] * rho_alpha[3] + f_alpha[3] * rho_alpha[0]  + f_beta[0] * rho_beta[3] + f_beta[3] * rho_beta[0];
            
        }
            break;
        case 3:
        {
            std::cout<< "3-phases system not implemented"<< std::endl;
            DebugStop();
            
        }
            break;
        default:
        {
            std::cout<< "4-phases system not implemented"<< std::endl;
            DebugStop();
        }
            break;
    }
    
    // Returning needed relationships
    
    props[0] = rho_alpha;
    props[1] = rho_beta;
    props[2] = rho_gamma;
    props[3] = f_alpha;
    props[4] = f_beta;
    props[5] = f_gamma;
    props[6] = l;
    props[7] = rho;
    props[8] = rhof;
    props[9] = Pc_beta;
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::Rho_alpha(TPZVec<REAL> P_alpha){
    
}

void TPZAxiSymmetricDarcyFlow::Rho_beta(TPZVec<REAL> P_beta){
    
    
    
}

void TPZAxiSymmetricDarcyFlow::Rho_gamma(TPZVec<REAL> P_gamma){
    
    
    
}

