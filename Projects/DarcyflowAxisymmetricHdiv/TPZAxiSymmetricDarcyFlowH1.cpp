/*
 *  TPZAxiSymmetricDarcyFlow.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZAxiSymmetricDarcyFlowH1.h"
#include "pzbndcond.h"
#include "pzaxestools.h"

TPZAxiSymmetricDarcyFlowH1::TPZAxiSymmetricDarcyFlowH1() : TPZDiscontinuousGalerkin()
{
    fReservoirdata=NULL;
}

TPZAxiSymmetricDarcyFlowH1::TPZAxiSymmetricDarcyFlowH1(int matid) : TPZDiscontinuousGalerkin(matid)
{
    fReservoirdata=NULL;
}


TPZAxiSymmetricDarcyFlowH1::TPZAxiSymmetricDarcyFlowH1(const TPZAxiSymmetricDarcyFlowH1 &mat) : TPZDiscontinuousGalerkin(mat)
{
    fReservoirdata = mat.fReservoirdata;
}

TPZAxiSymmetricDarcyFlowH1::~TPZAxiSymmetricDarcyFlowH1()
{
    
}

void TPZAxiSymmetricDarcyFlowH1::FillDataRequirements(TPZMaterialData &data)
{
        data.SetAllRequirements(false);
        data.fNeedsSol = true;
}

void TPZAxiSymmetricDarcyFlowH1::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
{

    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;

}

void TPZAxiSymmetricDarcyFlowH1::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

int TPZAxiSymmetricDarcyFlowH1::VariableIndex(const std::string &name) {
    if (!strcmp("Pressure", name.c_str())) return 0;
    if (!strcmp("Velocity", name.c_str())) return 1;
    if (!strcmp("Density", name.c_str())) return 2;
    if (!strcmp("Porosity", name.c_str())) return 3;
    if (!strcmp("DivofVeclocity", name.c_str())) return 4;
    std::cout  << " Var index not implemented " << std::endl;
    DebugStop();
    return 0;
}

int TPZAxiSymmetricDarcyFlowH1::NSolutionVariables(int var) {
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
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return 0;
}

void TPZAxiSymmetricDarcyFlowH1::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout) {
    

    
    TPZVec<REAL> P = data.sol[0];
    TPZFMatrix<STATE> dPdx = data.dsol[0];

    
    TPZFNMatrix<660> GradofP;
    TPZAxesTools<REAL>::Axes2XYZ(dPdx, GradofP, data.axes);
    
    // Computing normal flux
    // Getting required Data
    TPZFMatrix<STATE> K = fReservoirdata->Kabsolute();
    STATE visosity, porosity, density, labmda;
    STATE dvisositydp, dporositydp, ddensitydp, dlambdadp;
    fReservoirdata->Viscosity(P[0], visosity, dvisositydp);
    fReservoirdata->Porosity(P[0], porosity, dporositydp);
    fReservoirdata->Density(P[0], density, ddensitydp);
    
    // Computing the fluid mobility
    labmda = density/visosity;
    dlambdadp = (ddensitydp/visosity) - ((1/(visosity*visosity))*(density*dvisositydp));
    
    TPZFMatrix<STATE> lambdaKGradofP(2,1);
    lambdaKGradofP(0,0) = -1.0 * (labmda)* (K(0,0)*GradofP(0,0) + K(0,1)*GradofP(1,0));
    lambdaKGradofP(1,0) = -1.0 * (labmda)* (K(1,0)*GradofP(0,0) + K(1,1)*GradofP(1,0));
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0:
        {
            Solout[0] = P[0];
        }
            break;
        case 1:
        {
            Solout[0] = lambdaKGradofP(0,0);
            Solout[1] = lambdaKGradofP(1,0);
        }
            break;
        case 2:
        {
            Solout[0] = density;
        }
            break;
        case 3:
        {
            Solout[0] = porosity;
        }
            break;
        case 4:
        {
            Solout[0] = 0.0;
        }
            break;
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
}

void TPZAxiSymmetricDarcyFlowH1::Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    
    // At each Integration Point.
    
    // Getting test and basis functions
    TPZFMatrix<STATE> PhiH1 = data.phi; // For H1   test functions
    TPZFMatrix<STATE> dPhiH1 = data.dphix; // Derivative For H1   test functions
    
    // Getting Linear combination of basis functions
    TPZManVector<STATE> P = data.sol[0];
    TPZFMatrix<STATE> dPdx = data.dsol[0];
    
    TPZFNMatrix<660> GradofP;
    TPZAxesTools<REAL>::Axes2XYZ(dPdx, GradofP, data.axes);
    
    TPZFNMatrix<660> GradofPhiH1;
    TPZAxesTools<REAL>::Axes2XYZ(dPhiH1, GradofPhiH1, data.axes);
    
    // Number of phis
    int nPhiH1 = PhiH1.Rows(); // For H1
    
    // Getting required Data
    TPZFMatrix<STATE> K = fReservoirdata->Kabsolute();
    STATE visosity, porosity, density, labmda;
    STATE dvisositydp, dporositydp, ddensitydp, dlambdadp;
    fReservoirdata->Viscosity(P[0], visosity, dvisositydp);
    fReservoirdata->Porosity(P[0], porosity, dporositydp);
    fReservoirdata->Density(P[0], density, ddensitydp);
    
    // Computing the fluid mobility
    labmda = density/visosity;
    dlambdadp = (ddensitydp/visosity) - ((1/(visosity*visosity))*(density*dvisositydp));
    
    // Defining local variables

    TPZFMatrix<STATE> dlambdadpKGradofP(2,1);
    dlambdadpKGradofP(0,0) = (dlambdadp)* (K(0,0)*GradofP(0,0) + K(0,1)*GradofP(1,0));
    dlambdadpKGradofP(1,0) = (dlambdadp)* (K(1,0)*GradofP(0,0) + K(1,1)*GradofP(1,0));
    
    TPZFMatrix<STATE> Gravity(2,1);
    Gravity(0,0) = -0.0;
    Gravity(1,0) = -0.0;
    
    TPZFMatrix<STATE> lambdaKGravity(2,1);
    lambdaKGravity(0,0) = (labmda)* (K(0,0)*Gravity(0,0) + K(0,1)*Gravity(1,0));
    lambdaKGravity(1,0) = (labmda)* (K(1,0)*Gravity(0,0) + K(1,1)*Gravity(1,0));
    
    TPZFMatrix<STATE> dlambdadpKGravity(2,1);
    dlambdadpKGravity(0,0) = (dlambdadp)* (K(0,0)*Gravity(0,0) + K(0,1)*Gravity(1,0));
    dlambdadpKGravity(1,0) = (dlambdadp)* (K(1,0)*Gravity(0,0) + K(1,1)*Gravity(1,0));
    
    TPZFMatrix<STATE> lambdaKGradofPhiH1(2,1);
    

    

    
    TPZFMatrix<STATE> iPhiHdiv(2,1);
    TPZFMatrix<STATE> GradofiPhiH1(2,1);
    TPZFMatrix<STATE> NormalVectTensorProGradofiPhiH1(2,2);
    
    
    for (int ip = 0; ip < nPhiH1; ip++)
    {
        for (int jp = 0; jp < nPhiH1; jp++)
        {
            lambdaKGradofPhiH1(0,0) = (labmda)* (K(0,0)*GradofPhiH1(0,jp) + K(0,1)*GradofPhiH1(1,jp));
            lambdaKGradofPhiH1(1,0) = (labmda)* (K(1,0)*GradofPhiH1(0,jp) + K(1,1)*GradofPhiH1(1,jp));
            
            ek(ip,jp) += weight * ( ( ( PhiH1(jp,0)*dlambdadpKGradofP(0,0) + lambdaKGradofPhiH1(0,0) ) * GradofPhiH1(0,ip) +  ( PhiH1(jp,0)* dlambdadpKGradofP(1,0) + lambdaKGradofPhiH1(1,0) ) * GradofPhiH1(1,ip))

                                   -  PhiH1(jp,0)*((( ddensitydp*lambdaKGravity(0,0)+dlambdadpKGravity(0,0))*GradofPhiH1(0,ip)+(ddensitydp*lambdaKGravity(1,0)+dlambdadpKGravity(1,0))* GradofPhiH1(1,ip))));
        }
    }
    

    this->Contribute(data,weight,ef);
    
}

void TPZAxiSymmetricDarcyFlowH1::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) {
    
    // At each Integration Point.
    
    // Getting test and basis functions
    TPZFMatrix<STATE> PhiH1 = data.phi; // For H1   test functions
    TPZFMatrix<STATE> dPhiH1 = data.dphix; // Derivative For H1   test functions
    
    // Getting Linear combination of basis functions
    TPZManVector<STATE> P = data.sol[0];
    TPZFMatrix<STATE> dPdx = data.dsol[0];
    
    TPZFNMatrix<660> GradofP;
    TPZAxesTools<REAL>::Axes2XYZ(dPdx, GradofP, data.axes);
    
    TPZFNMatrix<660> GradofPhiH1;
    TPZAxesTools<REAL>::Axes2XYZ(dPhiH1, GradofPhiH1, data.axes);
    
    // Number of phis
    int nPhiH1 = PhiH1.Rows(); // For H1
    
    // Getting required Data
    TPZFMatrix<STATE> K = fReservoirdata->Kabsolute();
    STATE visosity, porosity, density, labmda;
    STATE dvisositydp, dporositydp, ddensitydp, dlambdadp;
    fReservoirdata->Viscosity(P[0], visosity, dvisositydp);
    fReservoirdata->Porosity(P[0], porosity, dporositydp);
    fReservoirdata->Density(P[0], density, ddensitydp);
    
    // Computing the fluid mobility
    labmda = density/visosity;
    dlambdadp = (ddensitydp/visosity) - ((1/(visosity*visosity))*(density*dvisositydp));
    
    // Defining local variables
    
    TPZFMatrix<STATE> lambdaKGradofP(2,1);
    lambdaKGradofP(0,0) = (labmda)* (K(0,0)*GradofP(0,0) + K(0,1)*GradofP(1,0));
    lambdaKGradofP(1,0) = (labmda)* (K(1,0)*GradofP(0,0) + K(1,1)*GradofP(1,0));
    
    TPZFMatrix<STATE> Gravity(2,1);
    Gravity(0,0) = -0.0;
    Gravity(1,0) = -0.0;
    
    TPZFMatrix<STATE> lambdaKGravity(2,1);
    lambdaKGravity(0,0) = (labmda)* (K(0,0)*Gravity(0,0) + K(0,1)*Gravity(1,0));
    lambdaKGravity(1,0) = (labmda)* (K(1,0)*Gravity(0,0) + K(1,1)*Gravity(1,0));
    
    TPZFMatrix<STATE> iPhiHdiv(2,1);
    TPZFMatrix<STATE> GradofiPhiH1(2,1);
    TPZFMatrix<STATE> NormalVectTensorProGradofiPhiH1(2,2);

    
    for (int ip = 0; ip < nPhiH1; ip++)
    {
        /* $ \underset{\Omega_{e}}{\int}\left(K\;\lambda\nabla\left(P\right)\right)\cdot\nabla\phi\;\partial\Omega_{e}-\underset{\Omega_{e}}{\int}\left(K\;\lambda\nabla\left(\rho_{f}g\; z\right)\right)\cdot\nabla\phi\;\partial\Omega_{e} $ */
        
        ef(ip) += weight * ((lambdaKGradofP(0,0)*GradofPhiH1(0,ip)+lambdaKGradofP(1,0)*GradofPhiH1(1,ip)) -  density * (lambdaKGravity(0,0)*GradofPhiH1(0,ip)+lambdaKGravity(1,0)*GradofPhiH1(1,ip)));
        
    }
    
    
}

void TPZAxiSymmetricDarcyFlowH1::ContributeBC(TPZMaterialData &data, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
 

    
    // Getting test and basis functions
    TPZFMatrix<STATE> PhiH1 = data.phi; // For H1   test functions
    
    TPZFMatrix<STATE> dPhiH1 = data.dphix; // Derivative For H1   test functions
    TPZManVector<REAL,3> &normal = data.normal; // does It make sense? normal
    
    // Getting Linear combination of basis functions
    TPZManVector<STATE> P = data.sol[0];
    TPZFMatrix<STATE> dPdx = data.dsol[0];
    
    TPZFNMatrix<660> GradofP;
    TPZAxesTools<REAL>::Axes2XYZ(dPdx, GradofP, data.axes);
    
    // Computing normal flux
    // Getting required Data
    TPZFMatrix<STATE> K = fReservoirdata->Kabsolute();
    STATE visosity, porosity, density, labmda;
    STATE dvisositydp, dporositydp, ddensitydp, dlambdadp;
    fReservoirdata->Viscosity(P[0], visosity, dvisositydp);
    fReservoirdata->Porosity(P[0], porosity, dporositydp);
    fReservoirdata->Density(P[0], density, ddensitydp);
    
    // Computing the fluid mobility
    labmda = density/visosity;
    dlambdadp = (ddensitydp/visosity) - ((1/(visosity*visosity))*(density*dvisositydp));
    
    TPZFMatrix<STATE> lambdaKGradofP(2,1);
    lambdaKGradofP(0,0) = -1.0 * (labmda)* (K(0,0)*GradofP(0,0) + K(0,1)*GradofP(1,0));
    lambdaKGradofP(1,0) = -1.0 * (labmda)* (K(1,0)*GradofP(0,0) + K(1,1)*GradofP(1,0));
    
    STATE Qn = lambdaKGradofP(0,0)*normal[0] + lambdaKGradofP(1,0)*normal[1];
    
    // Number of phis
    int nPhiH1 = PhiH1.Rows();  // For Hdiv
    
    TPZFMatrix<STATE> iPhiHdiv(2,1);
    TPZFMatrix<STATE> jPhiHdiv(2,1);
    
    
    STATE Value[1];
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD
        {
            Value[0] = bc.Val2()(0,0);         //  Pressure
            for (int ip = 0; ip < nPhiH1; ip++)
            {
                ef(ip) += weight * gBigNumber * ( P[0] - Value[0]) * PhiH1(ip,0);
                
                for (int jp = 0; jp < nPhiH1; jp++)
                {
                    ek(ip,jp) += weight * gBigNumber * ( PhiH1(jp,0) ) * PhiH1(ip,0);
                }
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN
        {
            Value[0] = bc.Val2()(0,0);         //  NormalFlux
            
            for (int ip = 0; ip < nPhiH1; ip++)
            {
                
                ef(ip) += weight * (Value[0]) * PhiH1(ip,0);
                
            }
            
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            
            DebugStop();
        }
            break;
    }
    
    
}

int TPZAxiSymmetricDarcyFlowH1::ClassId() const {
    return -6379;
}

// -------------------------------------------------------------------------------------------

void TPZAxiSymmetricDarcyFlowH1::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    buf.Write(&fReservoirdata->fPref);
    buf.Write(&fReservoirdata->fKab(0,0));
    
}

// -------------------------------------------------------------------------------------------

void TPZAxiSymmetricDarcyFlowH1::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    buf.Read(&fReservoirdata->fPref);
    buf.Read(&fReservoirdata->fKab(0,0));
    
}
