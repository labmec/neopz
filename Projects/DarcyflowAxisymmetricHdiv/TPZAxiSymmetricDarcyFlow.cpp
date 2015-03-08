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

TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow() : TPZDiscontinuousGalerkin()
{
    fReservoirdata=NULL;
}

TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow(int matid) : TPZDiscontinuousGalerkin(matid)
{
    fReservoirdata=NULL;
}


TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow(const TPZAxiSymmetricDarcyFlow &mat) : TPZDiscontinuousGalerkin(mat)
{
    fReservoirdata = mat.fReservoirdata;
}

TPZAxiSymmetricDarcyFlow::~TPZAxiSymmetricDarcyFlow()
{
    
}

void TPZAxiSymmetricDarcyFlow::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

void TPZAxiSymmetricDarcyFlow::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
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
    if (!strcmp("Pressure", name.c_str())) return 0;
    if (!strcmp("Velocity", name.c_str())) return 1;
    if (!strcmp("Density", name.c_str())) return 2;
    if (!strcmp("Porosity", name.c_str())) return 3;
    if (!strcmp("DivofVeclocity", name.c_str())) return 4;
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
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return 0;
}

void TPZAxiSymmetricDarcyFlow::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int Qblock = 0;
    int Pblock = 1;
    
    TPZVec<REAL> Q = datavec[Qblock].sol[0];
    TPZVec<REAL> P = datavec[Pblock].sol[0];
    
    TPZFMatrix<STATE> dQdx = datavec[Qblock].dsol[0];
    TPZFMatrix<STATE> dPdx = datavec[Pblock].dsol[0];
    
    STATE visosity, porosity, density;
    STATE dvisositydp, dporositydp, ddensitydp;
    fReservoirdata->Viscosity(P[0], visosity, dvisositydp);
    fReservoirdata->Porosity(P[0], porosity, dporositydp);
    fReservoirdata->Density(P[0], density, ddensitydp);
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0:
        {
            Solout[0] = P[0];
        }
            break;
        case 1:
        {
            Solout[0] = Q[0];
            Solout[1] = Q[1];
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
            Solout[0] = dQdx(0,0) + dQdx(1,1) + dQdx(2,2);
        }
            break;
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
}

void TPZAxiSymmetricDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    
    // At each Integration Point.
    
    int Qblock = 0;
    int Pblock = 1;
    
    // Getting test and basis functions
    TPZFMatrix<STATE> PhiH1 = datavec[Qblock].phi; // For H1   test functions
    TPZFMatrix<STATE> WL2   = datavec[Pblock].phi; // For HL2  test functions
    
    TPZFMatrix<STATE> dPhiH1 = datavec[Qblock].dphix; // Derivative For H1   test functions
    TPZFMatrix<STATE> dWL2   = datavec[Pblock].dphix; // Derivative For HL2  test functions
    
    // Getting Linear combinations of basis functions
    TPZManVector<STATE> Q = datavec[Qblock].sol[0];
    TPZManVector<STATE> P = datavec[Pblock].sol[0];
    
    TPZFMatrix<STATE> dQdaxes = datavec[Qblock].dsol[0];
    TPZFMatrix<STATE> dPdaxes = datavec[Pblock].dsol[0];
    
    // Computing Gradients of basis functions
    TPZFNMatrix<660> GradofWL2;
    TPZAxesTools<REAL>::Axes2XYZ(dWL2, GradofWL2, datavec[Pblock].axes);
    
    // Computing Gradient of P solution
    TPZFNMatrix<660> GradofP;
    TPZAxesTools<REAL>::Axes2XYZ(dPdaxes, GradofP, datavec[Pblock].axes);
    
    // Number of phis
    int nPhiHdiv = datavec[Qblock].fVecShapeIndex.NElements();  // For Hdiv
    int nPhiL2   = WL2.Rows();                                  // For L2
    
    // Getting required Data
    TPZFMatrix<STATE> KInverse = fReservoirdata->KabsoluteInv();
    STATE visosity, porosity, density, labmda;
    STATE dvisositydp, dporositydp, ddensitydp, dlambdadp;
    fReservoirdata->Viscosity(P[0], visosity, dvisositydp);
    fReservoirdata->Porosity(P[0], porosity, dporositydp);
    fReservoirdata->Density(P[0], density, ddensitydp);
    
    // Computing the fluid mobility
    labmda = density/visosity;
    dlambdadp = (ddensitydp/visosity) - ((1/(visosity*visosity))*(density*dvisositydp));
    
    // Defining local variables
    
    TPZFMatrix<STATE> lambdainvKinvjPhiHdiv(2,1);
    TPZFMatrix<STATE> lambdainvKinvQ(2,1);
    lambdainvKinvQ(0,0) = (1.0/labmda)* (KInverse(0,0)*Q[0] + KInverse(0,1)*Q[1]);
    lambdainvKinvQ(1,0) = (1.0/labmda)* (KInverse(1,0)*Q[0] + KInverse(1,1)*Q[1]);
    
    TPZFMatrix<STATE> Gravity(2,1);
    Gravity(0,0) = -0.0;
    Gravity(1,0) = -0.0 * density;
    
    int ishapeindex, jshapeindex;
    int ivectorindex, jvectorindex;
    
    TPZFMatrix<STATE> iPhiHdiv(2,1);
    TPZFMatrix<STATE> jPhiHdiv(2,1);
    TPZFMatrix<STATE> GradofiPhiH1(2,1);
    TPZFMatrix<STATE> NormalVectTensorProGradofiPhiH1(2,2);
    STATE divofiPhiHdiv = 0.0;
    TPZFMatrix<STATE> GradofjPhiH1(2,1);
    TPZFMatrix<STATE> NormalVectTensorProGradofjPhiH1(2,2);
    STATE divofjPhiHdiv = 0.0;
    
    for (int iq = 0; iq < nPhiHdiv; iq++)
    {
        ivectorindex = datavec[Qblock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[Qblock].fVecShapeIndex[iq].second;
        
        iPhiHdiv(0,0) = PhiH1(ishapeindex,0) * datavec[Qblock].fNormalVec(0,ivectorindex);
        iPhiHdiv(1,0) = PhiH1(ishapeindex,0) * datavec[Qblock].fNormalVec(1,ivectorindex);
        
        GradofiPhiH1(0,0) = dPhiH1(0,ishapeindex)*datavec[Qblock].axes(0,0) + dPhiH1(1,ishapeindex)*datavec[Qblock].axes(1,0);
        GradofiPhiH1(1,0) = dPhiH1(0,ishapeindex)*datavec[Qblock].axes(0,1) + dPhiH1(1,ishapeindex)*datavec[Qblock].axes(1,1);
        
        NormalVectTensorProGradofiPhiH1(0,0) = datavec[Qblock].fNormalVec(0,ivectorindex)*GradofiPhiH1(0,0);
        NormalVectTensorProGradofiPhiH1(0,1) = datavec[Qblock].fNormalVec(0,ivectorindex)*GradofiPhiH1(1,0);
        NormalVectTensorProGradofiPhiH1(1,0) = datavec[Qblock].fNormalVec(1,ivectorindex)*GradofiPhiH1(0,0);
        NormalVectTensorProGradofiPhiH1(1,1) = datavec[Qblock].fNormalVec(1,ivectorindex)*GradofiPhiH1(1,0);
        
        divofiPhiHdiv = NormalVectTensorProGradofiPhiH1(0,0) + NormalVectTensorProGradofiPhiH1(1,1);
        
        lambdainvKinvjPhiHdiv(0,0) = (1.0/labmda) * ( KInverse(0,0)*iPhiHdiv(0,0) + KInverse(0,1)*iPhiHdiv(1,0) );
        lambdainvKinvjPhiHdiv(1,0) = (1.0/labmda) * ( KInverse(1,0)*iPhiHdiv(0,0) + KInverse(1,1)*iPhiHdiv(1,0) );
        
        
        for (int jq = 0; jq < nPhiHdiv; jq++)
        {
            
            jvectorindex = datavec[Qblock].fVecShapeIndex[jq].first;
            jshapeindex = datavec[Qblock].fVecShapeIndex[jq].second;
            
            jPhiHdiv(0,0) = PhiH1(jshapeindex,0) * datavec[Qblock].fNormalVec(0,jvectorindex);
            jPhiHdiv(1,0) = PhiH1(jshapeindex,0) * datavec[Qblock].fNormalVec(1,jvectorindex);
            
            ek(iq,jq) += weight * (lambdainvKinvjPhiHdiv(0,0)*jPhiHdiv(0,0) + lambdainvKinvjPhiHdiv(1,0)*jPhiHdiv(1,0));
            
        }
        
        for (int jp = 0; jp < nPhiL2; jp++)
        {
            ek(iq,jp + nPhiHdiv) += weight * (- (dlambdadp/labmda)*(lambdainvKinvQ(0,0)*jPhiHdiv(0,0) + lambdainvKinvQ(1,0)*jPhiHdiv(1,0)) - divofiPhiHdiv - ddensitydp * (Gravity(0,0)*iPhiHdiv(0,0) + Gravity(1,0)*iPhiHdiv(1,0))) * WL2(jp,0);
            
        }
        
    }
    
    for (int ip = 0; ip < nPhiL2; ip++)
    {
        
        for (int jq = 0; jq < nPhiHdiv; jq++)
        {
            
            jvectorindex = datavec[Qblock].fVecShapeIndex[jq].first;
            jshapeindex = datavec[Qblock].fVecShapeIndex[jq].second;
            
            jPhiHdiv(0,0) = PhiH1(jshapeindex,0) * datavec[Qblock].fNormalVec(0,jvectorindex);
            jPhiHdiv(1,0) = PhiH1(jshapeindex,0) * datavec[Qblock].fNormalVec(1,jvectorindex);
            
            
            GradofjPhiH1(0,0) = dPhiH1(0,jshapeindex)*datavec[Qblock].axes(0,0) + dPhiH1(1,jshapeindex)*datavec[Qblock].axes(1,0);
            GradofjPhiH1(1,0) = dPhiH1(0,jshapeindex)*datavec[Qblock].axes(0,1) + dPhiH1(1,jshapeindex)*datavec[Qblock].axes(1,1);
            
            NormalVectTensorProGradofjPhiH1(0,0) = datavec[Qblock].fNormalVec(0,jvectorindex)*GradofjPhiH1(0,0);
            NormalVectTensorProGradofjPhiH1(0,1) = datavec[Qblock].fNormalVec(0,jvectorindex)*GradofjPhiH1(1,0);
            NormalVectTensorProGradofjPhiH1(1,0) = datavec[Qblock].fNormalVec(1,jvectorindex)*GradofjPhiH1(0,0);
            NormalVectTensorProGradofjPhiH1(1,1) = datavec[Qblock].fNormalVec(1,jvectorindex)*GradofjPhiH1(1,0);
            
            divofjPhiHdiv = NormalVectTensorProGradofjPhiH1(0,0) + NormalVectTensorProGradofjPhiH1(1,1);
            
            ek(ip + nPhiHdiv,jq) += -1.0 * weight * divofjPhiHdiv  * WL2(ip,0);
            
        }
        
    }
    
    this->Contribute(datavec,weight,ef);
    
}

void TPZAxiSymmetricDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef) {
    
    // At each Integration Point.
    
    int Qblock = 0;
    int Pblock = 1;
    
    // Getting test and basis functions
    TPZFMatrix<STATE> PhiH1 = datavec[Qblock].phi; // For H1   test functions
    TPZFMatrix<STATE> WL2   = datavec[Pblock].phi; // For HL2  test functions
    
    TPZFMatrix<STATE> dPhiH1 = datavec[Qblock].dphix; // Derivative For H1   test functions
    TPZFMatrix<STATE> dWL2   = datavec[Pblock].dphix; // Derivative For HL2  test functions
    
    // Getting Linear combinations of basis functions
    TPZManVector<STATE> Q = datavec[Qblock].sol[0];
    TPZManVector<STATE> P = datavec[Pblock].sol[0];
    
    TPZFMatrix<STATE> dQdx = datavec[Qblock].dsol[0];
    TPZFMatrix<STATE> dPdx = datavec[Pblock].dsol[0];
    
    TPZFNMatrix<660> GradofP;
    TPZAxesTools<REAL>::Axes2XYZ(dPdx, GradofP, datavec[Pblock].axes);
    
    // Number of phis
    int nPhiHdiv = datavec[Qblock].fVecShapeIndex.NElements();  // For Hdiv
    int nPhiL2   = WL2.Rows();                                  // For L2
    
    // Getting required Data
    TPZFMatrix<STATE> KInverse = fReservoirdata->KabsoluteInv();
    STATE visosity, porosity, density, labmda;
    STATE dvisositydp, dporositydp, ddensitydp, dlambdadp;
    fReservoirdata->Viscosity(P[0], visosity, dvisositydp);
    fReservoirdata->Porosity(P[0], porosity, dporositydp);
    fReservoirdata->Density(P[0], density, ddensitydp);
    
    // Computing the fluid mobility
    labmda = density/visosity;
    dlambdadp = (ddensitydp/visosity) - ((1/(visosity*visosity))*(density*dvisositydp));
    
    // Defining local variables
    
    TPZFMatrix<STATE> oneoverlambdaKinvQ(2,1);
    oneoverlambdaKinvQ(0,0) = (1.0/labmda)* (KInverse(0,0)*Q[0] + KInverse(0,1)*Q[1]);
    oneoverlambdaKinvQ(1,0) = (1.0/labmda)* (KInverse(1,0)*Q[0] + KInverse(1,1)*Q[1]);
    
    TPZFMatrix<STATE> Gravity(2,1);
    Gravity(0,0) = -0.0;
    Gravity(1,0) = -0.0;
    
    int ishapeindex;
    int ivectorindex;
    
    TPZFMatrix<STATE> iPhiHdiv(2,1);
    TPZFMatrix<STATE> GradofiPhiH1(2,1);
    TPZFMatrix<STATE> NormalVectTensorProGradofiPhiH1(2,2);
    STATE divofiPhiHdiv = 0.0;
    STATE divofQ = 0.0;
    
    for (int iq = 0; iq < nPhiHdiv; iq++)
    {
        ivectorindex = datavec[Qblock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[Qblock].fVecShapeIndex[iq].second;
        
        iPhiHdiv(0,0) = PhiH1(ishapeindex,0) * datavec[Qblock].fNormalVec(0,ivectorindex);
        iPhiHdiv(1,0) = PhiH1(ishapeindex,0) * datavec[Qblock].fNormalVec(1,ivectorindex);
        
        GradofiPhiH1(0,0) = dPhiH1(0,ishapeindex)*datavec[Qblock].axes(0,0) + dPhiH1(1,ishapeindex)*datavec[Qblock].axes(1,0);
        GradofiPhiH1(1,0) = dPhiH1(0,ishapeindex)*datavec[Qblock].axes(0,1) + dPhiH1(1,ishapeindex)*datavec[Qblock].axes(1,1);
        
        NormalVectTensorProGradofiPhiH1(0,0) = datavec[Qblock].fNormalVec(0,ivectorindex)*GradofiPhiH1(0,0);
        NormalVectTensorProGradofiPhiH1(0,1) = datavec[Qblock].fNormalVec(0,ivectorindex)*GradofiPhiH1(1,0);
        NormalVectTensorProGradofiPhiH1(1,0) = datavec[Qblock].fNormalVec(1,ivectorindex)*GradofiPhiH1(0,0);
        NormalVectTensorProGradofiPhiH1(1,1) = datavec[Qblock].fNormalVec(1,ivectorindex)*GradofiPhiH1(1,0);
        
        divofiPhiHdiv = NormalVectTensorProGradofiPhiH1(0,0) + NormalVectTensorProGradofiPhiH1(1,1);
        
        /* $ \underset{\Omega_{e}}{\int}\left(K\lambda\right)^{-1}\mathbf{q}\cdot\mathbf{v}\;\partial\Omega_{e}-\underset{\Omega_{e}}{\int}P\; div\left(\mathbf{v}\right)\partial\Omega-\underset{\Omega_{e}}{\int}\nabla\left(\rho_{f}g\; z\right)\cdot\mathbf{v} $ */
        
        ef(iq) += weight * ((oneoverlambdaKinvQ(0,0)*iPhiHdiv(0,0) + oneoverlambdaKinvQ(1,0)*iPhiHdiv(1,0)) - (P[0]*divofiPhiHdiv) - density * (Gravity(0,0)*iPhiHdiv(0,0) + Gravity(1,0)*iPhiHdiv(1,0)) );
        
    }
    
    divofQ = dQdx(0,0) + dQdx(1,1) + dQdx(2,2);
    
    /* $ - \underset{\Omega}{\int}w\; div\left(\mathbf{q}\right)\partial\Omega $ */
    
    for (int ip = 0; ip < nPhiL2; ip++)
    {
        ef(ip + nPhiHdiv) += -1.0 * weight * (divofQ * WL2(ip,0));
        
    }
    
}

void TPZAxiSymmetricDarcyFlow::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
    // At each Integration Point.
    
    int Qblock = 0;
    int Pblock = 1;
    
    // Getting test and basis functions
    TPZFMatrix<STATE> PhiH1 = dataleftvec[Qblock].phi; // For H1   test functions
    TPZFMatrix<STATE> WL2   = dataleftvec[Pblock].phi; // For HL2  test functions
    
    TPZFMatrix<STATE> dPhiH1 = dataleftvec[Qblock].dphix; // Derivative For H1   test functions
    TPZFMatrix<STATE> dWL2   = dataleftvec[Pblock].dphix; // Derivative For HL2  test functions
    TPZManVector<REAL,3> &normal = data.normal; // does It make sense? normal
    
    // Getting Linear combinations of basis functions
    TPZManVector<STATE> Q = dataleftvec[Qblock].sol[0];
    TPZManVector<STATE> P = dataleftvec[Pblock].sol[0];
    
    // Computing normal flux
    STATE Qn = Q[0]*normal[0] + Q[1]*normal[1];
    
    TPZFMatrix<STATE> dQdx = dataleftvec[Qblock].dsol[0];
    TPZFMatrix<STATE> dPdx = dataleftvec[Pblock].dsol[0];
    
    // Number of phis
    int nPhiHdiv = dataleftvec[Qblock].fVecShapeIndex.NElements();  // For Hdiv
    int nPhiL2   = WL2.Rows();                                  // For L2
    
    int ishapeindex, jshapeindex;
    int ivectorindex, jvectorindex;
    
    TPZFMatrix<STATE> iPhiHdiv(2,1);
    TPZFMatrix<STATE> jPhiHdiv(2,1);
    
    STATE iphinormal;
    STATE jphinormal;
    
    
    STATE Value[1];
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD
        {
            Value[0] = bc.Val2()(0,0);         //  Pressure
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                
                ivectorindex = dataleftvec[Qblock].fVecShapeIndex[iq].first;
                ishapeindex = dataleftvec[Qblock].fVecShapeIndex[iq].second;
                
                iPhiHdiv(0,0) = PhiH1(ishapeindex,0) * dataleftvec[Qblock].fNormalVec(0,ivectorindex);
                iPhiHdiv(1,0) = PhiH1(ishapeindex,0) * dataleftvec[Qblock].fNormalVec(1,ivectorindex);
                
                iphinormal = normal[0]*iPhiHdiv(0,0) + normal[1]*iPhiHdiv(1,0);
                
                ef(iq) += weight * ( Value[0] * iphinormal);
                
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN
        {
            Value[0] = bc.Val2()(0,0);         //  NormalFlux
            
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                
                ivectorindex = dataleftvec[Qblock].fVecShapeIndex[iq].first;
                ishapeindex = dataleftvec[Qblock].fVecShapeIndex[iq].second;
                
                iPhiHdiv(0,0) = PhiH1(ishapeindex,0) * dataleftvec[Qblock].fNormalVec(0,ivectorindex);
                iPhiHdiv(1,0) = PhiH1(ishapeindex,0) * dataleftvec[Qblock].fNormalVec(1,ivectorindex);
                
                iphinormal = normal[0]*iPhiHdiv(0,0) + normal[1]*iPhiHdiv(1,0);
                
                ef(iq) += weight * (gBigNumber * (Qn - Value[0]) + P[0]) * iphinormal;
                
                for (int jq = 0; jq < nPhiHdiv; jq++)
                {
                    
                    jvectorindex = dataleftvec[Qblock].fVecShapeIndex[jq].first;
                    jshapeindex = dataleftvec[Qblock].fVecShapeIndex[jq].second;
                    
                    jPhiHdiv(0,0) = PhiH1(jshapeindex,0) * dataleftvec[Qblock].fNormalVec(0,jvectorindex);
                    jPhiHdiv(1,0) = PhiH1(jshapeindex,0) * dataleftvec[Qblock].fNormalVec(1,jvectorindex);
                    
                    jphinormal = normal[0]*jPhiHdiv(0,0) + normal[1]*jPhiHdiv(1,0);
                    
                    ek(iq,jq) += gBigNumber * weight * (jphinormal) * iphinormal;
                }
                
                for (int jp = 0; jp < nPhiL2; jp++)
                {
                    ek(iq, jp + nPhiHdiv) += weight * (WL2(jp,0)) * iphinormal;
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
    
}

void TPZAxiSymmetricDarcyFlow::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
    
    return;
    
}

int TPZAxiSymmetricDarcyFlow::ClassId() const {
    return -6378;
}

// -------------------------------------------------------------------------------------------

void TPZAxiSymmetricDarcyFlow::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    buf.Write(&fReservoirdata->fPref);
    buf.Write(&fReservoirdata->fKab(0,0));
    
}

// -------------------------------------------------------------------------------------------

void TPZAxiSymmetricDarcyFlow::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    buf.Read(&fReservoirdata->fPref);
    buf.Read(&fReservoirdata->fKab(0,0));
    
}
