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

TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow() : TPZMaterial()
{
    fReservoirdata=NULL;
}

TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow(int matid) : TPZMaterial(matid)
{
	fReservoirdata=NULL;
}


TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow(const TPZAxiSymmetricDarcyFlow &mat) : TPZMaterial(mat)
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
		default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
	}
    return 0;
}

void TPZAxiSymmetricDarcyFlow::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
	
//	TPZFMatrix<STATE> &dQsol = datavec[0].dsol[0];
//	TPZFMatrix<STATE> &dPsol = datavec[1].dsol[0];

    TPZVec<REAL> Q = datavec[0].sol[0];
    TPZVec<REAL> P = datavec[1].sol[0];
    
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
    TPZVec<STATE> Q = datavec[Qblock].sol[0];
    TPZVec<STATE> P = datavec[Pblock].sol[0];
    
    TPZFMatrix<STATE> dQdx = datavec[Qblock].dsol[0];
    TPZFMatrix<STATE> dPdx = datavec[Pblock].dsol[0];
    
    // Number of phis
    int nPhiHdiv = datavec[Qblock].fVecShapeIndex.NElements();  // For Hdiv
    int nPhiL2   = WL2.Rows();                                  // For L2
    
    // Getting required Data
    TPZFMatrix<STATE> KInverse = fReservoirdata->KabsoluteInv();
    STATE visosity, porosity, density;
    STATE dvisositydp, dporositydp, ddensitydp;
    fReservoirdata->Viscosity(P[0], visosity, dvisositydp);
    fReservoirdata->Porosity(P[0], porosity, dporositydp);
    fReservoirdata->Density(P[0], density, ddensitydp);
    
    // Defining local variables
    
    TPZFMatrix<STATE> etaKinvQi(2,1);

    
    TPZFMatrix<STATE> Gravity(2,1);
    Gravity(0,0) = -0.0;
    Gravity(1,0) = -9.8 * density;
    
    int ishapeindex, jshapeindex;
    int ivectorindex, jvectorindex;
    
    TPZFMatrix<STATE> iPhiHdiv(2,1);
    TPZFMatrix<STATE> jPhiHdiv(2,1);
    TPZFMatrix<STATE> GradofPhiH1(2,1);
    TPZFMatrix<STATE> NormalVectTensorProGradofPhiH1(2,2);
    STATE divofPhiHdiv = 0.0;
    
    for (int iq = 0; iq < nPhiHdiv; iq++)
    {
        ivectorindex = datavec[Qblock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[Qblock].fVecShapeIndex[iq].second;
        
        iPhiHdiv(0,0) = PhiH1(ishapeindex,0) * datavec[Qblock].fNormalVec(0,ivectorindex);
        iPhiHdiv(1,0) = PhiH1(ishapeindex,0) * datavec[Qblock].fNormalVec(1,ivectorindex);
        
        etaKinvQi(0,0) = (visosity) * ( KInverse(0,0)*iPhiHdiv(0,0) + KInverse(0,1)*iPhiHdiv(1,0) );
        etaKinvQi(1,0) = (visosity) * ( KInverse(1,0)*iPhiHdiv(0,0) + KInverse(1,1)*iPhiHdiv(1,0) );
        
        GradofPhiH1(0,0) = dPhiH1(0,ishapeindex)*datavec[Qblock].axes(0,0) + dPhiH1(1,ishapeindex)*datavec[Qblock].axes(1,0);
        GradofPhiH1(1,0) = dPhiH1(0,ishapeindex)*datavec[Qblock].axes(0,1) + dPhiH1(1,ishapeindex)*datavec[Qblock].axes(1,1);
        
        NormalVectTensorProGradofPhiH1(0,0) = datavec[Qblock].fNormalVec(0,ivectorindex)*GradofPhiH1(0,0);
        NormalVectTensorProGradofPhiH1(0,1) = datavec[Qblock].fNormalVec(0,ivectorindex)*GradofPhiH1(1,0);
        NormalVectTensorProGradofPhiH1(1,0) = datavec[Qblock].fNormalVec(2,ivectorindex)*GradofPhiH1(0,0);
        NormalVectTensorProGradofPhiH1(1,1) = datavec[Qblock].fNormalVec(1,ivectorindex)*GradofPhiH1(1,0);

        divofPhiHdiv = NormalVectTensorProGradofPhiH1(0,0) + NormalVectTensorProGradofPhiH1(1,1);
        
        for (int jq = 0; jq < nPhiHdiv; jq++)
        {
        
            jvectorindex = datavec[Qblock].fVecShapeIndex[jq].first;
            jshapeindex = datavec[Qblock].fVecShapeIndex[jq].second;
            
            jPhiHdiv(0,0) = PhiH1(jshapeindex,0) * datavec[Qblock].fNormalVec(0,jvectorindex);
            jPhiHdiv(1,0) = PhiH1(jshapeindex,0) * datavec[Qblock].fNormalVec(1,jvectorindex);
            

            
            ek(iq,jq) += weight * ((etaKinvQi(0,0)*jPhiHdiv(0,0) + etaKinvQi(1,0)*jPhiHdiv(1,0)));
        
        }
        
        for (int jp = 0; jp < nPhiL2; jp++)
        {

            ek(iq,jp + nPhiHdiv) += -1.0 * weight * (WL2(jp,0)*divofPhiHdiv);
            ek(jp + nPhiHdiv,iq) += -1.0 * weight * (WL2(jp,0)* density * divofPhiHdiv);
            
        }
            
    }

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
    TPZVec<STATE> Q = datavec[Qblock].sol[0];
    TPZVec<STATE> P = datavec[Pblock].sol[0];
    
    TPZFMatrix<STATE> dQdx = datavec[Qblock].dsol[0];
    TPZFMatrix<STATE> dPdx = datavec[Pblock].dsol[0];
    
    // Number of phis
    int nPhiHdiv = datavec[Qblock].fVecShapeIndex.NElements();  // For Hdiv
    int nPhiL2   = WL2.Rows();                                  // For L2
    
    // Getting required Data
    TPZFMatrix<STATE> KInverse = fReservoirdata->fKabinv;
    STATE visosity, porosity, density;
    STATE dvisositydp, dporositydp, ddensitydp;
    fReservoirdata->Viscosity(P[0], visosity, dvisositydp);
    fReservoirdata->Porosity(P[0], porosity, dporositydp);
    fReservoirdata->Density(P[0], density, ddensitydp);
    
    // Defining local variables
    
    TPZFMatrix<STATE> etaKinvQ(2,1);
    etaKinvQ(0,0) = (visosity)* (KInverse(0,0)*Q[0] + KInverse(0,1)*Q[1]);
    etaKinvQ(1,0) = (visosity)* (KInverse(1,0)*Q[0] + KInverse(1,1)*Q[1]);
    
    TPZFMatrix<STATE> Gravity(2,1);
    Gravity(0,0) = -0.0;
    Gravity(1,0) = -0.0 * density;
    
    int ishapeindex;
    int ivectorindex;
    
    TPZFMatrix<STATE> PhiHdiv(2,1);
    TPZFMatrix<STATE> GradofPhiH1(2,1);
    TPZFMatrix<STATE> NormalVectTensorProGradofPhiH1(2,2);
    STATE divofPhiHdiv = 0.0;
    STATE divofPhiHdiv2 = 0.0;
    STATE divofQ = 0.0;
    
    for (int iq = 0; iq < nPhiHdiv; iq++)
    {
        ivectorindex = datavec[Qblock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[Qblock].fVecShapeIndex[iq].second;
        
        PhiHdiv(0,0) = PhiH1(ishapeindex,0) * datavec[Qblock].fNormalVec(0,ivectorindex);
        PhiHdiv(1,0) = PhiH1(ishapeindex,0) * datavec[Qblock].fNormalVec(1,ivectorindex);
        
        GradofPhiH1(0,0) = dPhiH1(0,ishapeindex)*datavec[Qblock].axes(0,0) + dPhiH1(1,ishapeindex)*datavec[Qblock].axes(1,0);
        GradofPhiH1(1,0) = dPhiH1(0,ishapeindex)*datavec[Qblock].axes(0,1) + dPhiH1(1,ishapeindex)*datavec[Qblock].axes(1,1);
        
        NormalVectTensorProGradofPhiH1(0,0) = datavec[Qblock].fNormalVec(0,ivectorindex)*GradofPhiH1(0,0);
        NormalVectTensorProGradofPhiH1(0,1) = datavec[Qblock].fNormalVec(0,ivectorindex)*GradofPhiH1(1,0);
        NormalVectTensorProGradofPhiH1(1,0) = datavec[Qblock].fNormalVec(1,ivectorindex)*GradofPhiH1(0,0);
        NormalVectTensorProGradofPhiH1(1,1) = datavec[Qblock].fNormalVec(1,ivectorindex)*GradofPhiH1(1,0);
        
        divofPhiHdiv = NormalVectTensorProGradofPhiH1(0,0) + NormalVectTensorProGradofPhiH1(1,1);
        
        ef(iq) += weight * ((etaKinvQ(0,0)*PhiHdiv(0,0) + etaKinvQ(1,0)*PhiHdiv(1,0)) - (P[0]*divofPhiHdiv) - (Gravity(0,0)*PhiHdiv(0,0) + Gravity(1,0)*PhiHdiv(1,0)) );
        
    }
    
    for (int ip = 0; ip < nPhiL2; ip++)
    {
        divofQ = dQdx(0,0) + dQdx(1,1);
        
        ef(ip + nPhiHdiv) += weight * (density * divofQ * WL2(ip,0));
        
    }
    
}

void TPZAxiSymmetricDarcyFlow::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
	
    // At each Integration Point.
    
    int Qblock = 0;
    int Pblock = 1;
    
    // Getting test and basis functions
    TPZFMatrix<STATE> PhiH1 = datavec[Qblock].phi; // For H1   test functions
    TPZFMatrix<STATE> WL2   = datavec[Pblock].phi; // For HL2  test functions
    
    TPZFMatrix<STATE> dPhiH1 = datavec[Qblock].dphix; // Derivative For H1   test functions
    TPZFMatrix<STATE> dWL2   = datavec[Pblock].dphix; // Derivative For HL2  test functions
//    TPZManVector<REAL,3> &normal = datavec[Qblock].normal; // does It make sense? normal
    
    // Getting Linear combinations of basis functions
    TPZManVector<REAL,3> Q = datavec[Qblock].sol[0];
    TPZManVector<REAL,3> P = datavec[Pblock].sol[0];
    
    TPZFMatrix<STATE> dQdx = datavec[Qblock].dsol[0];
    TPZFMatrix<STATE> dPdx = datavec[Pblock].dsol[0];
    
    // Number of phis
    int nPhiHdiv = PhiH1.Rows();  // For Hdiv
    int nPhiL2   = WL2.Rows();                                  // For L2
    
    STATE Value[1];
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD
        {
            Value[0] = bc.Val2()(0,0);         //  Pressure
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
               
                ef(iq) += -1.0*weight * (Value[0] * PhiH1(iq,0));
                
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN
        {
            Value[0] = bc.Val2()(0,0);         //  NormalFlux
            
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                
                ef(iq) += gBigNumber * weight * (Q[0] - Value[0]) * PhiH1(iq,0);
                
                for (int jq = 0; jq < nPhiHdiv; jq++)
                {
                    ek(iq,jq) += gBigNumber * weight * (PhiH1(jq,0)) * PhiH1(iq,0);
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

int TPZAxiSymmetricDarcyFlow::ClassId() const {
	return -6378;
}

// -------------------------------------------------------------------------------------------

void TPZAxiSymmetricDarcyFlow::Write(TPZStream &buf, int withclassid) {

    TPZMaterial::Write(buf, withclassid);
    buf.Write(&fReservoirdata->fPref);
    buf.Write(&fReservoirdata->fKab(0,0));
	
}

// -------------------------------------------------------------------------------------------

void TPZAxiSymmetricDarcyFlow::Read(TPZStream &buf, void *context) {
    TPZMaterial::Read(buf, context);
    buf.Read(&fReservoirdata->fPref);
    buf.Read(&fReservoirdata->fKab(0,0));
	
}
