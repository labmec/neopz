/**
 * \file
 * @brief Contains implementations of the TPZViscoelastic methods.
 */
/*
 *  pzviscoelastic.cpp
 *  pos_processamento
 *
 *  Created by Nathan Shauer on 07/16/11.
 *  Copyright 2011 LabMec. All rights reserved.
 *
 */


#include "pzviscoelastic.h"

TPZViscoelastic::TPZViscoelastic(int id,REAL ElaE,REAL poissonE, REAL lambdaV, REAL muV, REAL alphaT, TPZVec <REAL> &force): TPZMatWithMem<TPZFMatrix, TPZElasticity3D>(id),flambdaV(lambdaV),fmuV(muV),falphaT(alphaT)
{
	flambdaE = (poissonE * ElaE)/((1+poissonE)*(1-2*poissonE));
	fmuE = ElaE/(2*(1+poissonE));
	REAL lambda = flambdaE-(falphaT*flambdaV)/(1+falphaT);
	REAL mu = fmuE -(falphaT*fmuV)/(1+falphaT); 
	fElaVE = mu*(3*lambda+2*mu)/(lambda+mu);
	fPoissonVE = lambda/(2*(lambda+mu));
	SetMaterialDataHook(fElaVE, fPoissonVE);
	SetForce(force);
	SetC();
}

void TPZViscoelastic::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix &ek,TPZFMatrix &ef)
{	
	int nstate = this->NStateVariables();
	int in;
	REAL val;
    if(fUpdateMem != 0)
    {
        UpdateQsi(data);
    }
    TPZFMatrix dphi = data.dphix;
    TPZFMatrix phi = data.phi;
    const int phr = phi.Rows();
    int index = data.intPtIndex;

    TPZFNMatrix<6>  qsi(6,1);
		int rows = this->MemItem(index).Rows();
    if (rows != 6) 
    {
        DebugStop(); //deve inicializar o qsi pelo SetDefaultMemory
    }
    else 
    {
        qsi = this->MemItem(index);
    }

    for(in = 0; in < phr; in++) 
    { 
        //in: test function index. First equation	
        val = 0.;
        val -= qsi(_XX_,0) * dphi(0,in); // |
        val -= qsi(_XY_,0) * dphi(1,in); // fk
        val -= qsi(_XZ_,0) * dphi(2,in); // |
        val *= 1./(1+falphaT);
        ef(in*nstate+0,0) += weight * val;
        
        //Second equation: fb and fk
        val = 0.;
        val -= qsi(_XY_,0) * dphi(0,in); // |
        val -= qsi(_YY_,0) * dphi(1,in); // fk
        val -= qsi(_YZ_,0) * dphi(2,in); // |
        val *= 1./(1+falphaT);
        ef(in*nstate+1,0) += weight * val;
        
        //third equation: fb and fk
        val = 0.;
        val -= qsi(_XZ_,0) * dphi(0,in); // |
        val -= qsi(_YZ_,0) * dphi(1,in); // fk
        val -= qsi(_ZZ_,0) * dphi(2,in); // |
        val *= 1./(1+falphaT);
        ef(in*nstate+2,0) += weight * val;
    }
    TPZMatWithMem<TPZFMatrix,TPZElasticity3D>::Contribute(data,weight,ek,ef);


}

void TPZViscoelastic::UpdateQsi(TPZMaterialData &data)
{
	TPZFNMatrix<6> qsi;
	TPZFNMatrix<6> qsin1(6,1,0.);
	TPZFNMatrix<9> DSolXYZ(3,3,0.);
	TPZFNMatrix<6> Strain(6,1);
	int index = data.intPtIndex;
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	DSolXYZ = data.dsol[0];
	qsi = MemItem(index);
	//data.axes.Multiply(data.dsol,DSolXYZ,1/*transpose*/);
	
	Strain.Redim(6,1);
	Strain(_XX_,0) = DSolXYZ(0,0);
	Strain(_YY_,0) = DSolXYZ(1,1);
	Strain(_ZZ_,0) = DSolXYZ(2,2);
	Strain(_XY_,0) = 0.5 * ( DSolXYZ(1,0) + DSolXYZ(0,1) );
	Strain(_XZ_,0) = 0.5 * ( DSolXYZ(2,0) + DSolXYZ(0,2) );
	Strain(_YZ_,0) = 0.5 * ( DSolXYZ(2,1) + DSolXYZ(1,2) );
	
	REAL tr;
	tr = Strain(_XX_,0)+Strain(_YY_,0)+Strain(_ZZ_,0);
	//Strain.Print("Strain: ");
	
	//REAL lambdaE,REAL muE, REAL lambdaV, REAL muV, REAL alphaT
	
	qsin1(_XX_,0) = (falphaT*(-(tr)*flambdaV - 2*Strain(_XX_,0)*fmuV) + qsi(_XX_,0))/(1 + falphaT);
	qsin1(_YY_,0) = (falphaT*(-(tr)*flambdaV - 2*Strain(_YY_,0)*fmuV) + qsi(_YY_,0))/(1 + falphaT);
	qsin1(_ZZ_,0) = (falphaT*(-(tr)*flambdaV - 2*Strain(_ZZ_,0)*fmuV) + qsi(_ZZ_,0))/(1 + falphaT);
	qsin1(_XY_,0) = (-2*falphaT*Strain(_XY_,0)*fmuV + qsi(_XY_,0))/(1 + falphaT);
	qsin1(_XZ_,0) = (-2*falphaT*Strain(_XZ_,0)*fmuV + qsi(_XZ_,0))/(1 + falphaT);
	qsin1(_YZ_,0) = (-2*falphaT*Strain(_YZ_,0)*fmuV + qsi(_YZ_,0))/(1 + falphaT);
	
	//qsin1.Print("qsin1",std::cout);
	MemItem(index) = qsin1;		
}


int TPZViscoelastic::VariableIndex(const std::string &name)
{
	if(!strcmp("Displacement",name.c_str()))  return TPZElasticity3D::EDisplacement;
	if(!strcmp("state",name.c_str()))  return TPZElasticity3D::EDisplacement;
	if(!strcmp("DisplacementX",name.c_str()))  return TPZElasticity3D::EDisplacementX;
	if(!strcmp("DisplacementY",name.c_str()))  return TPZElasticity3D::EDisplacementY;
	if(!strcmp("DisplacementZ",name.c_str()))  return TPZElasticity3D::EDisplacementZ;
	if(!strcmp("PrincipalStrain", name.c_str()))  return TPZElasticity3D::EPrincipalStrain;
	if(!strcmp("ViscoStressX",name.c_str()))  return TPZViscoelastic::EViscoStressX;
	if(!strcmp("ViscoStressY",name.c_str()))  return TPZViscoelastic::EViscoStressY;
	if(!strcmp("ViscoStressZ",name.c_str()))  return TPZViscoelastic::EViscoStressZ;
	return -1;
}

int TPZViscoelastic::NSolutionVariables(int var)
{
	if(var == TPZElasticity3D::EDisplacement)        return 3;
	if(var == TPZElasticity3D::EDisplacementX)       return 1;
	if(var == TPZElasticity3D::EDisplacementY)       return 1;
	if(var == TPZElasticity3D::EDisplacementZ)       return 1;
	if(var == TPZElasticity3D::EPrincipalStrain)     return 3;
	if(var == TPZViscoelastic::EViscoStressX)        return 1;
	if(var == TPZViscoelastic::EViscoStressY)        return 1;
	if(var == TPZViscoelastic::EViscoStressZ)        return 1;
	PZError << "TPZViscoelastic::NSolutionVariables Error\n";
	return -1;
}

void TPZViscoelastic::ComputeStressTensor(TPZFMatrix &Stress, TPZMaterialData &data)
{
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZFNMatrix<9> Dsol = data.dsol[0];
	TPZMatWithMem<TPZFMatrix,TPZElasticity3D>::ComputeStressTensor(Stress,Dsol);
	TPZFNMatrix<6> qsi;
	const int index = data.intPtIndex;
	qsi = MemItem(index);
	qsi(_XX_,0)*= 1/(1+falphaT);
	qsi(_XY_,0)*= 1/(1+falphaT);
	qsi(_XZ_,0)*= 1/(1+falphaT);
	qsi(_YY_,0)*= 1/(1+falphaT);
	qsi(_YZ_,0)*= 1/(1+falphaT);
	qsi(_ZZ_,0)*= 1/(1+falphaT);
	
	Stress(0,0) += qsi(_XX_,0);
	Stress(1,1) += qsi(_YY_,0);
	Stress(2,2) += qsi(_ZZ_,0);
	Stress(0,1) += qsi(_XY_,0);
	Stress(0,2) += qsi(_XZ_,0);
	Stress(1,0) += qsi(_XY_,0);
	Stress(1,2) += qsi(_YZ_,0);
	Stress(2,0) += qsi(_XZ_,0);
	Stress(2,1) += qsi(_YZ_,0);
}

void TPZViscoelastic::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	if(var == TPZElasticity3D::EDisplacement){
		int i;
		for(i = 0; i < 3; i++){
			TPZVec<REAL> Sol(data.sol[0]); 
			Solout[i] = Sol[i];
		}//for
		return;
	}//TPZElasticity3D::EDisplacement
	
	if(var == TPZElasticity3D::EDisplacementX){
		//    int i;
		TPZVec<REAL> Sol(data.sol[0]); 
		Solout[0] = Sol[0];
		return;
	}//TPZElasticity3D::EDisplacementX
	
	if(var == TPZElasticity3D::EDisplacementY){
		//    int i;
		TPZVec<REAL> Sol(data.sol[0]); 
		Solout[0] = Sol[1];
		return;
	}//TPZElasticity3D::EDisplacementY  
	
	if(var == TPZElasticity3D::EDisplacementZ){
		//    int i;
		TPZVec<REAL> Sol(data.sol[0]); 
		Solout[0] = Sol[2];
		return;
	}//TPZElasticity3D::EDisplacementZ  

	
	if(var == TPZElasticity3D::EPrincipalStrain){
		TPZFNMatrix<9> StrainTensor(3,3);
		TPZFMatrix DSol(data.dsol[0]);
		TPZMatWithMem<TPZFMatrix,TPZElasticity3D>::ComputeStrainTensor(StrainTensor, DSol);
		int numiterations = 1000;
		REAL tol = TPZElasticity3D::gTolerance;
		bool result = StrainTensor.SolveEigenvaluesJacobi(numiterations, tol, &Solout);
#ifdef DEBUG    
		if (result == false){
			PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << tol << std::endl;
		}
#endif
	}//TPZElasticity3D::EPrincipalStrain
	
	
	if(var == TPZViscoelastic::EViscoStressX){
		TPZFNMatrix<9> Stress(3,3);
		this->ComputeStressTensor(Stress, data);
		Solout[0] = Stress(0,0);
		return;
	}
	if(var == TPZViscoelastic::EViscoStressY){
		TPZFMatrix Stress(3,3);
		this->ComputeStressTensor(Stress, data);
		Solout[0] = Stress(1,1);
		return;
	}
	if(var == TPZViscoelastic::EViscoStressZ){
		TPZFMatrix Stress(3,3);
		this->ComputeStressTensor(Stress, data);
		Solout[0] = Stress(2,2);
		return;
	}
}

