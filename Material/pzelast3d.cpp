/**
 * \file
 * @brief Contains implementations of the TPZElasticity3D methods.
 */
//$Id: pzelast3d.cpp,v 1.17 2011-05-11 02:24:19 phil Exp $

#include "pzelast3d.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
//#include "pztempmat.h"
#include "pzmanvector.h"
#include <math.h>
#include <fstream>

REAL TPZElasticity3D::gTolerance = 1.e-11;

TPZElasticity3D::TPZElasticity3D(int nummat, REAL E, REAL poisson, TPZVec<REAL> &force) : TPZMaterial(nummat){
	this->fE = E;
	this->fPoisson = poisson;
#ifdef DEBUG
	if (force.NElements() != 3) PZError << __PRETTY_FUNCTION__ << " - error!" << std::endl;
#endif  
	int i;  
	this->fForce.Resize(3);
	for(i = 0; i < 3; i++) this->fForce[i] = force[i];
	//Default directions is {1,0,0}
	this->fPostProcessDirection.Resize(3);
	this->fPostProcessDirection.Fill(0.);
	this->fPostProcessDirection[0] = 1.;
	this->SetYieldingStress(1.);
#ifndef CODE1
	C1 = E / (2.+ 2.*poisson);
	C2 = E * poisson / (-1. + poisson + 2.*poisson*poisson);
	C3 = E * (poisson - 1.) / (-1. + poisson +2. * poisson * poisson);
#endif
	
}//method

TPZElasticity3D::TPZElasticity3D(int nummat): TPZMaterial(nummat),fE(0.), fPoisson(0.),fForce(3,0.),fPostProcessDirection(3,0.), fFy(0.)
{
    SetC();
}

TPZElasticity3D::TPZElasticity3D() : TPZMaterial(0),fE(0.), fPoisson(0.),fForce(3,0.),fPostProcessDirection(3,0.), 
fFy(0.){}

TPZElasticity3D::~TPZElasticity3D(){}

TPZElasticity3D::TPZElasticity3D(const TPZElasticity3D &cp) : TPZMaterial(cp), fE(cp.fE), fPoisson(cp.fPoisson),fForce(cp.fForce),fPostProcessDirection(cp.fPostProcessDirection), 
fFy(cp.fFy)
{
    SetC();
}

void TPZElasticity3D::Print(std::ostream & out){
	out << "\nTPZElasticity3D material:\n";
	out << "\tfE       = " << this->fE << std::endl;
	out << "\tfPoisson = " << this->fPoisson << std::endl;
	out << "\tfForce   = " << this->fForce << std::endl;
	out << "\tBase class print\n";
	TPZMaterial::Print(out);
	out << "End of TPZElasticity3D::Print\n";
}

void TPZElasticity3D::Contribute(TPZMaterialData &data,
                                 REAL weight,
                                 TPZFMatrix &ek,
                                 TPZFMatrix &ef){
	TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	TPZFMatrix &phi = data.phi;
	// TPZFMatrix &phiL = data.phil;
	// TPZFMatrix &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &jacinv = data.jacinv;
	
	
	
	const int phr = phi.Rows();
	if(this->fForcingFunction){
		this->fForcingFunction->Execute(x,this->fForce);
	}
	
	//this matrix will store {{dvdx*dudx, dvdx*dudy, dvdx*dudz},
	//{dvdy*dudx, dvdy*dudy, dvdy*dudz},
	//{dvdz*dudx, dvdz*dudy, dvdz*dudz}}
	
#ifdef CODE1
	TPZFNMatrix<9> Deriv(3,3);
	const REAL E  = this->fE;
	const REAL nu = this->fPoisson;
	const REAL C1 = E / (2.+ 2.*nu);
	const REAL C2 = E * nu / (-1. + nu + 2.*nu*nu);
	const REAL C3 = E * (nu - 1.) / (-1. + nu +2. * nu * nu);
	
	int in;
	for(in = 0; in < phr; in++) {
		int kd;
		for(kd = 0; kd < 3; kd++){
			ef(in*3+kd, 0) += weight* fForce[kd] * phi(in,0);
		}//kd
		REAL val;
		for( int jn = 0; jn < phr; jn++ ) {
			
			//Compute Deriv matrix
			for(int ud = 0; ud < 3; ud++){
				for(int vd = 0; vd < 3; vd++){
					Deriv(vd,ud) = dphi(vd,in)*dphi(ud,jn);
				}//ud
			}//vd
			
			//First equation Dot[Sigma1, gradV1]
			val = ( Deriv(1,1) + Deriv(2,2) ) * C1 + Deriv(0,0) * C3;
			ek(in*3+0,jn*3+0) += weight * val;
			
			val = Deriv(1,0) * C1 - Deriv(0,1) * C2;
			ek(in*3+0,jn*3+1) += weight * val;
			
			val = Deriv(2,0) * C1 - Deriv(0,2) * C2;
			ek(in*3+0,jn*3+2) += weight * val;
			
			//Second equation Dot[Sigma2, gradV2]
			val = Deriv(0,1) * C1 - Deriv(1,0) * C2;
			ek(in*3+1,jn*3+0) += weight * val;
			
			val = ( Deriv(0,0) + Deriv(2,2) ) * C1 + Deriv(1,1) * C3;
			ek(in*3+1,jn*3+1) += weight * val;
			
			val = Deriv(2,1) * C1 - Deriv(1,2) * C2;
			ek(in*3+1,jn*3+2) += weight * val;
			
			//Third equation Dot[Sigma3, gradV3]
			val = Deriv(0,2) * C1 - Deriv(2,0) * C2;
			ek(in*3+2,jn*3+0) += weight * val;
			
			val = Deriv(1,2) * C1 - Deriv(2,1) * C2;
			ek(in*3+2,jn*3+1) += weight * val;
			
			val = ( Deriv(0,0) + Deriv(1,1) ) * C1 + Deriv(2,2) * C3;
			ek(in*3+2,jn*3+2) += weight * val;
			
		}//jn
	}//in
#else
	REAL Deriv[3][3];
	int in;
	for(in = 0; in < phr; in++) {
		int kd;
		for(kd = 0; kd < 3; kd++){
			ef(in*3+kd, 0) += weight* fForce[kd] * phi(in,0);
		}//kd
		for( int jn = 0; jn < phr; jn++ ) {
			
			//Compute Deriv matrix
			for(int ud = 0; ud < 3; ud++){
				for(int vd = 0; vd < 3; vd++){
					Deriv[vd][ud] = dphi(vd,in)*dphi(ud,jn);
				}//ud
			}//vd
			
			//First equation Dot[Sigma1, gradV1]
			REAL *ptr1 = &ek(in*3,jn*3);
			/*ek(in*3+0,jn*3+0)*/ *ptr1++ += weight * (( Deriv[1][1] + Deriv[2][2] ) * C1 + Deriv[0][0] * C3);
			
			//Second equation Dot[Sigma2, gradV2]
			/*ek(in*3+1,jn*3+0)*/ *ptr1++ += weight * (Deriv[0][1] * C1 - Deriv[1][0] * C2);
			
			//Third equation Dot[Sigma3, gradV3]
			/*ek(in*3+2,jn*3+0)*/ *ptr1 += weight * (Deriv[0][2] * C1 - Deriv[2][0] * C2);
			
			REAL *ptr2 = &ek(in*3+0, jn*3+1);
			/*ek(in*3+0,jn*3+1)*/ *ptr2++ += weight * (Deriv[1][0] * C1 - Deriv[0][1] * C2);
			
			/*ek(in*3+1,jn*3+1)*/ *ptr2++ += weight * (( Deriv[0][0] + Deriv[2][2] ) * C1 + Deriv[1][1] * C3);
			
			/*ek(in*3+2,jn*3+1)*/ *ptr2 += weight * (Deriv[1][2] * C1 - Deriv[2][1] * C2);
			
			REAL *ptr3 = &ek(in*3+0, jn*3+2);
			/*ek(in*3+0,jn*3+2)*/ *ptr3++ += weight * (Deriv[2][0] * C1 - Deriv[0][2] * C2);
			
			/*ek(in*3+1,jn*3+2)*/ *ptr3++ += weight * (Deriv[2][1] * C1 - Deriv[1][2] * C2);
			
			/*ek(in*3+2,jn*3+2)*/ *ptr3 += weight * (( Deriv[0][0] + Deriv[1][1] ) * C1 + Deriv[2][2] * C3);
			
		}//jn
	}//in
	
#endif
#ifdef DEBUG   
	if ( !ek.VerifySymmetry( 1.e-8 ) ) PZError << __PRETTY_FUNCTION__ << "\nERROR - NON SYMMETRIC MATRIX" << std::endl;
#endif
}//method

void TPZElasticity3D::ContributeBC(TPZMaterialData &data,
                                   REAL weight,
                                   TPZFMatrix &ek,
                                   TPZFMatrix &ef,
                                   TPZBndCond &bc){
	// TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	TPZFMatrix &phi = data.phi;
	// TPZFMatrix &phiL = data.phil;
	// TPZFMatrix &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	// TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &jacinv = data.jacinv;
	
	
	const REAL BIGNUMBER  = 1.e12;
	
	const int phr = phi.Rows();
	int in,jn,idf,jdf;
	REAL v2[3];
	v2[0] = bc.Val2()(0,0);
	v2[1] = bc.Val2()(1,0);
	v2[2] = bc.Val2()(2,0);
	TPZFMatrix &v1 = bc.Val1();
    
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
	switch (bc.Type()) {
		case 0: // Dirichlet condition
			for(in = 0 ; in < phr; in++) {
				ef(3*in+0,0) += BIGNUMBER * v2[0] * phi(in,0) * weight;
				ef(3*in+1,0) += BIGNUMBER * v2[1] * phi(in,0) * weight;        
				ef(3*in+2,0) += BIGNUMBER * v2[2] * phi(in,0) * weight;        
				
				for (jn = 0 ; jn < phr; jn++) {
					ek(3*in+0,3*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
					ek(3*in+1,3*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
					ek(3*in+2,3*jn+2) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
				}//jn
			}//in
			break;
			
		case 1: // Neumann condition
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(3*in+0,0) += v2[0] * phi(in,0) * weight;
				ef(3*in+1,0) += v2[1] * phi(in,0) * weight;
				ef(3*in+2,0) += v2[2] * phi(in,0) * weight;
			}//in
			break;
		case 2: // Mixed condition
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(3*in+0,0) += v2[0] * phi(in,0) * weight;
				ef(3*in+1,0) += v2[1] * phi(in,0) * weight;
				ef(3*in+2,0) += v2[2] * phi(in,0) * weight;
				for(jn=0; jn<phi.Rows(); jn++)
				{
					for(idf=0; idf<3; idf++) for(jdf=0; jdf<3; jdf++)
					{
						ek(3*in+idf,3*jn+jdf) += bc.Val1()(idf,jdf);
					}
				}
			}//in
			break;
		case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
			for(in = 0 ; in < phr; in++) {
				ef(3*in+0,0) += BIGNUMBER * (0. - data.sol[0][0]) * v2[0] * phi(in,0) * weight;
				ef(3*in+1,0) += BIGNUMBER * (0. - data.sol[0][1]) * v2[1] * phi(in,0) * weight;        
				ef(3*in+2,0) += BIGNUMBER * (0. - data.sol[0][2]) * v2[2] * phi(in,0) * weight;        
				for (jn = 0 ; jn < phr; jn++) {
					ek(3*in+0,3*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[0];
					ek(3*in+1,3*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[1];
					ek(3*in+2,3*jn+2) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[2];
				}//jn
			}//in
			break;
			
		case 4: // stressField Neumann condition
			for(in = 0; in < 3; in ++)
				v2[in] = - ( v1(in,0) * data.normal[0] +
							v1(in,1) * data.normal[1] +
							v1(in,2) * data.normal[2] );
			// The normal vector points towards the neighbour. The negative sign is there to 
			// reflect the outward normal vector.
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(3*in+0,0) += v2[0] * phi(in,0) * weight;
				ef(3*in+1,0) += v2[1] * phi(in,0) * weight;
				ef(3*in+2,0) += v2[2] * phi(in,0) * weight;
				//cout << "normal:" << data.normal[0] << ' ' << data.normal[1] << ' ' << data.normal[2] << endl;
				//cout << "val2:  " << v2[0]          << ' ' << v2[1]          << ' ' << v2[2]          << endl;
			}
			break;
		default:
			PZError << "TPZElastitity3D::ContributeBC error - Wrong boundary condition type" << std::endl;
	}//switch
	
}//method

int TPZElasticity3D::VariableIndex(const std::string &name){
	if(!strcmp("Displacement",name.c_str()))  return TPZElasticity3D::EDisplacement;
	if(!strcmp("state",name.c_str()))  return TPZElasticity3D::EDisplacement;
	if(!strcmp("DisplacementX",name.c_str()))  return TPZElasticity3D::EDisplacementX;
	if(!strcmp("DisplacementY",name.c_str()))  return TPZElasticity3D::EDisplacementY;
	if(!strcmp("DisplacementZ",name.c_str()))  return TPZElasticity3D::EDisplacementZ;
	if(!strcmp("PrincipalStress", name.c_str()))  return TPZElasticity3D::EPrincipalStress;
	if(!strcmp("PrincipalStrain", name.c_str()))  return TPZElasticity3D::EPrincipalStrain;
	if(!strcmp("VonMises",    name.c_str()))  return TPZElasticity3D::EVonMisesStress;
	if(!strcmp("Stress",     name.c_str()))  return TPZElasticity3D::EStress;
	if(!strcmp("Strain",     name.c_str()))  return TPZElasticity3D::EStrain;
	if(!strcmp("Stress1",     name.c_str()))  return TPZElasticity3D::EStress1;
	if(!strcmp("Strain1",     name.c_str()))  return TPZElasticity3D::EStrain1;  
	if(!strcmp("NormalStress",name.c_str()))  return TPZElasticity3D::ENormalStress;
	if(!strcmp("NormalStrain",name.c_str()))  return TPZElasticity3D::ENormalStrain;
	if(!strcmp("StressX",name.c_str()))  return TPZElasticity3D::EStressX;
	if(!strcmp("StressY",name.c_str()))  return TPZElasticity3D::EStressY;
	if(!strcmp("StressZ",name.c_str()))  return TPZElasticity3D::EStressZ;
	PZError << "TPZElasticity3D::VariableIndex Error\n";
	return -1;
}

int TPZElasticity3D::NSolutionVariables(int var){
	if(var == TPZElasticity3D::EDisplacement)        return 3;
	if(var == TPZElasticity3D::EDisplacementX)       return 1;
	if(var == TPZElasticity3D::EDisplacementY)       return 1;
	if(var == TPZElasticity3D::EDisplacementZ)       return 1;
	if(var == TPZElasticity3D::EPrincipalStress)     return 3;
	if(var == TPZElasticity3D::EPrincipalStrain)     return 3;
	if(var == TPZElasticity3D::EPrincipalDirection1) return 3;
	if(var == TPZElasticity3D::EPrincipalDirection2) return 3;
	if(var == TPZElasticity3D::EPrincipalDirection3) return 3;
	if(var == TPZElasticity3D::EVonMisesStress)      return 1;  
	if(var == TPZElasticity3D::EStress)              return 3;  
	if(var == TPZElasticity3D::EStrain)              return 3;  
	if(var == TPZElasticity3D::EStrain1)             return 1;  
	if(var == TPZElasticity3D::EStress1)             return 1; 
	if(var == TPZElasticity3D::ENormalStress)        return 3;  
	if(var == TPZElasticity3D::ENormalStrain)        return 3; 
	if(var == TPZElasticity3D::EStressX)						 return 1; 
    if(var == TPZElasticity3D::EStressY)						 return 1;
    if(var == TPZElasticity3D::EStressZ)						 return 1;
	PZError << "TPZElasticity3D::NSolutionVariables Error\n";
	return -1;
}

void TPZElasticity3D::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){
	
	if(var == TPZElasticity3D::EDisplacement){
		int i;
		for(i = 0; i < 3; i++){
			Solout[i] = Sol[i];
		}//for
		return;
	}//TPZElasticity3D::EDisplacement
	
	if(var == TPZElasticity3D::EDisplacementX){
		//    int i;
		Solout[0] = Sol[0];
		return;
	}//TPZElasticity3D::EDisplacementX
	
	if(var == TPZElasticity3D::EDisplacementY){
		//    int i;
		Solout[0] = Sol[1];
		return;
	}//TPZElasticity3D::EDisplacementY  
	
	if(var == TPZElasticity3D::EDisplacementZ){
		//    int i;
		Solout[0] = Sol[2];
		return;
	}//TPZElasticity3D::EDisplacementZ  
	
	if(var == TPZElasticity3D::EPrincipalStress){
		TPZFNMatrix<9> StressTensor(3,3);
		this->ComputeStressTensor(StressTensor, DSol);
		int numiterations = 1000;
		REAL tol = TPZElasticity3D::gTolerance;
		bool result = StressTensor.SolveEigenvaluesJacobi(numiterations, tol, &Solout);
#ifdef DEBUG        
		if (result == false){
			PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << tol << std::endl;
		}
#endif
	}//TPZElasticity3D::EPrincipalStress
	
	if(var == TPZElasticity3D::EStress1){
		TPZFNMatrix<9> StressTensor(3,3);
		TPZManVector<REAL, 3> PrincipalStress(3);
		this->ComputeStressTensor(StressTensor, DSol);
		int numiterations = 1000;
		REAL tol = TPZElasticity3D::gTolerance;
		bool result = StressTensor.SolveEigenvaluesJacobi(numiterations, tol, &PrincipalStress);
		Solout[0] = PrincipalStress[0];
#ifdef DEBUG        
		if (result == false){
			PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << tol << std::endl;
		}
#endif
	}//TPZElasticity3D::EStress1  
	
	
	if(var == TPZElasticity3D::EPrincipalStrain){
		TPZFNMatrix<9> StrainTensor(3,3);
		this->ComputeStrainTensor(StrainTensor, DSol);
		int numiterations = 1000;
		REAL tol = TPZElasticity3D::gTolerance;
		bool result = StrainTensor.SolveEigenvaluesJacobi(numiterations, tol, &Solout);
#ifdef DEBUG    
		if (result == false){
			PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << tol << std::endl;
		}
#endif
	}//TPZElasticity3D::EPrincipalStrain
	
	if(var == TPZElasticity3D::EStrain1){
		TPZFNMatrix<9> StrainTensor(3,3);
		TPZManVector<REAL, 3> PrincipalStrain(3);
		this->ComputeStrainTensor(StrainTensor, DSol);
		int numiterations = 1000;
		REAL tol = TPZElasticity3D::gTolerance;
		bool result = StrainTensor.SolveEigenvaluesJacobi(numiterations, tol, &PrincipalStrain);
		Solout[0] = PrincipalStrain[0];
#ifdef DEBUG    
		if (result == false){
			PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << tol << std::endl;
		}
#endif
	}//TPZElasticity3D::EPrincipalStrain  
	
	if(var == TPZElasticity3D::EPrincipalDirection1){
		this->PrincipalDirection(DSol, Solout, 0);    
	}//TPZElasticity3D::EPrincipalDirection1
	
	if(var == TPZElasticity3D::EPrincipalDirection2){
		this->PrincipalDirection(DSol, Solout, 1);    
	}//TPZElasticity3D::EPrincipalDirection2
	
	if(var == TPZElasticity3D::EPrincipalDirection3){
		this->PrincipalDirection(DSol, Solout, 2);    
	}//TPZElasticity3D::EPrincipalDirection3    
	
	if(var == TPZElasticity3D::EVonMisesStress){
		TPZManVector<REAL,3> PrincipalStress(3);
		TPZFNMatrix<9> StressTensor(3,3);
		this->ComputeStressTensor(StressTensor, DSol);
		int numiterations = 1000;
		REAL tol = TPZElasticity3D::gTolerance;
		bool result = StressTensor.SolveEigenvaluesJacobi(numiterations, tol, &PrincipalStress);
#ifdef DEBUG        
		if (result == false){
			PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << tol << std::endl;
		}
#endif    
		
		Solout.Resize(1);
		Solout[0] = ( PrincipalStress[0] - PrincipalStress[1] ) * ( PrincipalStress[0] - PrincipalStress[1] ) 
		+ ( PrincipalStress[1] - PrincipalStress[2] ) * ( PrincipalStress[1] - PrincipalStress[2] )
		+ ( PrincipalStress[2] - PrincipalStress[0] ) * ( PrincipalStress[2] - PrincipalStress[0] );
		Solout[0] = Solout[0] / (2. * this->fFy * this->fFy);    
	}//TPZElasticity3D::EVonMisesStress
	
	if(var == TPZElasticity3D::EStress){
		TPZFNMatrix<6> Stress(6,1);
		this->ComputeStressVector(Stress, DSol);
		this->ApplyDirection(Stress, Solout);
		return;
	}//TPZElasticity3D::EStress
	
	if(var == TPZElasticity3D::EStrain){
		TPZFNMatrix<6> Strain(6,1);
		this->ComputeStrainVector(Strain, DSol);
		this->ApplyDirection(Strain, Solout);
		return;
	}//TPZElasticity3D::EStrain
	
	if(var == TPZElasticity3D::ENormalStress){
		TPZFNMatrix<6> Stress(6,1);
		this->ComputeStressVector(Stress, DSol);
		Solout[0] = Stress(0,0);
		Solout[1] = Stress(1,0);
		Solout[2] = Stress(2,0);
		return;
	}//TPZElasticity3D::ENormalStress
	
	if(var == TPZElasticity3D::ENormalStrain){
		TPZFNMatrix<6> Strain(6,1);
		this->ComputeStrainVector(Strain, DSol);
		Solout[0] = Strain(0,0);
		Solout[1] = Strain(1,0);
		Solout[2] = Strain(2,0);
		return;
	}//TPZElasticity3D::ENormalStrain
	
	if(var == TPZElasticity3D::EStressX){
		TPZFNMatrix<9> Stress(3,3);
		this->ComputeStressTensor(Stress, DSol);
		Solout[0] = Stress(0,0);
		return;
	}//TPZElasticity3D::EStressX
    
    if(var == TPZElasticity3D::EStressY){
		TPZFNMatrix<9> Stress(3,3);
		this->ComputeStressTensor(Stress, DSol);
		Solout[0] = Stress(1,1);
		return;
	}//TPZElasticity3D::EStressY
    
    if(var == TPZElasticity3D::EStressZ){
		TPZFNMatrix<9> Stress(3,3);
		this->ComputeStressTensor(Stress, DSol);
		Solout[0] = Stress(2,2);
		return;
	}//TPZElasticity3D::EStressZ
	
}//Solution

void TPZElasticity3D::Errors(TPZVec<REAL> &x,TPZVec<REAL> &u, TPZFMatrix &dudx, 
							 TPZFMatrix &axes, TPZVec<REAL> &flux,
							 TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values){
	int i, j;
	
	/** L2 norm */
	REAL L2 = 0.;
	for(i = 0; i < 3; i++) L2 += (u[i] - u_exact[i]) * (u[i] - u_exact[i]);
	
	/** H1 semi-norm */
	REAL SemiH1 = 0.;
	for(i = 0; i < 3; i++) for(j = 0; j < 3; j++) SemiH1 += (dudx(i,j) - du_exact(i,j)) * (dudx(i,j) - du_exact(i,j));
	
	/** H1 norm */
	REAL H1 = L2 + SemiH1;
	
	//values[1] : eror em norma L2
	values[1]  = L2;
	
	//values[2] : erro em semi norma H1
	values[2] = SemiH1;
	
	//values[0] : erro em norma H1 <=> norma Energia
	values[0]  = H1;
	
}

void TPZElasticity3D::ComputeStrainTensor(TPZFMatrix &Strain, TPZFMatrix &DSol) const{
	Strain.Redim(3,3);
	Strain(0,0) = DSol(0,0);
	Strain(0,1) = 0.5 * ( DSol(1,0) + DSol(0,1) );
	Strain(0,2) = 0.5 * ( DSol(2,0) + DSol(0,2) );
	Strain(1,0) = Strain(0,1);
	Strain(1,1) = DSol(1,1);
	Strain(1,2) = 0.5 * ( DSol(2,1) + DSol(1,2) );
	Strain(2,0) = Strain(0,2);
	Strain(2,1) = Strain(1,2);
	Strain(2,2) = DSol(2,2);
}

void TPZElasticity3D::ComputeStressTensor(TPZFMatrix &Stress, TPZMaterialData &data) const{
    ComputeStressTensor(Stress, data.dsol[0]);
}

void TPZElasticity3D::ComputeStressTensor(TPZFMatrix &Stress, TPZFMatrix &dsol) const{
	TPZFNMatrix<6> Vec(6,1);
	this->ComputeStressVector(Vec, dsol);
	
	Stress.Redim(3,3);
	Stress(0,0) = Vec(0,0);
	Stress(0,1) = Vec(3,0);
	Stress(0,2) = Vec(4,0);
	Stress(1,0) = Vec(3,0);
	Stress(1,1) = Vec(1,0);
	Stress(1,2) = Vec(5,0);
	Stress(2,0) = Vec(4,0);
	Stress(2,1) = Vec(5,0);
	Stress(2,2) = Vec(2,0);
}

void TPZElasticity3D::ComputeStrainVector(TPZFMatrix &Strain, TPZFMatrix &DSol) const{
	Strain.Redim(6,1);
	Strain(0,0) = DSol(0,0);
	Strain(1,0) = DSol(1,1);
	Strain(2,0) = DSol(2,2);
	Strain(3,0) = 0.5 * ( DSol(1,0) + DSol(0,1) );
	Strain(4,0) = 0.5 * ( DSol(2,0) + DSol(0,2) );
	Strain(5,0) = 0.5 * ( DSol(2,1) + DSol(1,2) );
}

void TPZElasticity3D::ComputeStressVector(TPZFMatrix &Stress, TPZFMatrix &DSol) const{
	REAL const1 = -1. + this->fPoisson;
	REAL const2 = -1. + this->fPoisson + 2. * fPoisson * fPoisson;
	const REAL E = this->fE;
	const REAL ni = this->fPoisson;
	Stress.Redim(6,1);
	Stress(0,0) = E * ( DSol(0,0) * const1 - ( DSol(1,1) + DSol(2,2) ) * ni ) / const2;
	Stress(1,0) = E * ( DSol(1,1) * const1 - ( DSol(0,0) + DSol(2,2) ) * ni ) / const2;
	Stress(2,0) = E * ( DSol(2,2) * const1 - ( DSol(0,0) + DSol(1,1) ) * ni ) / const2;
	
	REAL const3 = 2. * ( 1. + ni );
	Stress(3,0) = E * ( DSol(1,0) + DSol(0,1) ) / const3;
	Stress(4,0) = E * ( DSol(2,0) + DSol(0,2) ) / const3;
	Stress(5,0) = E * ( DSol(2,1) + DSol(1,2) ) / const3;
}

void TPZElasticity3D::ApplyDirection(TPZFMatrix &StrVec, TPZVec<REAL> &Out) const{
	Out.Resize(3);
	const TPZVec<REAL> &Dir = this->fPostProcessDirection;
	Out[0] = Dir[0] * StrVec(0,0) + Dir[1] * StrVec(3,0) + Dir[2] * StrVec(4,0);
	Out[1] = Dir[1] * StrVec(1,0) + Dir[0] * StrVec(3,0) + Dir[2] * StrVec(5,0);
	Out[2] = Dir[2] * StrVec(2,0) + Dir[0] * StrVec(4,0) + Dir[1] * StrVec(5,0);
}

void TPZElasticity3D::PrincipalDirection(TPZFMatrix &DSol, TPZVec< REAL > &Solout, int direction) const{
	
	TPZFNMatrix<9> StrainTensor(3,3);
	TPZManVector<REAL,3> Eigenvalues;
	TPZFNMatrix<9> Eigenvectors(3,3);
	
	this->ComputeStrainTensor(StrainTensor, DSol);
	int numiterations = 1000;
	REAL tol = TPZElasticity3D::gTolerance;
	bool result = StrainTensor.SolveEigensystemJacobi(numiterations, tol, Solout, Eigenvectors); //Solout is used to store Eigenvaleus, but its values will be replaced below
#ifdef DEBUG        
	if (result == false){
		PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << tol << std::endl;
	}    
#endif
	Solout.Resize(3);
	for(int i = 0; i < 3; i++) Solout[i] = Eigenvectors(direction,i);
	
}

/**
 Save the element data to a stream
 */
void TPZElasticity3D::Write(TPZStream &buf, int withclassid)
{
	TPZMaterial::Write(buf,withclassid);
	buf.Write(&fE,1);
	buf.Write(&fForce[0],3);
	buf.Write(&fFy,1);
	buf.Write(&fPoisson,1);
	buf.Write(&fPostProcessDirection[0],3);
	
}

/**
 Read the element data from a stream
 */
void TPZElasticity3D::Read(TPZStream &buf, void *context)
{
	TPZMaterial::Read(buf,context);
	buf.Read(&fE,1);
	fForce.Resize(3);
	buf.Read(&fForce[0],3);
	buf.Read(&fFy,1);
	buf.Read(&fPoisson,1);
	fPostProcessDirection.Resize(3);
	buf.Read(&fPostProcessDirection[0],3);
}

int TPZElasticity3D::ClassId() const
{
	return TPZELASTICITY3DMATERIALID;
}

void TPZElasticity3D::FillDataRequirements(TPZMaterialData &data){
  	
	TPZMaterial::FillDataRequirements(data);
	data.fNeedsSol = false;
	
}

#ifndef BORLAND
template class TPZRestoreClass<TPZElasticity3D,TPZELASTICITY3DMATERIALID>;
#endif
