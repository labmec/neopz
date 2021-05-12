/**
 * @file
 * @brief Contains implementations of the TPZViscoelastic methods.
 */

#include "pzviscoelastic.h"

TPZViscoelastic::TPZViscoelastic() : 
TPZRegisterClassId(&TPZViscoelastic::ClassId),
TPZMatWithMem<TPZFMatrix<STATE>, TPZElasticity3D>(), fAlpha(-1), fDeltaT(-1), fLambdaV(-1.), fmuV(-1.)
{
    
}

TPZViscoelastic::TPZViscoelastic(int id) : 
TPZRegisterClassId(&TPZViscoelastic::ClassId),
TPZMatWithMem<TPZFMatrix<STATE>, TPZElasticity3D>(id), fAlpha(-1), fDeltaT(-1), fLambdaV(-1.), fmuV(-1.)
{
	
}

TPZViscoelastic::TPZViscoelastic(int id,STATE ElaE,STATE poissonE, STATE lambdaV, STATE muV, STATE alpha, STATE deltaT, TPZVec <STATE> &force): TPZRegisterClassId(&TPZViscoelastic::ClassId),
TPZMatWithMem<TPZFMatrix<STATE>, TPZElasticity3D>(id), fAlpha(alpha), fDeltaT(deltaT), fLambdaV(lambdaV), fmuV(muV)
																																																																																																																	 
{
	STATE lambdaE = (poissonE * ElaE)/((1+poissonE)*(1-2*poissonE));
	STATE muE = ElaE/(2*(1+poissonE));
	STATE lambdaVE = lambdaE-(fDeltaT*lambdaV)/(1+fAlpha*fDeltaT);
	STATE muVE = muE -(fDeltaT*muV)/(1+fAlpha*fDeltaT); 
	STATE ElaVE = muVE*(3*lambdaVE+2*muVE)/(lambdaVE+muVE);
	STATE PoissonVE = lambdaVE/(2*(lambdaVE+muVE));	
	if (lambdaVE < 0 || muVE < 0)
	{
		PZError << "lambdaVE and muVE must be positive. Check your constants values\n";
		DebugStop();
	}
	TPZMatWithMem<TPZFMatrix<STATE>,TPZElasticity3D>::SetMaterialDataHook(ElaVE, PoissonVE); // Creating an elastic material with the viscoelastic properties
	SetForce(force);
	//SetC(); already set in SetMaterialDataHook from elastic material
}

void TPZViscoelastic::SetMaterialDataHooke(STATE ElaE, STATE poissonE, STATE ElaV, STATE poissonV, STATE alpha, STATE deltaT, TPZVec <STATE> &force)
{
	fAlpha = alpha;
	fDeltaT = deltaT;
	STATE lambdaE = (poissonE * ElaE)/((1+poissonE)*(1-2*poissonE));
	STATE muE = ElaE/(2*(1+poissonE));
	STATE lambdaV = (poissonV * ElaV)/((1+poissonV)*(1-2*poissonV));
	STATE muV = ElaV/(2*(1+poissonV));
	fLambdaV = lambdaV; 
	fmuV = muV;
	STATE lambdaVE = lambdaE-(fDeltaT*lambdaV)/(1+fAlpha*fDeltaT);
	STATE muVE = muE -(fDeltaT*muV)/(1+fAlpha*fDeltaT); 
	STATE ElaVE = muVE*(3*lambdaVE+2*muVE)/(lambdaVE+muVE);
	STATE PoissonVE = lambdaVE/(2*(lambdaVE+muVE));
	if (lambdaVE < 0 || muVE < 0)
	{
		PZError << "lambdaVE and muVE must be positive. Check your constants values\n";
		DebugStop();
	}
	TPZMatWithMem<TPZFMatrix<STATE>,TPZElasticity3D>::SetMaterialDataHook(ElaVE, PoissonVE); // Creating an elastic material with the viscoelastic properties
	SetForce(force);
	//SetC(); already set in SetMaterialDataHook from elastic material
}

void TPZViscoelastic::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{	
	int nstate = this->NStateVariables();
	int in;
	STATE val;
    if(fUpdateMem != 0)
    {
        UpdateQsi(data);
    }
    TPZFMatrix<REAL> dphi = data.dphix;
    TPZFMatrix<REAL> phi = data.phi;
    const int phr = phi.Rows();
    int index = data.intGlobPtIndex;

    TPZFNMatrix<6,STATE>  qsi(6,1);
    int rows = this->MemItem(index).Rows();
    if (rows < 6 || rows > 6)
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
        val *= 1./(1+fAlpha*fDeltaT);
        ef(in*nstate+0,0) += weight * val;
        
        //Second equation: fb and fk
        val = 0.;
        val -= qsi(_XY_,0) * dphi(0,in); // |
        val -= qsi(_YY_,0) * dphi(1,in); // fk
        val -= qsi(_YZ_,0) * dphi(2,in); // |
        val *= 1./(1+fAlpha*fDeltaT);
        ef(in*nstate+1,0) += weight * val;
        
        //third equation: fb and fk
        val = 0.;
        val -= qsi(_XZ_,0) * dphi(0,in); // |
        val -= qsi(_YZ_,0) * dphi(1,in); // fk
        val -= qsi(_ZZ_,0) * dphi(2,in); // |
        val *= 1./(1+fAlpha*fDeltaT);
        ef(in*nstate+2,0) += weight * val;
    }
    TPZMatWithMem<TPZFMatrix<STATE>,TPZElasticity3D>::Contribute(data,weight,ek,ef);
}

void TPZViscoelastic::UpdateQsi(TPZMaterialData &data)
{
	TPZFNMatrix<6,STATE> qsi;
	TPZFNMatrix<6,STATE> qsin1(6,1,0.);
	TPZFNMatrix<9,STATE> DSolXYZ(3,3,0.);
	TPZFNMatrix<6,STATE> Strain(6,1);
	int index = data.intGlobPtIndex;
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	DSolXYZ = data.dsol[0];
	qsi = this->MemItem(index);
	
	Strain.Redim(6,1);
	Strain(_XX_,0) = DSolXYZ(0,0);
	Strain(_YY_,0) = DSolXYZ(1,1);
	Strain(_ZZ_,0) = DSolXYZ(2,2);
	Strain(_XY_,0) = 0.5 * ( DSolXYZ(1,0) + DSolXYZ(0,1) );
	Strain(_XZ_,0) = 0.5 * ( DSolXYZ(2,0) + DSolXYZ(0,2) );
	Strain(_YZ_,0) = 0.5 * ( DSolXYZ(2,1) + DSolXYZ(1,2) );
	
	STATE tr;
	tr = Strain(_XX_,0)+Strain(_YY_,0)+Strain(_ZZ_,0);
	
	qsin1(_XX_,0) = (fDeltaT*(-(tr)*fLambdaV - 2*Strain(_XX_,0)*fmuV) + qsi(_XX_,0))/(1 + fAlpha*fDeltaT);
	qsin1(_YY_,0) = (fDeltaT*(-(tr)*fLambdaV - 2*Strain(_YY_,0)*fmuV) + qsi(_YY_,0))/(1 + fAlpha*fDeltaT);
	qsin1(_ZZ_,0) = (fDeltaT*(-(tr)*fLambdaV - 2*Strain(_ZZ_,0)*fmuV) + qsi(_ZZ_,0))/(1 + fAlpha*fDeltaT);
	qsin1(_XY_,0) = (-2*fDeltaT*Strain(_XY_,0)*fmuV + qsi(_XY_,0))/(1 + fAlpha*fDeltaT);
	qsin1(_XZ_,0) = (-2*fDeltaT*Strain(_XZ_,0)*fmuV + qsi(_XZ_,0))/(1 + fAlpha*fDeltaT);
	qsin1(_YZ_,0) = (-2*fDeltaT*Strain(_YZ_,0)*fmuV + qsi(_YZ_,0))/(1 + fAlpha*fDeltaT);
	
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

void TPZViscoelastic::ComputeStressTensor(TPZFMatrix<STATE> &Stress, TPZMaterialData &data) const
{
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZFNMatrix<9,STATE> Dsol = data.dsol[0];
	TPZMatWithMem<TPZFMatrix<STATE>,TPZElasticity3D>::ComputeStressTensor(Stress,Dsol);
	TPZFNMatrix<6,STATE> qsi;
	const int index = data.intGlobPtIndex;
	qsi = *MemItem(index);
	const REAL denominador = 1 + fAlpha*fDeltaT;
	qsi(_XX_,0)*= 1/(1+denominador);
	qsi(_XY_,0)*= 1/(1+denominador);
	qsi(_XZ_,0)*= 1/(1+denominador);
	qsi(_YY_,0)*= 1/(1+denominador);
	qsi(_YZ_,0)*= 1/(1+denominador);
	qsi(_ZZ_,0)*= 1/(1+denominador);
	
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

void TPZViscoelastic::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	if(var == TPZElasticity3D::EDisplacement){
		int i;
		for(i = 0; i < 3; i++){
			TPZVec<STATE> Sol(data.sol[0]); 
			Solout[i] = Sol[i];
		}//for
		return;
	}//TPZElasticity3D::EDisplacement
	
	if(var == TPZElasticity3D::EDisplacementX){
		//    int i;
		TPZVec<STATE> Sol(data.sol[0]); 
		Solout[0] = Sol[0];
		return;
	}//TPZElasticity3D::EDisplacementX
	
	if(var == TPZElasticity3D::EDisplacementY){
		//    int i;
		TPZVec<STATE> Sol(data.sol[0]); 
		Solout[0] = Sol[1];
		return;
	}//TPZElasticity3D::EDisplacementY  
	
	if(var == TPZElasticity3D::EDisplacementZ){
		//    int i;
		TPZVec<STATE> Sol(data.sol[0]); 
		Solout[0] = Sol[2];
		return;
	}//TPZElasticity3D::EDisplacementZ  

	
	if(var == TPZElasticity3D::EPrincipalStrain){
		TPZFNMatrix<9,STATE> StrainTensor(3,3);
		TPZFMatrix<STATE> DSol(data.dsol[0]);
		TPZMatWithMem<TPZFMatrix<STATE>,TPZElasticity3D>::ComputeStrainTensor(StrainTensor, DSol);
		int64_t numiterations = 1000;
		REAL tol = TPZElasticity3D::gTolerance;
        TPZManVector<STATE,3> eigv(3,0.);
		bool result;
        result = StrainTensor.SolveEigenvaluesJacobi(numiterations, tol, &eigv);
        for (int i=0; i<3; i++) {
            Solout[i] = eigv[i];
        }
#ifdef PZDEBUG    
		if (result == false){
			PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << tol << std::endl;
		}
#endif
	}//TPZElasticity3D::EPrincipalStrain
	
	
	if(var == TPZViscoelastic::EViscoStressX){
		TPZFNMatrix<9,STATE> Stress(3,3);
		this->ComputeStressTensor(Stress, data);
		Solout[0] = Stress(0,0);
		return;
	}
	if(var == TPZViscoelastic::EViscoStressY){
		TPZFMatrix<STATE> Stress(3,3);
		this->ComputeStressTensor(Stress, data);
		Solout[0] = Stress(1,1);
		return;
	}
	if(var == TPZViscoelastic::EViscoStressZ){
		TPZFMatrix<STATE> Stress(3,3);
		this->ComputeStressTensor(Stress, data);
		Solout[0] = Stress(2,2);
		return;
	}
}

void TPZViscoelastic::FillDataRequirements(TPZMaterialData &data){
	
	TPZMaterial::FillDataRequirements(data);
	data.fNeedsSol = true;
}

/** Save the element data to a stream */
void TPZViscoelastic::Write(TPZStream &buf, int withclassid) const
{
	TPZMatWithMem<TPZFMatrix<STATE>,TPZElasticity3D>::Write(buf,withclassid);
	buf.Write(&fAlpha, 1);	
	buf.Write(&fDeltaT, 1);		
	buf.Write(&fLambdaV, 1);	
	buf.Write(&fmuV, 1);	
}

/** Read the element data from a stream */
void TPZViscoelastic::Read(TPZStream &buf, void *context)
{
	TPZMatWithMem<TPZFMatrix<STATE>,TPZElasticity3D>::Read(buf,context);
	buf.Read(&fAlpha, 1);
	buf.Read(&fDeltaT, 1);
	buf.Read(&fLambdaV, 1);
	buf.Read(&fmuV, 1);

}

int TPZViscoelastic::ClassId() const{
    return Hash("TPZViscoelastic") ^ TPZMatWithMem<TPZFMatrix<STATE>, TPZElasticity3D>::ClassId() << 1;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZViscoelastic>;
#endif
