/**
 * @file
 * @brief Contains the implementation of the TPZFlowCompMesh methods.
 */

#include "pzflowcmesh.h"
#include "ConsLaw/TPZConsLaw.h"
#include "TPZMaterialDataT.h"
#include "TPZCompElDisc.h"
#include "pzintel.h"


TPZFlowCompMesh::TPZFlowCompMesh(TPZGeoMesh* gr) : TPZRegisterClassId(&TPZFlowCompMesh::ClassId),
TPZCompMesh(gr) {
	
}

TPZFlowCompMesh::TPZFlowCompMesh() : TPZRegisterClassId(&TPZFlowCompMesh::ClassId),
TPZCompMesh() {
	
}

REAL TPZFlowCompMesh::MaxVelocityOfMesh(){
	
	int nel = ElementVec().NElements(), i, nstate, dim, elDim;
	TPZManVector<STATE> density(1), velocity(1);// sol(nstate);
	//TPZManVector<STATE> density(1), sol, velocity(1);// sol(nstate);
	REAL maxvel = 0.0, veloc, sound, gamma;
	STATE press;
	TPZVec<REAL> param(3,0.);
    TPZFNMatrix<9,REAL> axes;
	// loop over all elements, computing the velocity for
	// the non-interface elements.	
	i = 0;
	for(i=0;i<nel;i++){
		
		TPZInterpolationSpace *pElComp =
			dynamic_cast<TPZInterpolationSpace *>(ElementVec()[i]);
		if(!pElComp)
			continue;
		
		TPZMaterial * mat = pElComp->Material();
		if(!mat)PZError << "TPZFlowCompMesh::MaxVelocityOfMesh ERROR: null material.\n";
		
		TPZGeoEl *pElGeo = pElComp->Reference();
		if(!pElGeo)PZError << "TPZFlowCompMesh::MaxVelocityOfMesh ERROR: null element.\n";
		
		dim = mat->Dimension();
		elDim = pElGeo->Dimension();

		TPZMaterialDataT<STATE> elCompData;
		pElComp->InitMaterialData(elCompData);
		if(elDim == dim /*&& pElDisc*/){ // if the dimension of the material fits the
			// dimension of the element and if the element is discontinuous.
			
			TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(mat);
			if(!law)PZError << "TPZFlowCompMesh::MaxVelocityOfMesh2 ERROR: non-fluid material.\n";
			// number of state variables for this material.
			nstate = law->NStateVariables();
			
			// getting the center point of the internal side.
			pElGeo->CenterPoint(pElGeo->NSides()-1,param);//com->Solution(sol,j+100,sol2);
			// verifying the density answer
			pElComp->Solution(param,1,density);
			if(density[0] < 0.0)PZError << "TPZFlowCompMesh::MaxVelocityOfMesh: Negative density\n";
			
			// Getting the velocity answers
			pElComp->Solution(param,6,velocity);
			constexpr bool hasPhi{false};
			pElComp->ComputeSolution(param, elCompData, hasPhi);
			TPZSolVec<STATE> &sol = elCompData.sol;
			press = law->Pressure(sol[0]);
			
#ifndef STATE_COMPLEX
			if(press < 0.0)PZError << "TPZFlowCompMesh::MaxVelocityOfMesh Negative pressure\n";
#endif			
			// retrieving the constant of gas.
			gamma = law->Gamma();
			
			sound = sqrt(gamma*press/density[0]);
			
			// maximal eigenvalue velocity
			veloc = velocity[0] + sound;
			if(veloc > maxvel) maxvel = veloc;
		}else{
			
		}
	}
	return maxvel;
	
}

void TPZFlowCompMesh::CollectFluidMaterials()
{
	std::map<int, TPZMaterial * >::iterator matit;
	
	// buffering the fluid materials
	for(matit = fMaterialVec.begin(); matit != fMaterialVec.end(); matit++)
	{
		TPZConservationLaw * pConsLaw;
		pConsLaw = dynamic_cast<TPZConservationLaw *>(matit->second);
		if(pConsLaw)
		{
			int index =  pConsLaw->Id();
			fFluidMaterial[index] = matit->second;
		}
	}
}

REAL TPZFlowCompMesh::ComputeTimeStep()
{
    REAL maxVel = MaxVelocityOfMesh();
    REAL deltax = LesserEdgeOfMesh();
	
    int nel = fElementVec.NElements();
    int iel;
    REAL meanTimeStep = 0.;
    int numcontr = 0;
	
    for(iel=0; iel<nel; iel++)
    {
		TPZCompEl *cel = fElementVec[iel];
		if(!cel)continue;
		TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(cel);
		int porder = 0;
		if(disc) 
		{
			porder = disc->Degree();
		}
		TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
		if(intel)
		{
			porder = intel->PreferredSideOrder(intel->Reference()->NSides()-1);
		}
		TPZConservationLaw *mat = dynamic_cast<TPZConservationLaw *>( cel->Material());
		if(!mat) continue;
		meanTimeStep += mat->SetTimeStep(maxVel, deltax, porder);
		numcontr++;
    }
	
    return meanTimeStep / (double)numcontr;
}

/**
 * @ingroup material
 * @brief Function for dynamic cast of the material based on map A (second data)
 */
#define FL(A) dynamic_cast<TPZConservationLaw *>(A->second)

/**
 * @ingroup material
 * @brief Maxime value to CFL coefficient
 */
#define MAXCFL 1.e6

void TPZFlowCompMesh::SetCFL(REAL CFL)
{
    if(CFL > MAXCFL) CFL = MAXCFL;
    std::map<int, TPZMaterial * >::iterator matit;
    for(matit=fFluidMaterial.begin(); matit!=fFluidMaterial.end(); matit++)
    {
		FL(matit)->SetCFL(CFL);
    }
}

void TPZFlowCompMesh::ScaleCFL(REAL scale)
{
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit=fFluidMaterial.begin(); matit!=fFluidMaterial.end(); matit++)
	{
		REAL newCFL = FL(matit)->CFL()*scale;
		if(newCFL > MAXCFL) newCFL = MAXCFL;
		std::cout << "CFL = " << newCFL << std::endl;
		FL(matit)->SetCFL(newCFL);
	}
}


void TPZFlowCompMesh::SetContributionTime(TPZContributeTime time)
{
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit=fFluidMaterial.begin(); matit!=fFluidMaterial.end(); matit++)
	{
		FL(matit)->SetContributionTime(time);
	}
}

void TPZFlowCompMesh::SetFlowForcingFunction(ForcingFunctionType<STATE> fp, const int pOrder)
{
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit=fFluidMaterial.begin(); matit!=fFluidMaterial.end(); matit++)
	{
		FL(matit)->SetForcingFunction(fp,pOrder);
	}
}

void TPZFlowCompMesh::AutoBuild()
{
	TPZCompMesh::AutoBuild();
	CollectFluidMaterials();
}


int TPZFlowCompMesh::NFlowMaterials()
{
	return fFluidMaterial.size();
}

/** Returns the first flow material in the mesh */
TPZMaterial * TPZFlowCompMesh::GetFlowMaterial()
{
	TPZMaterial * result = NULL;
	if(fFluidMaterial.size())
	{
		return fFluidMaterial.begin()->second;
	} else
	{
		return result;
	}
}

void TPZFlowCompMesh::SetResidualType(TPZResidualType type)
{
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit=fFluidMaterial.begin(); matit!=fFluidMaterial.end(); matit++)
	{
		FL(matit)->SetResidualType(type);
	}
}

int TPZFlowCompMesh::ClassId() const{
    return Hash("TPZFlowCompMesh") ^ TPZCompMesh::ClassId() << 1;
}

#ifndef BORLAND
template class
TPZRestoreClass< TPZFlowCompMesh>;
#endif

void TPZFlowCompMesh::Write(TPZStream &buf, int withclassid) const
{
	TPZCompMesh::Write(buf,0);
}

void TPZFlowCompMesh::Read(TPZStream &buf, void *context)
{
	TPZCompMesh::Read(buf, context);
	CollectFluidMaterials();
}

void TPZFlowCompMesh::ExpandSolution2()
{
	//TODOCOMPLEX
	ExpandSolution2Internal<STATE>(fSolution);
}

template<class TVar>
void TPZFlowCompMesh::ExpandSolution2Internal(TPZFMatrix<TVar>& mysol)
{
	int nFluid = fFluidMaterial.size();
	
	if(nFluid != 1)
	{
		PZError << "\nTPZFlowCompMesh::ExpandSolution - Wrong number of fluid materials\n";
	}
	
	int nstate = fFluidMaterial.begin()->second->NStateVariables();
	
	fBlock.Resequence();
	int ibl,nblocks = fBlock.NBlocks();
	
	int ic, cols = fSolution.Cols();
	for(ic=0; ic<cols; ic++) {
		for(ibl = 0;ibl<nblocks;ibl++) {
			int size = fSolutionBlock.Size(ibl);
			int position = fSolutionBlock.Position(ibl);
			int lastStatePos = position + size - nstate;
			int ieq;
			//REAL temp;
			STATE temp;
			for(ieq=0; ieq<nstate; ieq++) {
				temp = mysol(position+ieq,ic);
				mysol(position+ieq,ic) = mysol(lastStatePos+ieq,ic);
				mysol(lastStatePos+ieq,ic) = temp;
			}
		}
	}
	
	TPZCompMesh::ExpandSolutionInternal(mysol);
	
	for(ic=0; ic<cols; ic++) {
		for(ibl = 0;ibl<nblocks;ibl++) {
			int size = fSolutionBlock.Size(ibl);
			int position = fSolutionBlock.Position(ibl);
			int lastStatePos = position + size - nstate;
			int ieq;
			//REAL temp;
			STATE temp;
			for(ieq=0; ieq<nstate; ieq++) {
				temp = mysol(position+ieq,ic);
				mysol(position+ieq,ic) = mysol(lastStatePos+ieq,ic);
				mysol(lastStatePos+ieq,ic) = temp;
			}
		}
	}
}
