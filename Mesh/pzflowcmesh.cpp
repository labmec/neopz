//$Id: pzflowcmesh.cpp,v 1.11 2004-05-25 12:58:55 erick Exp $

#include "pzflowcmesh.h"
#include "TPZCompElDisc.h"
#include "TPZConservationLaw.h"
//#include "TPZInterfaceEl.h"

TPZFlowCompMesh::TPZFlowCompMesh(TPZGeoMesh* gr) : TPZCompMesh(gr) {

}

TPZFlowCompMesh::TPZFlowCompMesh() : TPZCompMesh() {

}

REAL TPZFlowCompMesh::MaxVelocityOfMesh(){

  int nel = ElementVec().NElements(), i, nstate, dim, elDim;
  TPZManVector<REAL> density(1), sol, velocity(1);// sol(nstate);
  REAL maxvel = 0.0, veloc, sound, press, gamma;
  TPZVec<REAL> param(3,0.);

  // loop over all elements, computing the velocity for
  // the non-interface elements.

  i = 0;
  while(i<nel){

    TPZCompEl *pElComp = ElementVec()[i];
    if(!pElComp)
    {
       cout << "TPZFlowCompMesh::MaxVelocityOfMesh: No associated computation element.\n";
       continue;
    }

    TPZMaterial *mat = pElComp->Material();
    if(!mat)PZError << "TPZFlowCompMesh::MaxVelocityOfMesh ERROR: null material.\n";

    TPZGeoEl *pElGeo = pElComp->Reference();
    if(!pElGeo)PZError << "TPZFlowCompMesh::MaxVelocityOfMesh ERROR: null element.\n";

    dim = mat->Dimension();
    elDim = pElGeo->Dimension();

    TPZCompElDisc *pElDisc = dynamic_cast<TPZCompElDisc *>(pElComp);

    if(elDim == dim && pElDisc){ // if the dimension of the material fits the
       // dimension of the element and if the element is discontinuous.

#ifdef _AUTODIFF
       TPZConservationLaw2 *law = dynamic_cast<TPZConservationLaw2 *>(mat);
       if(!law)PZError << "TPZFlowCompMesh::MaxVelocityOfMesh2 ERROR: non-fluid material.\n";
#else
       TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(mat);
       if(!law)PZError << "TPZFlowCompMesh::MaxVelocityOfMesh2 ERROR: non-fluid material.\n";
#endif
       // number of state variables for this material.
       nstate = law->NStateVariables();

       // getting the center point of the internal side.
       pElGeo->CenterPoint(pElGeo->NSides()-1,param);//com->Solution(sol,j+100,sol2);
       // verifying the density answer
       pElComp->Solution(param,1,density);
       if(density[0] < 0.0)PZError << "TPZFlowCompMesh::MaxVelocityOfMesh: Negative density\n";

       // Getting the velocity answers
       pElComp->Solution(param,6,velocity);
       // getting the whole vector of solutions
       sol.Resize(nstate);
       pElComp->Solution(param,5,sol);

       press = law->Pressure(sol);
       if(press < 0.0)PZError << "TPZFlowCompMesh::MaxVelocityOfMesh Negative pressure\n";

       // retrieving the constant of gas.
       gamma = law->Gamma();

       sound = sqrt(gamma*press/density[0]);

       // maximal eigenvalue velocity
       veloc = velocity[0] + sound;
       if(veloc > maxvel) maxvel = veloc;
    }else{

    }
    i++;
  }
  return maxvel;

}

void TPZFlowCompMesh::CollectFluidMaterials()
{
   int i, NumMat;

   NumMat = fMaterialVec.NElements();

   // buffering the fluid materials
   for(i = 0; i < NumMat; i++)
   {
      TPZConservationLaw2 * pConsLaw;
      pConsLaw = dynamic_cast<TPZConservationLaw2 *>(fMaterialVec[i]);
      if(pConsLaw)
      {
         int index =  fFluidMaterial.AllocateNewElement();
	 fFluidMaterial[index] = pConsLaw;
      }
   }
}


void TPZFlowCompMesh::ComputeTimeStep()
{
    REAL maxVel = MaxVelocityOfMesh();
    REAL deltax = LesserEdgeOfMesh();

    TPZCompElDisc *disc;
    int degree = disc->gDegree;

    int i, NumFluid = fFluidMaterial.NElements();
    for(i = 0; i < NumFluid; i++)
    {
       fFluidMaterial[i]->SetTimeStep(maxVel, deltax, degree);
    }
}

void TPZFlowCompMesh::SetCFL(REAL CFL)
{
    int i, NumFluid = fFluidMaterial.NElements();
    for(i = 0; i < NumFluid; i++)
    {
       fFluidMaterial[i]->SetCFL(CFL);
    }
}

void TPZFlowCompMesh::ScaleCFL(REAL scale)
{
    int i, NumFluid = fFluidMaterial.NElements();
    for(i = 0; i < NumFluid; i++)
    {
       fFluidMaterial[i]->SetCFL(fFluidMaterial[i]->CFL()*scale);
    }
}


void TPZFlowCompMesh::SetContributionTime(TPZContributeTime time)
{
   int i, NumFluid;
   NumFluid = fFluidMaterial.NElements();
   for(i = 0; i < NumFluid; i++)
   {
      fFluidMaterial[i]->SetContributionTime(time);
   }
}

void TPZFlowCompMesh::SetFlowforcingFunction(void (*fp)(TPZVec<REAL> &loc,
					 TPZVec<REAL> &result))
{
   int i, NumFluid;
   NumFluid = fFluidMaterial.NElements();
   for(i = 0; i < NumFluid; i++)
   {
      fFluidMaterial[i]->SetForcingFunction(fp);
   }
}

void TPZFlowCompMesh::AutoBuild()
{
   TPZCompMesh::AutoBuild();
   CollectFluidMaterials();
}


TPZMaterial * TPZFlowCompMesh::GetFlowMaterial(int i)
{
   return fFluidMaterial[i];
}

int TPZFlowCompMesh::NFlowMaterials()
{
   return fFluidMaterial.NElements();
}

void TPZFlowCompMesh::SetResidualType(TPZResidualType type)
{
   int i, NumFluid;
   NumFluid = fFluidMaterial.NElements();
   for(i = 0; i < NumFluid; i++)
   {
      fFluidMaterial[i]->SetResidualType(type);
   }
}

int TPZFlowCompMesh::ClassId() const 
{
   return TPZFLOWCOMPMESHID;
}

void TPZFlowCompMesh::Write(TPZStream &buf, int withclassid)
{
   TPZSaveable::Write(buf,withclassid);
   TPZCompMesh::Write(buf,0);
}

void TPZFlowCompMesh::Read(TPZStream &buf, void *context)
{
   TPZSaveable::Read(buf, context);
   TPZCompMesh::Read(buf, context);
   CollectFluidMaterials();
}


void TPZFlowCompMesh::ExpandSolution2()
{
   int nFluid = fFluidMaterial.NElements();

   if(nFluid != 1)
   {
      PZError << "\nTPZFlowCompMesh::ExpandSolution - Wrong number of fluid materials\n";
   }

   int nstate = fFluidMaterial[0]->NStateVariables();

   fBlock.Resequence();
   int ibl,nblocks = fBlock.NBlocks();

   int ic, cols = fSolution.Cols();
   for(ic=0; ic<cols; ic++) {
      for(ibl = 0;ibl<nblocks;ibl++) {
         int size = fSolutionBlock.Size(ibl);
         int position = fSolutionBlock.Position(ibl);
	 int lastStatePos = position + size - nstate;
         int ieq;
	 double temp;
         for(ieq=0; ieq<nstate; ieq++) {
            temp = fSolution(position+ieq,ic);
            fSolution(position+ieq,ic) = fSolution(lastStatePos+ieq,ic);
	    fSolution(lastStatePos+ieq,ic) = temp;
         }
      }
   }

   TPZCompMesh::ExpandSolution();

   for(ic=0; ic<cols; ic++) {
      for(ibl = 0;ibl<nblocks;ibl++) {
         int size = fSolutionBlock.Size(ibl);
         int position = fSolutionBlock.Position(ibl);
	 int lastStatePos = position + size - nstate;
         int ieq;
	 double temp;
         for(ieq=0; ieq<nstate; ieq++) {
            temp = fSolution(position+ieq,ic);
            fSolution(position+ieq,ic) = fSolution(lastStatePos+ieq,ic);
	    fSolution(lastStatePos+ieq,ic) = temp;
         }
      }
   }
}
