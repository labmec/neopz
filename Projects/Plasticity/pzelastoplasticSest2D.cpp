//
//  pzelastoplasticSest2D.cpp
//  PZ
//
//  Created by Diogo Cecilio on 9/23/14.
//
//

#include "pzelastoplasticSest2D.h"
#include "TPZElasticResponse.h"
#include "TPZLadeKim.h"
#include "TPZSandlerDimaggio.h"
#include "TPZYCDruckerPrager.h"
#include "TPZThermoForceA.h"
#include "poroelastoplasticid.h"
#include "TPBrBiotForce.h"


#define SANDLERDIMAGGIOSTEP1 TPZPlasticStep<TPZYCSandlerDimaggioL, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>
#define SANDLERDIMAGGIOSTEP2 TPZPlasticStep<TPZYCSandlerDimaggioL2, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>

/**
 * Default constructor
 */
template <class T, class TMEM>
TPZMatElastoPlasticSest2D<T,TMEM>::TPZMatElastoPlasticSest2D():TPZMatElastoPlastic2D<T, TMEM>(), fZDeformation(0.), fbiot(0.)
{
}

/** Creates a material object and inserts it in the vector of
 *  material pointers of the mesh. Upon return vectorindex
 *  contains the index of the material object within the
 *  vector
 */
template <class T, class TMEM>
TPZMatElastoPlasticSest2D<T,TMEM>::TPZMatElastoPlasticSest2D(int id ,  int PlaneStrainOrPlaneStress):TPZMatElastoPlastic2D<T, TMEM>(id,PlaneStrainOrPlaneStress), fZDeformation(0.), fbiot(0.)
{
  
}

/** Creates a material object based on the referred object and
 *  inserts it in the vector of material pointers of the mesh.
 *  Upon return vectorindex contains the index of the material
 *  object within the vector
 */
template <class T, class TMEM>
TPZMatElastoPlasticSest2D<T,TMEM>::TPZMatElastoPlasticSest2D(const TPZMatElastoPlasticSest2D<T,TMEM> &mat):TPZMatElastoPlastic2D<T, TMEM>(mat), fZDeformation(mat.fZDeformation), fbiot(mat.fbiot)
{
  TPBrBiotForce * func = dynamic_cast<TPBrBiotForce *>(this->fForcingFunction.operator->());
    
  if (func) {
    TPZAutoPointer<TPZFunction<STATE> > clonefunc = new TPBrBiotForce(*func);
    this->fForcingFunction = clonefunc;
  }
    else
    {
    }
}

template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::ComputeDeltaStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &DeltaStrain)
{
  TPZFNMatrix<9> DSolXYZ(3,3,0.);
  data.axes.Multiply(data.dsol[0],DSolXYZ,1/*transpose*/);
  if (DeltaStrain.Rows() != 6) {
    DebugStop();
  }
  //  DeltaStrain.Redim(3,1);
  DeltaStrain(_XX_,0) = DSolXYZ(0,0);
  DeltaStrain(_YY_,0) = DSolXYZ(1,1);
  DeltaStrain(_XY_,0) = 0.5 * ( DSolXYZ(1,0) + DSolXYZ(0,1) );
  DeltaStrain(_XZ_,0) = 0.;
  DeltaStrain(_YZ_,0) = 0.;
  DeltaStrain(_ZZ_,0) = fZDeformation;
}

template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::ApplyDeltaStrainComputeDep(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep)
{
  
  if (DeltaStrain.Rows() != 6) {
    DebugStop();
  }
  TPZMatElastoPlastic<T,TMEM>::ApplyDeltaStrainComputeDep(data,DeltaStrain,Stress,Dep);//
  if (this->fPlaneStrain) //
  {//
    
  }
  else//PlaneStress
  {
    DebugStop();
  }
}

/**
 * Save the element data to a stream
 */
template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::Write(TPZStream &buf, int withclassid)
{
    TPZMatElastoPlastic2D<T,TMEM>::Write(buf,withclassid);
    buf.Write(&fZDeformation);
    buf.Write(&fbiot);
}

/**
 * Read the element data from a stream
 */
template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::Read(TPZStream &buf, void *context)
{
    TPZMatElastoPlastic2D<T,TMEM>::Read(buf,context);
    buf.Read(&fZDeformation);
    buf.Read(&fbiot);
    
}



template <class T, class TMEM>
int TPZMatElastoPlasticSest2D<T,TMEM>::ClassId() const
{
  DebugStop();
  return TPZSANDLERDIMAGGIOL_ID + NUMPLASTICMODELS+20;
}


template <>
int TPZMatElastoPlasticSest2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>, TPZElastoPlasticMem>::ClassId() const
{
  return TPZSANDLERDIMAGGIOL_ID + NUMPLASTICMODELS+20;
}

template <>
int TPZMatElastoPlasticSest2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZElastoPlasticMem>::ClassId() const
{
  return TPZSANDLERDIMAGGIOL2_ID + NUMPLASTICMODELS+20;
}

#include "pzsandlerextPV.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"

template<>
int TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem>::ClassId() const
{
  return TPZSANDLERDIMAGGIOPV_ID + NUMPLASTICMODELS+20;
}

template<>
int TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem>::ClassId() const
{
  return TPZMOHRCOULOMBPV_ID + NUMPLASTICMODELS+20;
}

template <class T, class TMEM>
int TPZMatElastoPlasticSest2D<T,TMEM>::VariableIndex(const std::string &name)
{
  if(!strcmp("StrainVol",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainVol;
  if(!strcmp("StrainXX",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXX;
  if(!strcmp("StrainYY",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainYY;
  if(!strcmp("StrainZZ",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainZZ;
  if(!strcmp("StrainXY",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXY;
  if(!strcmp("StrainXZ",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXZ;
  if(!strcmp("StrainYZ",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainYZ;
  if(!strcmp("ElStrainVol",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainVol;
  if(!strcmp("ElStrainXX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXX;
  if(!strcmp("ElStrainYY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainYY;
  if(!strcmp("ElStrainZZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainZZ;
  if(!strcmp("ElStrainXY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXY;
  if(!strcmp("ElStrainXZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXZ;
  if(!strcmp("ElStrainYZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainYZ;
  if(!strcmp("PlStrainVol",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainVol;
  if(!strcmp("PlStrainXX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXX;
  if(!strcmp("PlStrainYY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainYY;
  if(!strcmp("PlStrainZZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainZZ;
  if(!strcmp("PlStrainXY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXY;
  if(!strcmp("PlStrainXZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXZ;
  if(!strcmp("PlStrainYZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainYZ;
  if(!strcmp("PlStrainSqJ2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainSqJ2;
  if(!strcmp("PlStrainSqJ2El",		name.c_str()))  return 100;//return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainSqJ2El;
  if(!strcmp("PlAlpha",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlAlpha;
  if(!strcmp("DisplacementX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementX;
  if(!strcmp("DisplacementY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementY;
  if(!strcmp("DisplacementZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementZ;
  if(!strcmp("DisplacementTotal",	name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementTotal;
  if(!strcmp("TotStressI1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressI1;
  if(!strcmp("TotStressJ2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressJ2;
  if(!strcmp("TotStressXX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXX;
  if(!strcmp("TotStressYY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressYY;
  if(!strcmp("TotStressZZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressZZ;
  if(!strcmp("TotStressXY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXY;
  if(!strcmp("TotStressXZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXZ;
  if(!strcmp("TotStressYZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressYZ;
  if(!strcmp("TotStress1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress1;
  if(!strcmp("TotStress2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress2;
  if(!strcmp("TotStress3",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress3;
  if(!strcmp("EffStressI1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressI1;
  if(!strcmp("EffStressJ2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressJ2;
  if(!strcmp("EffStressXX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXX;
  if(!strcmp("EffStressYY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressYY;
  if(!strcmp("EffStressZZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressZZ;
  if(!strcmp("EffStressXY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXY;
  if(!strcmp("EffStressXZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXZ;
  if(!strcmp("EffStressYZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressYZ;
  if(!strcmp("EffStress1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress1;
  if(!strcmp("EffStress2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress2;
  if(!strcmp("EffStress3",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress3;
  if(!strcmp("YieldSurface1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface1;
  if(!strcmp("YieldSurface2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface2;
  if(!strcmp("YieldSurface3",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface3;
  if(!strcmp("POrder",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPOrder;
  if(!strcmp("NSteps",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ENSteps;
  if(!strcmp("PorePressure",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPorePressure;
  if(!strcmp("MatPorosity",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EMatPorosity;
  if(!strcmp("MatE",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EMatE;
  if(!strcmp("MatPoisson",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EMatPoisson;
  
  //return TPZMatWithMem<TMEM>::VariableIndex(name);
  PZError << "TPZMatElastoPlastic::VariableIndex Error\n";
  DebugStop();
  return -1;
}

template <class T, class TMEM>
int TPZMatElastoPlasticSest2D<T,TMEM>::NSolutionVariables(int var)
{
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainVol)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXX	)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainYY)			 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainZZ)			 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXY)			 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXZ)			 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainYZ)			 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainVol)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXX)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainYY)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainZZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXY)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainYZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainVol)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXX)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainYY)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainZZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXY)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainYZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainSqJ2)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainSqJ2El)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlAlpha)			 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementX)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementY)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementTotal)	 return 2;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressI1)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressJ2)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXX)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressYY)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressZZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXY)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressYZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress1)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress2)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress3)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressI1)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressJ2)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXX)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressYY)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressZZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXY)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressYZ)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress1)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress2)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress3)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface1)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface2)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface3)		 return 1; // Should never be called
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPOrder)			 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ENSteps)			 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPorePressure)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EMatPorosity)		 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EMatE)			 return 1;
  if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EMatPoisson)		 return 1;
  
  if(var == 100) return 1;
  return TPZMatWithMem<TMEM>::NSolutionVariables(var);
}

template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
  int intPt = data.intGlobPtIndex;
  TMEM &Memory = TPZMatWithMem<TMEM>::fMemory[intPt];
  T plasticloc(this->fPlasticity);
  plasticloc.SetState(Memory.fPlasticState);
    TPZTensor<STATE> Sigma = Memory.fSigma;
    STATE normdsol = Norm(data.dsol[0]);
    if (normdsol != 0.) {
        TPZTensor<REAL> EpsT;
        TPZFNMatrix<6,STATE> deltastrain(6,1,0.);
        ComputeDeltaStrainVector(data, deltastrain);
        
        EpsT.CopyFrom(deltastrain);
        EpsT.Add(plasticloc.GetState().fEpsT, 1.);
        
        plasticloc.ApplyStrainComputeSigma(EpsT, Sigma);
    }
    TPZPlasticState<STATE> PState = plasticloc.GetState();
    TPZTensor<REAL> totalStrain = PState.fEpsT;
    TPZTensor<REAL> plasticStrain = PState.fEpsP;
    
    

  //Elastic Strain
  TPZTensor<REAL> elasticStrain = totalStrain; // Look at line below
  elasticStrain -= plasticStrain; // here it becomes elasticStrain
  

  //Total Stress
  TPZTensor<REAL> totalStress = Sigma;
  
  TPZManVector<STATE,3> AlphagradP(2,0.);
  TPZFNMatrix<1,STATE> AlphaP(1,1,0.);
  if (this->fForcingFunction) {
    this->fForcingFunction->Execute(data.x,AlphagradP,AlphaP);
    totalStress.XX() -= AlphaP(0,0);
    totalStress.YY() -= AlphaP(0,0);
    totalStress.ZZ() -= AlphaP(0,0);
  }
  
    STATE ux = Memory.fDisplacement[0];
    STATE uy = Memory.fDisplacement[1];
    
  switch (var) {
    // Total Strain
    case EStrainVol:
      Solout[0] = totalStrain.XX() + totalStrain.YY() + totalStrain.ZZ();
      break;
    case EStrainXX:
      Solout[0] = totalStrain.XX();
      break;
    case EStrainYY:
      Solout[0] = totalStrain.YY();
      break;
    case EStrainZZ:
      Solout[0] = totalStrain.ZZ();
      break;
    case EStrainXY:
      Solout[0] = totalStrain.XY();
      break;
    case EStrainXZ:
      Solout[0] = totalStrain.XZ();
      break;
    case EStrainYZ:
      Solout[0] = totalStrain.YZ();
      break;
    // Elastic Strain
    case EElStrainVol:
      Solout[0] = elasticStrain.XX() + elasticStrain.YY() + elasticStrain.ZZ();
      break;
    case EElStrainXX:
      Solout[0] = elasticStrain.XX();
      break;
    case EElStrainYY:
      Solout[0] = elasticStrain.YY();
      break;
    case EElStrainZZ:
      Solout[0] = elasticStrain.ZZ();
      break;
    case EElStrainXY:
      Solout[0] = elasticStrain.XY();
      break;
    case EElStrainXZ:
      Solout[0] = elasticStrain.XZ();
      break;
    case EElStrainYZ:
      Solout[0] = elasticStrain.YZ();
      break;
    // Plastic Strain
    case EPlStrainVol:
      Solout[0] = plasticStrain.XX() + plasticStrain.YY() + plasticStrain.ZZ();
      break;
    case EPlStrainXX:
      Solout[0] = plasticStrain.XX();
      break;
    case EPlStrainYY:
      Solout[0] = plasticStrain.YY();
      break;
    case EPlStrainZZ:
      Solout[0] = plasticStrain.ZZ();
      break;
    case EPlStrainXY:
      Solout[0] = plasticStrain.XY();
      break;
    case EPlStrainXZ:
      Solout[0] = plasticStrain.XZ();
      break;
    case EPlStrainYZ:
      Solout[0] = plasticStrain.YZ();
      break;
    // SqJ2 and alpha
    case EPlStrainSqJ2:
      Solout[0] = sqrt(plasticStrain.J2());
      break;
    case EPlStrainSqJ2El:
      DebugStop();
      break;
    case EPlAlpha:
   		Solout[0] = TPZMatWithMem<TMEM>::fMemory[intPt].fPlasticState.fAlpha;
      break;
    // Displacement
    case EDisplacementX:
      Solout[0] = TPZMatWithMem<TMEM>::fMemory[intPt].fDisplacement[0];
      break;
    case EDisplacementY:
      Solout[0] = TPZMatWithMem<TMEM>::fMemory[intPt].fDisplacement[1];
      break;
    case EDisplacementZ:
      Solout[0] = 0.; //DUVIDA
      break;
    case EDisplacementTotal:
          

          
          Solout[0] = ux;
          Solout[1] = uy;
      break;
      // Total Stress
    case ETotStressI1:
      Solout[0] = totalStress.I1();
      break;
    case ETotStressJ2:
      Solout[0] = totalStress.J2();
      break;
    case ETotStressXX:
      Solout[0] = totalStress.XX();
      break;
    case ETotStressYY:
      Solout[0] = totalStress.YY();
      break;
    case ETotStressZZ:
      Solout[0] = totalStress.ZZ();
      break;
    case ETotStressXY:
      Solout[0] = totalStress.XY();
      break;
    case ETotStressXZ:
      Solout[0] = totalStress.XZ();
      break;
    case ETotStressYZ:
      Solout[0] = totalStress.YZ();
      break;
    case ETotStress1:
      {
        TPZTensor<STATE> eigenval;
        totalStress.EigenValue(eigenval);
        Solout[0] = eigenval.XX();
      }
      break;
    case ETotStress2:
      {
        TPZTensor<STATE> eigenval;
        totalStress.EigenValue(eigenval);
        Solout[0] = eigenval.YY();
      }
      break;
    case ETotStress3:
      {
        TPZTensor<STATE> eigenval;
        totalStress.EigenValue(eigenval);
        Solout[0] = eigenval.ZZ();
      }
      break;
    // Effective stress
    case EEffStressI1:
      Solout[0] = Sigma.I1();
      break;
    case EEffStressJ2:
      Solout[0] = Sigma.J2();
      break;
    case EEffStressXX:
      Solout[0] = Sigma.XX();
      break;
    case EEffStressYY:
      Solout[0] = Sigma.YY();
      break;
    case EEffStressZZ:
      Solout[0] = Sigma.ZZ();
      break;
    case EEffStressXY:
      Solout[0] = Sigma.XY();
      break;
    case EEffStressXZ:
      Solout[0] = Sigma.XZ();
      break;
    case EEffStressYZ:
      Solout[0] = Sigma.YZ();
      break;
    case EEffStress1:
      {
        TPZTensor<STATE> eigenval;
        Sigma.EigenValue(eigenval);
        Solout[0] = eigenval.XX();
      }
      break;
    case EEffStress2:
      {
        TPZTensor<STATE> eigenval;
        Sigma.EigenValue(eigenval);
        Solout[0] = eigenval.YY();
      }
      break;
    case EEffStress3:
      {
        TPZTensor<STATE> eigenval;
        Sigma.EigenValue(eigenval);
        Solout[0] = eigenval.ZZ();
      }
      break;
    // Yield Surface
    case EYieldSurface1:
      {
        TPZManVector<STATE,3> yieldVal(3,0.);
        plasticloc.Phi(elasticStrain,yieldVal);
        Solout[0] = yieldVal[0];
      }
      break;
    case EYieldSurface2:
      {
        TPZManVector<STATE,3> yieldVal(3,0.);
        plasticloc.Phi(elasticStrain,yieldVal);
        Solout[0] = yieldVal[1];
      }
      break;
    case EYieldSurface3:
      Solout[0] = 0.;
      break;
    // Simulation
    case EPOrder:
      Solout[0] = data.p;
      break;
    case ENSteps:
   		Solout[0] = TPZMatWithMem<TMEM>::fMemory[intPt].fPlasticSteps;
      break;
      // Pore pressure
    case EPorePressure:
      Solout[0] = AlphaP(0,0)/fbiot;
      break;
      // Material
    case EMatPorosity:
      Solout[0] = 6378; // AQUINATHAN raio da terra
      break;
    case EMatE:
      Solout[0] = plasticloc.fER.E();
      break;
    case EMatPoisson:
      Solout[0] = plasticloc.fER.Poisson();
      break;
    default:
      DebugStop();
      break;
  }
}


template class TPZMatElastoPlasticSest2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>, TPZPoroElastoPlasticMem>;
template class TPZMatElastoPlasticSest2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZPoroElastoPlasticMem>;
//template class TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZPoroElastoPlasticMem>;
template class TPZMatElastoPlasticSest2D<TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> , TPZPoroElastoPlasticMem>;

template class TPZRestoreClass< TPZMatElastoPlasticSest2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>, TPZElastoPlasticMem>, TPZSANDLERDIMAGGIOL_ID + NUMPLASTICMODELS +20>;
template class TPZRestoreClass< TPZMatElastoPlasticSest2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>, TPZElastoPlasticMem>,TPZSANDLERDIMAGGIOL2_ID + NUMPLASTICMODELS+20>;


template class TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem>;
template class TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem>;


template class TPZRestoreClass< TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem>, TPZMOHRCOULOMBPV_ID + NUMPLASTICMODELS+20>;
template class TPZRestoreClass< TPZMatElastoPlasticSest2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse>, TPZElastoPlasticMem>,TPZSANDLERDIMAGGIOPV_ID + NUMPLASTICMODELS+20>;


