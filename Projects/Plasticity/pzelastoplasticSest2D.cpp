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


#define SANDLERDIMAGGIOSTEP1 TPZPlasticStep<TPZYCSandlerDimaggioL, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>
#define SANDLERDIMAGGIOSTEP2 TPZPlasticStep<TPZYCSandlerDimaggioL2, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>

/**
 * Default constructor
 */
template <class T, class TMEM>
TPZMatElastoPlasticSest2D<T,TMEM>::TPZMatElastoPlasticSest2D():TPZMatElastoPlastic2D<T, TMEM>(), fZDeformation(0.)
{
}

/** Creates a material object and inserts it in the vector of
 *  material pointers of the mesh. Upon return vectorindex
 *  contains the index of the material object within the
 *  vector
 */
template <class T, class TMEM>
TPZMatElastoPlasticSest2D<T,TMEM>::TPZMatElastoPlasticSest2D(int id ,  int PlaneStrainOrPlaneStress, STATE sigmaZ):TPZMatElastoPlastic2D<T, TMEM>(id,PlaneStrainOrPlaneStress, sigmaZ), fZDeformation(0.)
{
}

/** Creates a material object based on the referred object and
 *  inserts it in the vector of material pointers of the mesh.
 *  Upon return vectorindex contains the index of the material
 *  object within the vector
 */
template <class T, class TMEM>
TPZMatElastoPlasticSest2D<T,TMEM>::TPZMatElastoPlasticSest2D(const TPZMatElastoPlasticSest2D<T,TMEM> &mat):TPZMatElastoPlastic2D<T, TMEM>(mat), fZDeformation(mat.fZDeformation)
{
    
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
	if(!strcmp("EStrainVol",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainVol;
	if(!strcmp("EStrainXX",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXX;
	if(!strcmp("EStrainYY",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainYY;
	if(!strcmp("EStrainZZ",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainZZ;
	if(!strcmp("EStrainXY",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXY;
	if(!strcmp("EStrainXZ",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXZ;
	if(!strcmp("EStrainYZ",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EStrainYZ;
	if(!strcmp("EElStrainVol",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainVol;
	if(!strcmp("EElStrainXX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXX;
	if(!strcmp("EElStrainYY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainYY;
	if(!strcmp("EElStrainZZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainZZ;
	if(!strcmp("EElStrainXY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXY;
	if(!strcmp("EElStrainXZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXZ;
	if(!strcmp("EElStrainYZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainYZ;
	if(!strcmp("EPlStrainVol",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainVol;
	if(!strcmp("EPlStrainXX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXX;
	if(!strcmp("EPlStrainYY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainYY;
	if(!strcmp("EPlStrainZZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainZZ;
	if(!strcmp("EPlStrainXY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXY;
	if(!strcmp("EPlStrainXZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXZ;
	if(!strcmp("EPlStrainYZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainYZ;
	if(!strcmp("EPlStrainSqJ2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainSqJ2;
	if(!strcmp("EPlStrainSqJ2El",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainSqJ2El;
	if(!strcmp("EPlAlpha",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPlAlpha;
	if(!strcmp("EDisplacementX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementX;
	if(!strcmp("EDisplacementY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementY;
	if(!strcmp("EDisplacementZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementZ;
	if(!strcmp("EDisplacementTotal",	name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementTotal;
	if(!strcmp("ETotStressI1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressI1;
	if(!strcmp("ETotStressJ2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressJ2;
	if(!strcmp("ETotStressXX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXX;
	if(!strcmp("ETotStressYY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressYY;
	if(!strcmp("ETotStressZZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressZZ;
	if(!strcmp("ETotStressXY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXY;
	if(!strcmp("ETotStressXZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXZ;
	if(!strcmp("ETotStressYZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressYZ;
	if(!strcmp("ETotStress1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress1;
	if(!strcmp("ETotStress2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress2;
	if(!strcmp("ETotStress3",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress3;
	if(!strcmp("EEffStressI1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressI1;
	if(!strcmp("EEffStressJ2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressJ2;
	if(!strcmp("EEffStressXX",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXX;
	if(!strcmp("EEffStressYY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressYY;
	if(!strcmp("EEffStressZZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressZZ;
	if(!strcmp("EEffStressXY",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXY;
	if(!strcmp("EEffStressXZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXZ;
	if(!strcmp("EEffStressYZ",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressYZ;
	if(!strcmp("EEffStress1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress1;
	if(!strcmp("EEffStress2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress2;
	if(!strcmp("EEffStress3",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress3;
	if(!strcmp("EYieldSurface1",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface1;
	if(!strcmp("EYieldSurface2",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface2;
	if(!strcmp("EYieldSurface3",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface3;
	if(!strcmp("EPOrder",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPOrder;
	if(!strcmp("ENSteps",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::ENSteps;
	if(!strcmp("EPorePressure",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EPorePressure;
	if(!strcmp("EMatPorosity",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EMatPorosity;
	if(!strcmp("EMatE",			name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EMatE;
	if(!strcmp("EMatPoisson",		name.c_str()))  return TPZMatElastoPlasticSest2D<T,TMEM>::EMatPoisson;

	//return TPZMatWithMem<TMEM>::VariableIndex(name);
	PZError << "TPZMatElastoPlastic::VariableIndex Error\n";
	return -1;
}

template <class T, class TMEM>
int TPZMatElastoPlasticSest2D<T,TMEM>::NSolutionVariables(int var)
{
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainVol)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXX	)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainYY)			 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainZZ)			 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXY)			 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainXZ)			 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EStrainYZ)			 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainVol)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXX)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainYY)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainZZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXY)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainXZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EElStrainYZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainVol)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXX)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainYY)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainZZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXY)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainXZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainYZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainSqJ2)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlStrainSqJ2El)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPlAlpha)			 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementX)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementY)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EDisplacementTotal)	 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressI1)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressJ2)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXX)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressYY)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressZZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXY)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressXZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStressYZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress1)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress2)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ETotStress3)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressI1)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressJ2)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXX)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressYY)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressZZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXY)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressXZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStressYZ)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress1)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress2)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EEffStress3)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface1)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface2)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EYieldSurface3)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPOrder)			 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::ENSteps)			 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EPorePressure)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EMatPorosity)		 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EMatE)			 return -1;
	if(var == TPZMatElastoPlasticSest2D<T,TMEM>::EMatPoisson)		 return -1;

	if(var == 100) return 1;
	return TPZMatWithMem<TMEM>::NSolutionVariables(var);
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


