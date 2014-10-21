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

//enum SOLUTIONVARS{ENone = -1};
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


