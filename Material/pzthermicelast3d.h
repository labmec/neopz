// -*- c++ -*-

//$Id: pzthermicelast3d.h,v 1.3 2009-09-01 19:44:48 phil Exp $

#ifndef PZTHERMICELAST3D
#define PZTHERMICELAST3D

#include "pzelast3d.h"

/** This class implements a 3D isotropic elasticity material with thermal stress.
 * The thermal field can be a constant value or a field obtained from a previous simmulation.
 * In case of using a field previously obtained it is input by using TPZReferredCompEl and the field
 * comes in the TPZVec<REAL> sol parameter of Contribute method.
 *  @since May 10, 2006.
 */
class TPZThermicElast3D : public TPZElasticity3D {

private:

/** Flag indicating a constant field of temperature (fReferred = false) or a field obtained in a previous simmulation (fReferred = true)
 * Default value is true.
 */
bool fReferred;

/** When solving problem with a constant field of temperature (fReferred = false) fFinalTemperature is the final temperature
 * of the body. Thermal stress is given by the difference of fFinalTemperature and fRefTemperature.
 * When fReferred = true fFinalTemperature does not make sense and is not used.
 */
REAL fFinalTemperature;

public:

/** Thermal expansion coefficient of the material
 */
REAL fThermalCoeff;

/** Temperature of reference. The thermal stress is given by the difference of the temperature field to this reference value.
 */
REAL fRefTemperature;

/** Class constructor.
 * @param nummat - material ID.
 * @param ThermalCoeff - thermal expansion coefficient of the material
 * @param RefTemp - temperature of reference
 * @param E - Young's modulus.
 * @param poisson - poisson's ratio 
 * @param force - external forces
 */ 
TPZThermicElast3D(int nummat, REAL ThermalCoeff, REAL RefTemp, REAL E, REAL poisson, TPZVec<REAL> &force);

/** Class destructor.
 */
virtual ~TPZThermicElast3D();

void SetConstantTemperatureField(REAL FinalTemp){
  this->fFinalTemperature = FinalTemp;
  this->fReferred = false;
}//method

void SetReferredTemperatureField(){
  this->fFinalTemperature = 0.; //any value, it will not be used
  this->fReferred = true;
}//method

bool IsReferredTemperatureField(){
  return this->fReferred;
}//method

void ContributeThermalStress(TPZVec<REAL> &sol, TPZFMatrix &phi, TPZFMatrix &dphi, REAL weight, TPZFMatrix &ef);

/** Contribute to stiff matrix and load vector.
 *  See base class to more informations.
 */
virtual void Contribute(TPZMaterialData &data,
						  REAL weight,
						  TPZFMatrix &ek,
						  TPZFMatrix &ef);

virtual void Contribute(TPZMaterialData &data,
						  REAL weight,
						  TPZFMatrix &ef)
{
	 TPZElasticity3D::Contribute(data,weight,ef);
}

protected:
virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
					  TPZFMatrix &axes, int var, TPZVec<REAL> &Solout);
public:
virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
    Solution(data.sol,data.dsol,data.axes,var,Solout);
}



};

#endif
