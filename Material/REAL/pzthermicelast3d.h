/**
 * @file
 * @brief Contains the TPZThermicElast3D class which implements a 3D isotropic elasticity material with thermal stress.
 */

#ifndef PZTHERMICELAST3D
#define PZTHERMICELAST3D

#include "pzelast3d.h"

/**
 * @ingroup material
 * @brief This class implements a 3D isotropic elasticity material with thermal stress.
 */
/**
 * The thermal field can be a constant value or a field obtained from a previous simmulation.
 * In case of using a field previously obtained it is input by using TPZReferredCompEl and the field
 * comes in the TPZVec<REAL> sol parameter of Contribute method.
 *  @since May 10, 2006.
 */
class TPZThermicElast3D : public TPZElasticity3D {
	
private:
	
	/** 
	 * @brief Flag indicating a constant field of temperature (fReferred = false) or a field obtained \n
	 * in a previous simmulation (fReferred = true). Default value is true
	 */
	bool fReferred;
	
	/**
	 * @brief When solving problem with a constant field of temperature (fReferred = false) fFinalTemperature \n
	 * is the final temperature of the body. Thermal stress is given by the difference of fFinalTemperature and fRefTemperature. \n
	 * When fReferred = true fFinalTemperature does not make sense and is not used
	 */
	STATE fFinalTemperature;
	
public:
	
	/** @brief Thermal expansion coefficient of the material */
	STATE fThermalCoeff;
	
	/** @brief Temperature of reference. The thermal stress is given by the difference of the temperature field to this reference value */
	STATE fRefTemperature;
	
	/**
	 * @brief Class constructor.
	 * @param nummat - material ID.
	 * @param ThermalCoeff - thermal expansion coefficient of the material
	 * @param RefTemp - temperature of reference
	 * @param E - Young's modulus.
	 * @param poisson - poisson's ratio 
	 * @param force - external forces
	 */ 
	TPZThermicElast3D(int nummat, STATE ThermalCoeff, STATE RefTemp, STATE E, STATE poisson, TPZVec<STATE> &force);
	
	/** @brief Default destructor */
	virtual ~TPZThermicElast3D();
	
	void SetConstantTemperatureField(STATE FinalTemp){
		this->fFinalTemperature = FinalTemp;
		this->fReferred = false;
	}
	
	void SetReferredTemperatureField(){
		this->fFinalTemperature = 0.; //any value, it will not be used
		this->fReferred = true;
	}
	
	bool IsReferredTemperatureField(){
		return this->fReferred;
	}
	
	void ContributeThermalStress(TPZVec<STATE> &sol, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, REAL weight, TPZFMatrix<STATE> &ef);
	
	/** @brief Contribute to stiff matrix and load vector */
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef) override;
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ef) override
	{
		TPZElasticity3D::Contribute(data,weight,ef);
	}
	
protected:
	virtual void Solution(TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol,
						  TPZFMatrix<REAL> &axes, int var, TPZVec<STATE> &Solout) override;
public:
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }
		Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
	}
	
	
    public:
int ClassId() const override;
 
};

#endif
