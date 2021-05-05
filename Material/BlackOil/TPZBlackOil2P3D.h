/**
 * \file
 * @brief Contains the TPZBlackOil2P3D class which implements a 3D two-phase (oil-water) black-oil flow.
 */

#ifndef PZBLACKOIL2P3D_H
#define PZBLACKOIL2P3D_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatInterfaceSingleSpace.h"
#include "fad.h"

/**
 * @ingroup material
 * @brief Implements a 3D two-phase (oil-water) black-oil flow.
 * @since November 11, 2008
 */
class TPZBlackOil2P3D : public
TPZMatBase<STATE,
		   TPZMatSingleSpaceT<STATE>,
		   TPZMatInterfaceSingleSpace<STATE>>
{
    using TBase = TPZMatBase<STATE,
							 TPZMatSingleSpaceT<STATE>,
							 TPZMatInterfaceSingleSpace<STATE>>;
	
public:
	using BFadREAL =  TFad<4, REAL>;
	static void SetLastState(){ gState = ELastState; }
	static void SetCurrentState(){ gState = ECurrentState; }
	
	/** 
	 * @brief Class constructor 
	 * @param id material id
	 * @param deltaT time step
	 */
	TPZBlackOil2P3D(int id, double deltaT);
	
	/** @brief Defines simulation time step. */
	void SetTimeStep(double timestep){
		this->fDeltaT = timestep;
	}
	
	/** @brief Creates another material of the same type  */
	TPZMaterial * NewMaterial() const override;
	
	/** @brief Returns problem dimension */
	int Dimension() const override { return 3; }
	
	/** @brief Returns number of state variables: oil pressure and oil saturation */
	int NStateVariables() const override{ return 2; }
	
	/** @brief Contribute method */
	void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight,
                    TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

	void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight,
                      TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;
   	
    void ContributeInterface(const TPZMaterialDataT<STATE> &data,
                             const TPZMaterialDataT<STATE> &dataleft,
                             const TPZMaterialDataT<STATE> &dataright,
                             REAL weight,
                             TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
	
    void ContributeBCInterface(const TPZMaterialDataT<STATE> &data,
                               const TPZMaterialDataT<STATE> &dataleft,
                               REAL weight,
                               TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                               TPZBndCondT<STATE> &bc) override;
	
	
	/**
	 * @name Solution methods
	 * @{
	 */
	
	/** @brief Solution indices of post-processing */
	enum ESolutionVars { ENone = 0, EWaterPressure = 1, EOilPressure, EWaterSaturation, EOilSaturation, EDarcyVelocity };

	int VariableIndex(const std::string &name) const override;

	int NSolutionVariables(int var) const override;
	
    void Solution(const TPZMaterialDataT<STATE>&data,
                  int var, TPZVec<STATE> &Solout) override;
    
    void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                      const TPZMaterialDataT<STATE> &dataleft,
                      const TPZMaterialDataT<STATE> &dataright,
                      const int var,
                      TPZVec<STATE> &Solout) override;
	/** @} */

    void GetSolDimensions(uint64_t &u_len,
                          uint64_t &du_row,
                          uint64_t &du_col) const override;
    
	void FillDataRequirements(TPZMaterialData &data) const override {
		data.SetAllRequirements(true);
		data.fNeedsNeighborSol = false;
		data.fNeedsNeighborCenter = false;
	}
    
    void FillDataRequirementsInterface(TPZMaterialData &data) const override {
		data.SetAllRequirements(true);
		data.fNeedsSol = false;
	}
	/**
	 * @name Getting Data
	 * @{
	 */
	
	/**
	 * @brief Oil relative permeability.
	 * \f$ Kro = Kro( Sw ) \f$
	 */
	void Kro(double So, double &Kro, double &dKroSo);
	void Kro(BFadREAL So, BFadREAL &Kro);
	
	/**
	 * @brief Water relative permeability.
	 * \f$ Krw = Krw( Sw ) \f$
	 */
	void Krw(double So, double &Krw, double &dKrwSo);
	void Krw(BFadREAL So, BFadREAL &Krw);
	
	/** 
	 * @brief \f$ Bo = Bo( po ) \f$
	 * @param po initial pressure
	 */
	void Bo(double po, double &Bo, double &dBoDpo);
	void Bo(BFadREAL po, BFadREAL &Bo);
	
	/** @brief Oil viscosity. \f$ ViscOleo = ViscOleo( po ) \f$ */
	void OilVisc(double po, double &ViscOleo, double &dViscOleoDpo);
	void OilVisc(BFadREAL po, BFadREAL &ViscOleo);
	
	/** @brief Capilar pressure. \f$ pc = pc( Sw ) \f$ */
	void CapilarPressure(double So, double &pc, double &DpcDSo);
	void CapilarPressure(BFadREAL So, BFadREAL &pc);
	
	/** @brief Porosity. \f$ Phi = Phi( pw ) - fizemos como Phi ( po ) \f$ */
	void Porosity(double po, double &poros, double &dPorosDpo);
	void Porosity(BFadREAL po, BFadREAL &poros);
	
	/** @brief Oil density on standard conditions - kg/m3 */
	STATE OilRhoSC() const;
	
	/** @brief Water density on standard conditions - kg/m3 */
	STATE WaterRhoSC() const;
	
	/** @brief Gravity */
	STATE g() const;
	
	/** @brief Bw = constante   */
	STATE Bw() const;
	
	/** @brief Water viscosity (constant) */
	STATE WaterVisc() const;
	
	/** @brief Absolute permeability. */
	void K(TPZFMatrix<REAL> &K);
	
	/** @} */
protected:
	
	/** @brief Interpolacao linear */
	void Interpolate(std::map<REAL,REAL> &dados, double x, double &y, double &dy);
	void Interpolate(std::map<REAL,REAL> &dados, BFadREAL x, BFadREAL &y);
	
	/** @brief Simulation time step */
	STATE fDeltaT;
	
	/** @brief State: one ou one+1 */
	enum EState { ELastState = 0, ECurrentState = 1 };
	
	static EState gState;
};

#endif

