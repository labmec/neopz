/**
 * \file
 * @brief Contains the TPZBlackOil2P3D class which implements a 3D two-phase (oil-water) black-oil flow.
 */
//$Id: pzblackoil2p3d.h,v 1.4 2011-02-04 08:53:02 fortiago Exp $

#ifndef PZBLACKOIL2P3D_H
#define PZBLACKOIL2P3D_H

#include "pzmaterial.h"
#include "pzdiscgal.h"

#ifdef _AUTODIFF

#include "fad.h"


/**
 * @ingroup material
 * @brief Implements a 3D two-phase (oil-water) black-oil flow.
 * @since November 11, 2008
 */
class TPZBlackOil2P3D : public TPZDiscontinuousGalerkin{
	
public:
	
	typedef TFad<4, REAL> BFadREAL;
	
protected:
	
	/** @brief Interpolacao linear */
	void Interpolate(std::map<REAL,REAL> &dados, double x, double &y, double &dy);
	void Interpolate(std::map<REAL,REAL> &dados, BFadREAL x, BFadREAL &y);
	
	/** @brief Simulation time step */
	double fDeltaT;
	
	/** @brief State: one ou one+1 */
	enum EState { ELastState = 0, ECurrentState = 1 };
	
	static EState gState;
	
	void testedoBo();
	
	void testeKrw();
	
public:
	
	static void SetLastState(){ gState = ELastState; }
	static void SetCurrentState(){ gState = ECurrentState; }
	
	/** 
	 * @brief Class constructor 
	 * @param id material id
	 * @param dim problem dimension
	 * @param nstate number of state variables
	 * @param sol constant solution vector
	 */
	TPZBlackOil2P3D(int id, double deltaT);
	
	/** @brief Class destructor */
	~TPZBlackOil2P3D();
	
	/** @brief Copy constructor */
	TPZBlackOil2P3D(const TPZBlackOil2P3D &cp);
	
	/** @brief Defines simulation time step. */
	void SetTimeStep(double timestep){
		this->fDeltaT = timestep;
	}
	
	/** @brief Creates another material of the same type  */
	virtual TPZAutoPointer<TPZMaterial> NewMaterial();
	
	/** @brief Returns problem dimension */
	virtual int Dimension(){ return 3; }
	
	/** @brief Returns number of state variables: oil pressure and oil saturation */
	virtual int NStateVariables(){ return 2; }
	
	/** @brief Contribute method */
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data[in] stores all input data
	 * @param weight[in] is the weight of the integration rule
	 * @param ek[out] is the stiffness matrix
	 * @param ef[out] is the load vector
	 * @param bc[in] is the boundary condition material
	 */
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc);
	
	/** @brief To satisfy base class interface. */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);
	
	/** @brief To satisfy base class interface. */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
	
	
	/**
	 * @name Solution methods
	 * @{
	 */
	
	/** @brief Solution indices of post-processing */
	enum ESolutionVars { ENone = 0, EWaterPressure = 1, EOilPressure, EWaterSaturation, EOilSaturation, EDarcyVelocity };
	
	/** @brief It returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name);
	
	/** 
	 * @brief It returns the number of variables associated with the variable indexed by var.  
	 * @param var is obtained by calling VariableIndex
	 */
	virtual int NSolutionVariables(int var);
	
	/** 
	 * @brief It returns the solution associated with the var index based on
	 * the finite element approximation
	 */
	virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
						  TPZFMatrix &axes, int var, TPZVec<REAL> &Solout);
	/** @} */
	
	/** @brief Fill material data parameter with necessary requirements for the Contribute method. */
	/**
	 * Here, in base class, all requirements are considered as necessary. \n
	 * Each derived class may optimize performance by selecting only the necessary data.
	 */
	virtual void FillDataRequirements(TPZMaterialData &data){
		data.SetAllRequirements(true);
		data.fNeedsNeighborSol = false;
		data.fNeedsNeighborCenter = false;
	}
	
	/** @brief Fill material data parameter with necessary requirements for the ContributeInterface method. */
	/**
	 * Here, in base class, all requirements are considered as necessary. \n
	 * Each derived class may optimize performance by selecting only the necessary data.
	 */
	virtual void FillDataRequirementsInterface(TPZMaterialData &data){
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
	void ViscOleo(double po, double &ViscOleo, double &dViscOleoDpo);
	void ViscOleo(BFadREAL po, BFadREAL &ViscOleo);
	
	/** @brief Capilar pressure. \f$ pc = pc( Sw ) \f$ */
	void PressaoCapilar(double So, double &pc, double &DpcDSo);
	void PressaoCapilar(BFadREAL So, BFadREAL &pc);
	
	/** @brief Porosity. \f$ Phi = Phi( pw ) - fizemos como Phi ( po ) \f$ */
	void Porosidade(double po, double &poros, double &dPorosDpo);
	void Porosidade(BFadREAL po, BFadREAL &poros);
	
	/** @brief Oil density on standard conditions - kg/m3 */
	double RhoOleoSC();
	
	/** @brief Water density on standard conditions - kg/m3 */
	double RhoAguaSC();
	
	/** @brief Gravity */
	double g();
	
	/** @brief Bw = constante   */
	double Bw();
	
	/** @brief Water viscosity (constant) */
	double ViscAgua();
	
	/** @brief Absolute permeability. */
	void K(TPZFMatrix &K);
	
	/** @} */
};

#endif

#endif
