/**
 * \file
 * @brief DEPRECATED FILE. This file contains the TPZConservationLawDEP class which implements the interface for conservation laws
 */
#ifndef CONSERVATIONLAWHPP
#define CONSERVATIONLAWHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"

/**
 * @deprecated DEPRECATED conservation law material CLASS.
 * @brief Implements the interface for conservation laws, keeping track of the timestep as well. \n THIS CLASS IS DEPRECATED BY TPZConservationLaw2.
 */
class TPZConservationLawDEP  : public TPZMaterial {
	
	int fDim;
	REAL fTimeStep;
	
	public :
	
	REAL fDelta;
	
	TPZConservationLawDEP(int nummat,REAL delta_t,int dim);
	
	/** @brief copy constructor*/
	TPZConservationLawDEP(TPZConservationLawDEP &copy);
	
	/** @brief To create another material of the same type*/
	TPZAutoPointer<TPZMaterial> NewMaterial();
	
	~TPZConservationLawDEP(){};
	
	/** @brief compute the boundary condition left solution */
	virtual void ComputeSolLeft(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcleft);
	
	/** @brief compute the boundary condition right solution */
	virtual void ComputeSolRight(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcright);
	
	/** @brief termodinamic pressure determined by the law of an ideal gas */
	virtual REAL Pressure(TPZVec<REAL> &U);
	
	virtual REAL Gamma();
	
	//virtual int IntegrationDegree() = 0;
	
	//virtual void SetIntegDegree(int degree) = 0;
	
	virtual void SetDelta(REAL delta){fDelta = delta;}
	
	virtual void SetDeltaTime(REAL maxveloc,REAL deltax,int degree);
	
	REAL Delta();
	
	virtual void SetTimeStep(REAL timestep){fTimeStep = timestep;}
	
	virtual REAL TimeStep(){return fTimeStep;}
	
	int Dimension() { return fDim;}
	
	int NStateVariables();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZConservationLawDEP"; }
	
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef);
	
	virtual void ContributeInterface(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix &ek,
                                     TPZFMatrix &ef);
	
	virtual void ContributeBC(TPZMaterialData &data,
                              REAL weight,
                              TPZFMatrix &ek,
                              TPZFMatrix &ef,
                              TPZBndCond &bc);
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ef,
							  TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ef);
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 1;}
	
protected:
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
public:
	
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
				TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
};


// inline void TPZConservationLawDEP::ContributeInterface(TPZMaterialData &data,
//                                                       REAL ,
//                                                       TPZFMatrix &){
//   PZError << "TPZConservationLawDEP::ContributeInterface it would never have to be called\n";
// }

inline void TPZConservationLawDEP::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){
	PZError << "TPZConservationLawDEP::Solution it would never have to be called\n";
}

inline void TPZConservationLawDEP::SetDeltaTime(REAL maxveloc,REAL deltax,int degree){
	PZError << "TPZConservationLawDEP::SetDeltaTime it would never have to be called\n";
}

inline REAL TPZConservationLawDEP::Pressure(TPZVec<REAL> &U) {
	PZError << "TPZConservationLawDEP::Pressure it would never have to be called\n";
	return 0.0;
}

inline REAL TPZConservationLawDEP::Gamma() {
	PZError << "TPZConservationLawDEP::Gamma it would never have to be called\n";
	return 0.0;
}

inline void TPZConservationLawDEP::ComputeSolLeft(TPZVec<REAL> &,TPZVec<REAL> &,TPZVec<REAL> &,TPZBndCond *){
	PZError << "TPZConservationLawDEP::ComputeSolLeft it would never have to be called\n";
}

inline void TPZConservationLawDEP::ComputeSolRight(TPZVec<REAL> &,TPZVec<REAL> &,TPZVec<REAL> &,TPZBndCond *){
	PZError << "TPZConservationLawDEP::ComputeSolRight it would never have to be called\n";
}
#endif
