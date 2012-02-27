/**
 * \file
 * @brief Contains the TPZBurger class which implements a linear convection equation using a burger flux.
 */
// -*- c++ -*-

//$Id: pzburger.h,v 1.7 2009-10-09 14:53:27 fortiago Exp $

#ifndef BURGERH
#define BURGERH

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

/**
 * @ingroup material
 * @brief This class implements a linear convection equation using
 * a burger flux instead of the linear flux.
 * @author Tiago Forti
 */
/** 
 * It's been developed for a Petrobras report where the water temperature is transported \n
 * into the reservoir following the Darcy's velocity field. \n
 * I apologise for the class name which is not exact.
 */
class TPZBurger : public TPZMatPoisson3dReferred {
	
public:
	/** @brief Enum for stabilization scheme into of the element */
	enum EStabilizationScheme{ESUPG = 1, EGRADIENT = 2};
	
	static int gStabilizationScheme;
	/** @brief Constructor with id of material and dimension of the space */
	TPZBurger(int nummat, int dim);
	/** @brief Copy constructor */
	TPZBurger(const TPZBurger &cp);
	/** @brief Destructor */
	virtual ~TPZBurger();
	
	//   virtual int HasForcingFunction() {
	//     return true;
	//   }
	
	bool IsReferred(){ return this->fIsReferred;}
	
	void SetReferred(bool Is){ this->fIsReferred = Is; }
	
	REAL fSolRef;
	
	/** @name Contribute methods */
	/** @{ */
	
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef) {
        int numbersol = data.sol.size();
        if(numbersol != 1) 
        {
            DebugStop();
        }
		if (TPZBurger::gStabilizationScheme == ESUPG){
			this->ContributeSUPG(data.x,data.jacinv,data.sol[0],data.dsol[0],weight,data.axes,data.phi,data.dphix,ek,ef);
		}
		if (TPZBurger::gStabilizationScheme == EGRADIENT){
			this->ContributeGradStab(data.x,data.jacinv,data.sol[0],data.dsol[0],weight,data.axes,data.phi,data.dphix,ek,ef);
		}
	}//Contribute
	
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
							TPZFMatrix &ef)
	{
		TPZMatPoisson3dReferred::Contribute(data,weight,ef);
	}
	
	void ContributeGradStab(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
							TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
	void ContributeSUPG(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
						TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ek,
							  TPZFMatrix &ef,
							  TPZBndCond &bc);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix &ek,
									 TPZFMatrix &ef);
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix &ek,
									   TPZFMatrix &ef,
									   TPZBndCond &bc);
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ef,
							  TPZBndCond &bc)
	{
		TPZMatPoisson3dReferred::ContributeBC(data,weight,ef,bc);
	}
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix &ef)
	{
		TPZMatPoisson3dReferred::ContributeInterface(data,dataleft,dataright,weight,ef);
	}
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix &ef,
									   TPZBndCond &bc)
	{
    	TPZMatPoisson3dReferred::ContributeBCInterface(data,dataleft,weight,ef,bc);
	}
	
	/** @} */
	
protected:
    bool fIsReferred;
	
};

#endif
