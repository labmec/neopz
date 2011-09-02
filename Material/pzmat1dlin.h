/**
 * \file
 * @brief Contains the TPZMat1dLin class which implements a one-dimensional linear problem.
 */
#ifndef MAT1DLINHPP
#define MAT1DLINHPP

#include <iostream>


#include "pzfmatrix.h"
#include "pzmaterial.h"
#include "pzvec.h"

struct TPZElementMatrix;
class TPZBndCond;
template<class T>
class TPZVec;

/**
 * @ingroup material
 * @brief Implements a one dimensional linear problem.
 */
class TPZMat1dLin : public TPZMaterial{
	
	TPZFMatrix		fXk,fXc,fXb,fXf;
	
	public :
	
	
    TPZMat1dLin(int num) : TPZMaterial(num) , fXk(), fXc(), fXb(), fXf() {
	}
	
	virtual int NStateVariables() { return fXk.Rows(); }
	
	int Dimension() { return 1;}
	
	void Print(std::ostream & out);
	
	void SetMaterial(TPZFMatrix &xkin,TPZFMatrix &xcin,TPZFMatrix &xbin,TPZFMatrix &xfin){
		fXk = xkin;
		fXc = xcin;
		fXb = xbin;
		fXf = xfin;
	}
	
	virtual std::string Name() { return "TPZMat1dLin"; }
	
	int NFluxes() { return NStateVariables(); }
	
	/** @brief Computes contribution to the stiffness matrix and right hand
	 * side at an integration point*/
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix &ek,
							TPZFMatrix &ef);
	
	/** @brief Computes contribution to the stiffness matrix and right hand
	 * side at an integration point*/
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	/** @brief Computes contribution to the stiffness matrix and right hand
	 * side at the integration point of a boundary*/
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ek,
							  TPZFMatrix &ef,
							  TPZBndCond &bc);
	
	/** @brief Computes contribution to the stiffness matrix and right hand
	 * side at the integration point of a boundary*/
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ef,
							  TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &fl);
	
	virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
						TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
};

#endif
