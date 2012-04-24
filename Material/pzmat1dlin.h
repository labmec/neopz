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
	
	TPZFMatrix<STATE>		fXk,fXc,fXb,fXf;
	
	public :
	
	
    TPZMat1dLin(int num) : TPZMaterial(num) , fXk(), fXc(), fXb(), fXf() {
	}
	
	virtual int NStateVariables() { return fXk.Rows(); }
	
	int Dimension() { return 1;}
	
	void Print(std::ostream & out);
	
	void SetMaterial(TPZFMatrix<STATE> &xkin,TPZFMatrix<STATE> &xcin,TPZFMatrix<STATE> &xbin,TPZFMatrix<STATE> &xfin){
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
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef);
	
	/** @brief Computes contribution to the stiffness matrix and right hand
	 * side at an integration point*/
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	/** @brief Computes contribution to the stiffness matrix and right hand
	 * side at the integration point of a boundary*/
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc);
	
	/** @brief Computes contribution to the stiffness matrix and right hand
	 * side at the integration point of a boundary*/
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &fl);
	
	virtual void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
						TPZVec<REAL> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<STATE> &values);
};

#endif
