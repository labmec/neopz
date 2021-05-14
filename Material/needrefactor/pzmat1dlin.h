/**
 * @file
 * @brief Contains the TPZMat1dLin class which implements a one-dimensional linear problem.
 */

#ifndef MAT1DLINHPP
#define MAT1DLINHPP

#include <iostream>


#include "pzfmatrix.h"
#include "TPZMaterial.h"
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
	
	
    TPZMat1dLin(int num) : TPZRegisterClassId(&TPZMat1dLin::ClassId), TPZMaterial(num) , fXk(1,1,0.), fXc(1,1,0.), fXb(1,1,0.), fXf(1,1,0.) {
	}
    
    TPZMat1dLin(const TPZMat1dLin &copy) : TPZRegisterClassId(&TPZMat1dLin::ClassId),
    TPZMaterial(copy), fXk(copy.fXk), fXc(copy.fXc), fXb(copy.fXb), fXf(copy.fXf)
    {
        
    }
    
    TPZMat1dLin &operator=(const TPZMat1dLin &copy)
    {
        TPZMaterial::operator=(copy);
        fXk = copy.fXk;
        fXc = copy.fXc;
        fXb = copy.fXb;
        fXf = copy.fXf;
        return *this;
    }
    
    virtual TPZMaterial *NewMaterial()  override 
    {
        return new TPZMat1dLin(*this);
    }
	
	virtual int NStateVariables() const override {
        return fXk.Rows();
    }
	
	virtual int Dimension() const  override {
        return 1;
    }
	
	void Print(std::ostream & out) override ;
	
	void SetMaterial(TPZFMatrix<STATE> &xkin,TPZFMatrix<STATE> &xcin,TPZFMatrix<STATE> &xbin,TPZFMatrix<STATE> &xfin){
		fXk = xkin;
		fXc = xcin;
		fXb = xbin;
		fXf = xfin;
	}
	
	virtual std::string Name()  override { return "TPZMat1dLin"; }
	
	/** @brief Computes contribution to the stiffness matrix and right hand
	 * side at an integration point*/
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef) override ;
	
	/** @brief Computes contribution to the stiffness matrix and right hand
	 * side at an integration point*/
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ef) override
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	/** @brief Computes contribution to the stiffness matrix and right hand
	 * side at the integration point of a boundary*/
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override ;
	
	/** @brief Computes contribution to the stiffness matrix and right hand
	 * side at the integration point of a boundary*/
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override 
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
protected:
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, 
						TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;
    public:
virtual int ClassId() const override ;
 
};

#endif
