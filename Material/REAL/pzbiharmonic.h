/**
 * \file
 * @brief Contains the TPZBiharmonic class which implements a discontinuous Galerkin formulation for the bi-harmonic equation.
 */

#ifndef  TPZBIHARMONICHPP
#define TPZBIHARMONICHPP

#include <iostream>
#include "TPZMaterial.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief Implements discontinuous Galerkin formulation for the bi-harmonic equation.
 * @since Nov 27, 2003
 * @author Igor Mozolevski e Paulo Bosing
 */
class TPZBiharmonic : public TPZMaterial {
	
protected:
	REAL  fXf;
	
	public :
	
	static REAL gLambda1, gLambda2, gSigmaA,gSigmaB, gL_alpha, gM_alpha, gL_betta, gM_betta;
	
	/** @brief Inicialisation of biharmonic material */
	TPZBiharmonic(int nummat, REAL f);
	/** @brief Destructor */
	virtual ~TPZBiharmonic();
	
	/** @brief Returns the number of norm errors. Default is 3: energy, L2,  H1, semi-norm H2 and H2. */
	virtual int NEvalErrors() override {return 8;}
	
	void SetMaterial(REAL &xfin) {
		fXf = xfin;
	}
	/** @brief Model dimension */
	int Dimension() const  override { return 2;}
	
	/** @brief Returns one because of scalar problem */
	virtual int NStateVariables() const override {
		return 1;
	};
	
	virtual void Print(std::ostream & out) override ;
	
	virtual std::string Name() override  { return "TPZBiharmonic"; }
	
	/** @brief Implements integral over  element's volume */
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
							TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override ;
	/** @brief Implements integral over  element's volume */
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
							TPZFMatrix<STATE> &ef) override 
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	/** @brief Implements boundary conditions for continuous Galerkin */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override ;
	
	/** @brief Implements boundary conditions for continuous Galerkin */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override 
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int VariableIndex(const std::string &name) override ;
	
	virtual int NSolutionVariables(int var) override ;
        
        public:
virtual int ClassId() const override ;

	
protected:
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override ;

public:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<STATE> &Solout)
	{
		TPZMaterial::SolutionDisc(data,dataleft,dataright,var,Solout);
	}
	
    /**
	 * @brief Compute the error due to the difference between the interpolated flux 
	 * and the flux computed based on the derivative of the solution
	 */
protected:
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;
public:
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix<STATE> &ek,
									 TPZFMatrix<STATE> &ef) override ;
	
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ek,
									   TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc) override ;
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix<STATE> &ef) override 
	{
		TPZMaterial::ContributeInterface(data,dataleft,dataright,weight,ef);
	}
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc) override 
	{
		TPZMaterial::ContributeBCInterface(data,dataleft,weight,ef,bc);
	}

};

#endif
