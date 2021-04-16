/**
 * @file
 * @brief Contains the TPZNonLinBiharmonic class which implements a discontinuous Galerkin formulation for the non-linear bi-harmonic equation.
 */

#ifndef TPZNONLINBIHARMONICHPP
#define TPZNONLINBIHARMONICHPP

#include <iostream>
#include "TPZMaterial.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief This class implements discontinuous Galerkin formulation for the non-linear bi-harmonic equation.
 * @since Jan 31, 2005
 * @author Igor Mozolevski e Paulo Bosing
 */
class TPZNonLinBiharmonic : public TPZMaterial {
	
private:
	STATE  fXf;
	
	public :
	
	static STATE gLambda1, gLambda2, gSigmaA,gSigmaB, gL_alpha, gM_alpha, gL_betta,
	gM_betta, g_teta, Re;
	static int NorP;
	
	/** @brief Inicialisation of biharmonic material */
	TPZNonLinBiharmonic(int nummat, STATE f);
	
	virtual ~TPZNonLinBiharmonic();
	
	/** @brief Returns the number of norm errors. Default is 3: energy, L2,  H1, semi-norm H2 and H2. */
	virtual int NEvalErrors()  override {return 8;}
	
	void SetMaterial(STATE &xfin){
		fXf = xfin;
	}
	
	virtual int Dimension() const  override { return 2;}
	
	/** @brief Returns one because of scalar problem */
	virtual int NStateVariables() const override {
		return 1;
	};
	
	virtual void Print(std::ostream & out) override;
	
	virtual std::string Name()  override { return "TPZBiharmonic"; }

	/**
	 * @name Contribute methods from weak formulation
	 * @{
	 */
	 
	/** @brief Implements integral over  element's volume */
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;
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
							  TPZBndCond &bc) override;
	
	/** @brief Implements boundary conditions for continuous Galerkin */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}

	/** @} */
	
	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;
	
protected:
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
public:
	
	virtual void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<STATE> &Solout)
	{
		TPZMaterial::SolutionDisc(data,dataleft,dataright,var,Solout);
	}
protected:
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;
public:

	/**
	 * @name Contribute interface methods
	 * @{
	 */
	 
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix<STATE> &ek,
									 TPZFMatrix<STATE> &ef) override;
	
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ek,
									   TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc) override;
	
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
    public:
virtual int ClassId() const override;
 
	/** @} */
	
};

#endif
