/**
 * @file
 * @brief Contains TPZBiharmonicEstimator class estimates error to biharmonic problem.
 */

#ifndef TPZBIHARMONICESTIMATOR_H
#define TPZBIHARMONICESTIMATOR_H
#include "pzbiharmonic.h"

/**
 * @brief Estimates error to biharmonic problem. Also computes the contributions on elements and interfaces. \ref analysis "Analysis"
 * @author Joao Luis Goncalves
 * @since Maio 16, 2008
 * @ingroup material
 * @note This class seems as material class, must to be put into the material group (Jorge??)
 */
class TPZBiharmonicEstimator: public  TPZBiharmonic
{
private:
	
	/** @brief Attributes required for goal oriented error estimation validation */
	void (*fPrimalExactSol)(TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv);
	void (*fDualExactSol)(TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv);
	
    public:
        int ClassId() const override;

	/** @brief Constructor */
    TPZBiharmonicEstimator(int nummat, STATE f);
	/** @brief Destructor */
    ~TPZBiharmonicEstimator();
	
	/** @brief Set the pointer of the solution function */
    void SetExactSolutions(void (*fp)(TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv),
                           void (*fd)(TPZVec<REAL> &locdual,TPZVec<STATE> &valdual,TPZFMatrix<STATE> &derivdual));
	
	/** @brief Returns the number of norm errors. Default is 3: energy, L2 and H1. */
	virtual int NEvalErrors() override {return 4;}
	
	/** @brief Implements integration of the internal part of an error estimator. */
	/** It performs nk[0] += weight * ( residuo(u) *(Z1-Z) ); \n
	 * where u is the current solution and Z and Z1 are the dual solution. */
	virtual void ContributeErrorsDual(TPZMaterialData &data,
									  REAL weight,
									  TPZVec<REAL> &nk);
	
	virtual void ContributeErrors(TPZMaterialData &data,
								  REAL weight,
								  TPZVec<REAL> &nk,
								  int &errorid)
	{
		if (errorid == 0) this->ContributeErrorsDual(data,weight,nk);
		if (errorid == 2) this->ContributeErrorsSimple(data,weight,nk);
	}
	
	/** @brief Implements integration of the interface part of an error estimator. */
	/** It performs \f$ nk[0] += weight * ( residuo(u )*(Z1-Z) ) \f$ ; \n
	 * where u is the current solution and Z and Z1 are the dual solution. */
	virtual void ContributeInterfaceErrors(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight,
										   TPZVec<STATE> &nkL, TPZVec<STATE> &nkR,
										   int &errorid)
	{
		if (errorid == 0) this->ContributeInterfaceErrorsDual(data,dataleft,dataright,weight,nkL,nkR);
		if (errorid == 2) this->ContributeInterfaceErrorsSimple(data,dataleft,dataright,weight,nkL,nkR);
	}
	
	virtual void ContributeInterfaceErrorsDual(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
											   REAL weight,
											   TPZVec<STATE> &nkL, 
											   TPZVec<STATE> &nkR);
	
	/** @brief Implements integration of the boundary interface part of an error estimator. */
	/** 
	 * It performs \f$ nk[0] += weight * ( residuo(u ) * (Z1-Z) ) \f$ ; \n
	 * where u is the current solution and Z and Z1 are the dual solution.
	 */
	virtual void ContributeInterfaceBCErrorsDual(TPZMaterialData &data, TPZMaterialData &dataleft,
												 REAL weight,
												 TPZVec<STATE> &nk, 
												 TPZBndCond &bc);
	
	virtual void ContributeInterfaceBCErrors(TPZMaterialData &data, TPZMaterialData &dataleft,
											 REAL weight,
											 TPZVec<STATE> &nk,
											 TPZBndCond &bc,
											 int &errorid)
	{
		if (errorid == 0) this->ContributeInterfaceBCErrorsDual(data,dataleft,weight,nk,bc);
		if (errorid == 2) this->ContributeInterfaceBCErrorsSimple(data,dataleft,weight,nk,bc);
	}
	
	virtual void ContributeErrorsSimple(TPZMaterialData &data,
										REAL weight,
										TPZVec<REAL> &nk);
	
	
	virtual void ContributeInterfaceErrorsSimple(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
												 REAL weight,
												 TPZVec<STATE> &nkL, 
												 TPZVec<STATE> &nkR);
	
	
	virtual void ContributeInterfaceBCErrorsSimple(TPZMaterialData &data, TPZMaterialData &dataleft,
												   REAL weight,
												   TPZVec<STATE> &nk,
												   TPZBndCond &bc);

    /**
	 * @brief Compute the error due to the difference between the interpolated flux \n
	 * and the flux computed based on the derivative of the solution
	 */
protected:
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u, TPZFMatrix<STATE> &dudx,
				TPZFMatrix<REAL> &axes,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,
				TPZVec<REAL> &values) override; 
public:
	/** @brief Kernel of the functional */
	void Psi(TPZVec<REAL> &x, TPZVec<STATE> &pisci);
	
	/** @brief Computes the primal and dual exact error */
	void OrderSolution(TPZMaterialData &data);
	/** @brief Computes the primal and dual exact error over dataright considering the solution in data */
	void OrderSolutionRight(TPZMaterialData &data, TPZMaterialData &dataright);
	/** @brief Computes the primal and dual exact error over dataleft considering the solution in data */
	void OrderSolutionLeft(TPZMaterialData &data, TPZMaterialData &dataleft);
	
};

#endif
