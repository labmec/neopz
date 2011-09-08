/**
 * @file
 * @brief DEPRECATED FILE. Contains the TPZNonLinElliptic class.
 */
//
// C++ Interface: pznonlinelliptic
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef PZNONLINELLIPTIC_H
#define PZNONLINELLIPTIC_H

#include "pzmaterial.h"

/**
 * @deprecated DEPRECATED non linear elliptic material CLASS. 
 * @brief Implements a non linear elliptic equation. \f$ Laplac( u ) + div(V u) + Sigma u = u^2 \f$
 * @author Roberto Heiderich.
 * @since Feb 23, 2005
 */
class TPZNonLinElliptic : public TPZMaterial
{
public:
	
    /** 
	 * @brief Constructor.
     * @param id = material identification
     * @param dimension = set problem dimension
     */
    TPZNonLinElliptic(int id, int dimension);
	
    /** @brief Defaul destructor. */
    ~TPZNonLinElliptic();
    
    /** @brief Returns the integrable dimension of the material*/
    virtual int Dimension(){ return fDim;}
    
    /** @brief Returns the number of state variables associated with the material*/
    virtual int NStateVariables(){return 1;}
    
    /** @brief Computes contribution to the stiffness matrix and right hand
     * side at an integration point*/
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix &ek,
							TPZFMatrix &ef);
	
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix &ef);
	
    /** @brief Computes contribution to the stiffness matrix and right hand
     * side at the integration point of a boundary*/
    virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight, 
							  TPZFMatrix &ek,
							  TPZFMatrix &ef,
							  TPZBndCond &bc);
	void Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,TPZVec<REAL> &sol, TPZFMatrix &dsol,
					REAL weight,TPZFMatrix &axes,TPZFMatrix &phi, TPZFMatrix &dphi,
					TPZFMatrix &ek, TPZFMatrix &ef);
	
	void Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,TPZVec<REAL> &sol, TPZFMatrix &dsol,
					REAL weight,TPZFMatrix &axes,TPZFMatrix &phi, TPZFMatrix &dphi,TPZFMatrix &ef);
	
    void SetParameters(REAL D, TPZVec<REAL> &V, REAL Sigma, REAL LambdaDivK, REAL F);
    
    void GetParameters(REAL &D, TPZVec<REAL> &V, REAL &Sigma, REAL &LambdaDivK, REAL &F);    
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &/*axes*/,int var,TPZVec<REAL> &Solout);
	
	virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
						TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
						TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
	
private:
	
    /** @brief Problem dimension. */
    int fDim;
    
    /** 
	 * @brief Populational dispersion coefficient. \n
     * Multiplies laplacian operator.
     */
    REAL fCoeffD;
    
    /** @brief Inside divergent operator. */
    TPZVec<REAL> fConvDir;
	
    /** @brief Reaction coefficient. */
	/** Hostilidade do meio - taxa de crescimento. */    
    REAL fSigma;
    
    /** @brief Increment rate divided by the support capacity. */
    REAL fLambdaDivK;
    
    /** @brief Source constant term. */
    REAL fSource;   
    
};

#endif
