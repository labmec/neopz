/*
 * @file Helmhotz1D.h
 * @brief Contains the declaration of the TPZHelmholtz1D class.
 */

#ifndef HELMHOLTZ1DHH
#define HELMHOLTZ1DHH

#include "pzmaterial.h"
#include "pzmat1dlin.h"
#include "pzfmatrix.h"
#include "pzvec.h"

/**
 * The associated material is the one-dimensional Helmholtz equations for steady-state oscillations (electromagnetics).
 *
 * \par Variables:
 * \f$ Ez =  \f$ ...\n
 * \f$ Hz = \f$ ... \n
 *
 * \par Equations
 * 
 */
/**
 * @brief Implements the interface for Helmholtz equation one-dimensional
 * @ingroup material
 */
class TPZHelmholtz1D  : public TPZMat1dLin
{
	/** @brief Pointer to alpha function, coefficient of the first derivative of u */
    TPZAutoPointer<TPZFunction> fAlpha;
	/** @brief Pointer to beta function, variable coefficient of the u */
    TPZAutoPointer<TPZFunction> fBeta;
        /** @brief Pointer to phi function */
    TPZAutoPointer<TPZFunction> fPhi;
    
public:
	/** @brief Simple constructor with material id and dimension of the spatial domain */
	TPZHelmholtz1D(int nummat,int dim);
	/** @brief Copy constructor */
	TPZHelmholtz1D(const TPZHelmholtz1D &cp);
	
	/** @brief Default destructor */
	~TPZHelmholtz1D();

    /** @brief Returns the name of the material */
    virtual std::string Name() { return "Helmholtz1D"; }

	
	// MUST TO BE IMPLEMENTED
	
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() { return 1; }

    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() { return 2; }

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
//    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc);
	
	/** @brief Sets the variable coefficient alpha */
	void SetAlphaFunction(TPZAutoPointer<TPZFunction> falpha) {
		fAlpha = falpha;
	}
	/** @brief Sets the variable coefficient beta */
	void SetBetaFunction(TPZAutoPointer<TPZFunction> fbeta) {
		fBeta = fbeta;
	}
        
        void SetPhiFunction(TPZAutoPointer<TPZFunction> fphi){
                fPhi = fphi;
        }
};

#endif