/*
 * @file Helmhotz1D.h
 * @brief Contains the declaration of the TPZHelmholtz1D class.
 */

#ifndef HELMHOLTZCOMPLEX1DHH
#define HELMHOLTZCOMPLEX1DHH

#include "TPZMaterial.h"
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
class TPZHelmholtzComplex1D  : public TPZMat1dLin
{
    /** @brief Pointer to alpha function, coefficient of the first derivative of u */
    TPZAutoPointer<TPZFunction<STATE> > fAlpha;
	/** @brief Pointer to beta function, variable coefficient of the u */
    TPZAutoPointer<TPZFunction<STATE> > fBeta;
        /** @brief Pointer to phi function */
    TPZAutoPointer<TPZFunction<STATE> > fPhi;
    
public:
    /** @brief Simple constructor with material id and dimension of the spatial domain */
    TPZHelmholtzComplex1D(int nummat, int dim);
    /** @brief Copy constructor */
    TPZHelmholtzComplex1D(const TPZHelmholtzComplex1D &cp);
    /** @brief Default destructor */
    ~TPZHelmholtzComplex1D();

    /** @brief Returns the name of the material */
    virtual std::string Name() { return "HelmholtzComplex1D"; }

    // MUST TO BE IMPLEMENTED
	
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() { return 1; }

    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() { return 1; }
    
    /** @brief Returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name);
    
    virtual int NSolutionVariables(int var);
    
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
    /** @brief Sets the variable coefficient alpha */
    void SetAlphaFunction(TPZAutoPointer<TPZFunction<STATE> > falpha) {
		fAlpha = falpha;
    }
    /** @brief Sets the variable coefficient beta */
    void SetBetaFunction(TPZAutoPointer<TPZFunction<STATE> > fbeta) {
		fBeta = fbeta;
    }

    void SetPhiFunction(TPZAutoPointer<TPZFunction<STATE> > fphi) {
		fPhi = fphi;
    }
};

#endif