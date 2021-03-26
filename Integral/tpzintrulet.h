/**
 * @file
 * @brief Contains the TPZIntRuleT class which defines integration rule for triangles based on Linbo Zhang's paper( http://www.jstor.org/stable/43693493).
 */

#ifndef TPZINTRULET_H
#define TPZINTRULET_H

#include "pzreal.h"
#include "pzvec.h"
#include "pzmanvector.h"

/**
 * @ingroup integral
 * @brief Integration rule (points and weights) for triangles. \ref integral "Numerical Integration"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZIntRuleT {
	
	/** @brief The list can to access the constructor of the current class. */
	friend class TPZIntRuleList;
	
	/** @brief Number of integration points for this object */
    int	   fNumInt;
	/** @brief Location of the integration point Ksi */
    TPZManVector<long double>	fLocationKsi;
	/** @brief Location of the integration point Eta */
    TPZManVector<long double>	fLocationEta;
	/** @brief Weight of the integration point */
    TPZManVector<long double>	fWeight;
    /// the polynomial order this integration rule can integrate
    int fOrder;

	/**
	 * @brief Constructor of integration rule for triangle.
	 * @param order Order of the polinomial will be exactly integrated.
	 */
	TPZIntRuleT(int order);
    
    TPZIntRuleT(const TPZIntRuleT &copy) : fNumInt(copy.fNumInt), fLocationKsi(copy.fLocationKsi), fLocationEta(copy.fLocationEta),
        fWeight(copy.fWeight), fOrder(copy.fOrder)
    {
        
    }
    
    TPZIntRuleT &operator=(const TPZIntRuleT &copy)
    {
        fNumInt = copy.fNumInt;
        fLocationKsi = copy.fLocationKsi;
        fLocationEta = copy.fLocationEta;
        fWeight = copy.fWeight;
        fOrder = copy.fOrder;
        return *this;
    }
	/** @brief Default destructor. */
	~TPZIntRuleT();

	/**
	 * @brief Computes the cubature rules following the symmetric construction presented at Linbo Zhang article.
	 * @param order Order of the polinomial will be exactly integrated.
	 */
	int ComputingSymmetricCubatureRule(int order);

	/**
	 * @brief Transforms barycentric coordinates (3 component) of the point in triange in cartesian coordinates (2 component).
	 * @param baryvec Vector of barycentric coordinates (3 components) for each point
	 * @param weightvec Vector of weights related for each point in baryvec vector
	 */
	void TransformBarycentricCoordInCartesianCoord(long double baryvec[],long double weightvec[]);

public:
    enum {NRULESTRIANGLE_ORDER = 21};

	/** @brief Returns number of integration points */
    int NInt() const { return fNumInt;}

	/** @brief Returns location of the ith point */
    void Loc(int i, TPZVec<REAL> &pos) const;

	/** @brief Return weight for the ith point */
    REAL W(int i) const;
    
    /// Order associated with the integration rule
    int Order() const
    {
        return fOrder;
    }

};

#endif
