/**
 * @file
 * @brief Contains the TPZIntRuleT class which defines integration rule for triangles.
 */

#ifndef TPZINTRULET_H
#define TPZINTRULET_H

#include "pzreal.h"
#include "pzvec.h"

/**
 * @ingroup integral
 * @brief Integration rule (points and weights) for triangles. \ref integral "Numerical Integration"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZIntRuleT {
	
	/** @brief The list can to access the constructor of the current class. */
	friend class TPZIntRuleList;
	
	/** @brief Number of integration points for this object */
    short	   fNumInt;
	/** @brief Location of the integration point Ksi */
    TPZVec<long double>	fLocationKsi;
	/** @brief Location of the integration point Eta */
    TPZVec<long double>	fLocationEta;
	/** @brief Weight of the integration point */
    TPZVec<long double>	fWeight;

	/**
	 * @brief Constructor of integration rule for triangle.
	 * @param order Order of the polinomial will be exactly integrated.
	 */
	TPZIntRuleT(int order);
	/** @brief Default destructor. */
	~TPZIntRuleT();

	/**
	 * @brief Computes the cubature rules following the symmetric construction presented at Linbo Zhang article.
	 * @param order Order of the polinomial will be exactly integrated.
	 */
	void ComputingSymmetricCubatureRule(int order);

	/**
	 * @brief Transforms barycentric coordinates (3 component) of the point in triange in cartesian coordinates (2 component).
	 * @param baryvec Vector of barycentric coordinates (3 components) for each point
	 * @param weightvec Vector of weights related for each point in baryvec vector
	 */
	void TransformBarycentricCoordInCartesianCoord(long double baryvec[],long double weightvec[]);

	/**
	 * @brief Checks sum of the weights is equal than measure of the master element, 
	 * and all of integration points belong to the master element.
	 * @return Returns false if one integration point is outside of the master element or the sum of weights is not one.
	 */
	bool CheckCubatureRule();

public:
    enum {NRULESTRIANGLE_ORDER = 21};

	/** @brief Returns number of integration points */
    short NInt() const { return fNumInt;}

	/** @brief Returns location of the ith point */
    void Loc(int i, TPZVec<REAL> &pos) const;

	/** @brief Return weight for the ith point */
    REAL W(int i) const;
};

#endif
