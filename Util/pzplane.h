/**
 * @file
 * @brief Contains declaration of the TPZPlane class which implements a plane.
 */

#ifndef TPZPLANE_H
#define TPZPLANE_H

#include <cmath>
#include <iostream>
#include "pzvec.h"
#include "pznumeric.h"

/**
 * @ingroup util
 * @brief Implements a plane (stores the plane equation). \ref util "Utility"
 * @author Renato Gomes Damas
 */
class TPZPlane {
	
public:
	/** @brief Simple constructor */
	TPZPlane();
	/** @brief Destructor */
	~TPZPlane();
	
	/** @brief Computes the equation of the plane from three given points */
	int SetPlane(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2, const TPZVec<REAL> &p3);

	/** @brief Check whether the point belongs at the plane */
	bool Belongs(const TPZVec<REAL> &point);

	/** @brief Checks whether the current plane coincides with the plane formed by three given points */
	bool Belongs(const TPZVec<REAL> &ponto1, const TPZVec<REAL> &ponto2, const TPZVec<REAL> &ponto3);
	
private:
	
	/** @brief Coefficients of the plane with equation: fCoef[0]*x + fCoef[1]*y + fCoef[2]*z + fCoef[3] = 0 */
	TPZVec<REAL> fCoef;

};

#endif
