/**
 * @file
 * @brief Contains declaration of the TPZLine class which implements a line.
 */

/**
 **************************************************************************
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 */
#ifndef TPZLINE_H
#define TPZLINE_H
#include "pzvec.h"
#include "pzreal.h"
#include <cmath>

/** 
 * @author Renato Gomes Damas 
 * @since Apr 16, 2002.
 * @ingroup util
 * @brief Implements a line. \ref util "Utility"
 */
class TPZLine {
public:
	TPZLine();
	~TPZLine();

	/** @brief Store a point and direction vector of the line */
	void SetLine(const TPZVec<REAL> &point, const TPZVec<REAL> &dir);

	/** @brief Verify whether the point belongs at the line */
	bool Belongs(const TPZVec<REAL> &point);

	/** @brief Returns the tolerance value considered to compute. */
	REAL GetTolerance();

	/** @brief Sets the tolerance value to compute. */
	void SetTolerance(const REAL &tol);

private:
	/** @brief Any point belongs at the line */
	TPZVec<REAL> fPoint;

	/** @brief Vector direction of the line */
	TPZVec<REAL> fDirection;

	/** @brief Tolerance value to computes. */
	REAL fTolerance;
};

#endif
