/**
 ************************************************************************** TPZLine.h  -  description -------------------
 * begin                : Tue Apr 16 2002 copyright            : (C) 2002 by Renato Gomes Damas
 * email                : rgdamas@fec.unicamp.br
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
 * @ingroup util
 * @brief Implements a line.
 */
class TPZLine {
public:
	TPZLine();
	~TPZLine();

	/** Armazena um ponto da reta e sua direcao. */
	/** @brief Store a point and direction vector of the line */
	void SetLine(const TPZVec<REAL> &point, const TPZVec<REAL> &dir);

	/** Verifica se o ponto pertence a reta. */
	/** @brief Verify whether the point belongs at the line */
	bool Belongs(const TPZVec<REAL> &point);

	/** Fornece a tolerancia do calculo, armazenando em tol. */
	/** @brief Returns the tolerance value considered to compute. */
	REAL GetTolerance();

	/** Especifica a tolerancia para os calculos */
	/** @brief Sets the tolerance value to compute. */
	void SetTolerance(const REAL &tol);

private:
	/** um ponto qualquer da reta. */
	/** @brief Any point belongs at the line */
	TPZVec<REAL> fPoint;

	/** Direcao da reta. */
	/** @brief Vector direction of the line */
	TPZVec<REAL> fDirection;

	/** Tolerancia para os calculos */
	/** @brief Tolerance value to computes. */
	REAL fTolerance;
};

#endif
