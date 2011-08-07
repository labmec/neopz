/**
 * @file
 * @brief Contains declaration of the TPZPlane class which implements a plane.
 */
/***************************************************************************
                          TPZPlane.h  -  description
                             -------------------
    begin                : Wed Mar 27 2002
    copyright            : (C) 2002 by Renato Gomes Damas
    email                : rgdamas@fec.unicamp.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

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
	/** Dado tres pontos calcula a equacao do plano que os contem. */
	/** @brief Computes the equation of the plane from three given points */
	int SetPlane(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2, const TPZVec<REAL> &p3);
	/** Verifica se o point[3] pertence ao plano. Se pertencer retorna 1, caso contrario 0. */
	/** @brief Check whether point belongs at the plane */
	bool Belongs(const TPZVec<REAL> &point);
	/** Verifica se o plano coincide com plano formado pelos tres pontos passados. Se pertencer retorna 1, caso contrario 0. */
	/** @brief Checks whether the current plane coincides with the plane formed by three given points */
	bool Belongs(const TPZVec<REAL> &ponto1, const TPZVec<REAL> &ponto2, const TPZVec<REAL> &ponto3);
	
private:
	
	/** Coeficientes da equação do plano: fCoef[0]*x + fCoef[1]*y + fCoef[2]*z + fCoef[3] = 0. */
	/** @brief Coefficients of the plane with equation: fCoef[0]*x + fCoef[1]*y + fCoef[2]*z + fCoef[3] = 0 */
	TPZVec<REAL> fCoef;
	
};

#endif
