/**
 * \file
 * @brief Contains the TPZViscoelastic class which implements an isotropic viscoelasticity material.
 */
#ifndef TPZVISCOELASTIC_H
#define TPZVISCOELASTIC_H

/** @brief First index to qsi vector for contribute method of viscoelasticity material */
const int _XX_ = 0;
/** @brief Second index to qsi vector for contribute method of viscoelasticity material */
const int _XY_ = 1;
/** @brief Third index to qsi vector for contribute method of viscoelasticity material */
const int _XZ_ = 2;
/** @brief Fourth index to qsi vector for contribute method of viscoelasticity material */
const int _YY_ = 3;
/** @brief Fifth index to qsi vector for contribute method of viscoelasticity material */
const int _YZ_ = 4;
/** @brief Sixth index to qsi vector for contribute method of viscoelasticity material */
const int _ZZ_ = 5;


#include <iostream>
#include "pzfmatrix.h"
#include "pzvec.h"
#include <vector>

#include "pzelast3d.h"
#include "pzmatwithmem.h"

/**
 * @ingroup material
 * @brief This class implements an isotropic viscoelasticity material.
 * @author Pamela Diaz
 * @since 7/16/2010.
 */
class TPZViscoelastic : public TPZMatWithMem<TPZFMatrix, TPZElasticity3D>
{
	
public:
	TPZViscoelastic( TPZMatWithMem<TPZFMatrix, TPZElasticity3D> &matwithmem, int id,REAL lambdaE,REAL muE, REAL lambdaV, REAL muV, REAL alphaT);
	
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix &ek,
							TPZFMatrix &ef);
	
	
protected:
	
	REAL flambdaE,fmuE,flambdaV,fmuV,falphaT;  
	//int fMemory;
	
};
#endif