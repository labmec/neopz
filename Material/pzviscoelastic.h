/**
 * \file
 * @brief Contains the TPZViscoelastic class which implements an isotropic viscoelasticity material.
 */
#ifndef TPZVISCOELASTIC_H
#define TPZVISCOELASTIC_H

const int _XX_ = 0;
const int _XY_ = 1;
const int _XZ_ = 2;
const int _YY_ = 3;
const int _YZ_ = 4;
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