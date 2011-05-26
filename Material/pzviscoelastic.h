/*
 *  pzviscoelastic.h
 *  pos_processamento
 *
 *  Created by Pamela Diaz on 7/16/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
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

/** This class implements an isotropic viscoelasticity material.
 *  @since Aug 31, 2005.
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