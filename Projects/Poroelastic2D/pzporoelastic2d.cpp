/*
 *  pzporoelastic2d.cpp
 *  PZ
 *
 *  Created by Agnaldo on 11/28/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "pzporoelastic2d.h"

TPZPoroElastic2d::TPZPoroElastic2d():TPZDiscontinuousGalerkin(), fDim(1){
	
}

TPZPoroElastic2d::TPZPoroElastic2d(int matid, int dim):TPZDiscontinuousGalerkin(matid),fDim(dim){
	
}

TPZPoroElastic2d::~TPZPoroElastic2d(){
}
