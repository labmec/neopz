/*
 *  pzmaterialcoupling.h
 *  PZ
 *
 *  Created by Denise de Siqueira on 8/1/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PZMATERIALCOUPLING_H
#define PZMATERIALCOUPLING_H

#include "pzpoisson3d.h"
#include "pzmaterialdata.h"
#include "pzfmatrix.h"

/**
 *Implemented a Poisson Problem coupling the interpolation spaces H(div) and H1
 */
class TPZMaterialCoupling : public TPZMatPoisson3d {
		
public:
		
		
		
		TPZMaterialCoupling();
		TPZMaterialCoupling(int nummat, int dim);
		
		//TPZMaterialCoupling(const TPZMaterialCoupling &copy);
		/**
		 * method to possibilite the coupling between H(div) and H1  
		 */
		virtual void ContributeInterface(TPZMaterialData &data,TPZMaterialData &dataleft,TPZMaterialData &dataright, 
                                         REAL weight,TPZFMatrix &ek,TPZFMatrix &ef);
		virtual void ContributeInterface2(TPZMaterialData &data, TPZMaterialData &dataleft,TPZMaterialData &dataright, REAL weight,TPZFMatrix &ek,TPZFMatrix &ef);
    virtual void InitMaterialData(TPZMaterialData &data);		
		
		virtual ~TPZMaterialCoupling();
		
		
};



#endif