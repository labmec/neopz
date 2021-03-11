/*
 *  TPZStokesMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZMatWithMem.h"

#include "pzfmatrix.h"
#include "pzbndcond.h"
#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZMaterial.h"
#include "TPZStokesMaterial.h"

#ifndef TPZHSTOKESMATERIAL
#define TPZHSTOKESMATERIAL

class TPZHStokesMaterial: public TPZStokesMaterial
{
  
   REAL fBeta = 0.6;
public:
      /**
     * Empty Constructor
     */
    TPZHStokesMaterial() : TPZStokesMaterial()
    {
    }
    
    /** Creates a material object and inserts it in the vector of
     *  material pointers of the mesh.
     */
    TPZHStokesMaterial(int matid, int dimension, int space, STATE viscosity, STATE theta, STATE Sigma) : TPZStokesMaterial(matid,dimension,space,viscosity,theta,Sigma)
      {
      }
    
    
    /** Creates a material object based on the referred object and
     *  inserts it in the vector of material pointers of the mesh.
     */
    TPZHStokesMaterial(const TPZStokesMaterial &mat) : TPZStokesMaterial(mat)
    {}
    
    /**
     * Destructor
     */
    ~TPZHStokesMaterial()
    {
    }
    
    virtual TPZMaterial *NewMaterial()
    {
      return new TPZHStokesMaterial(*this);
    }
    
        /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

};


#endif
