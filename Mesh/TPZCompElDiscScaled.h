//
//  TPZCompElDiscScaled.h
//  MixedElasticityGirkmann
//
//  Created by Philippe Devloo on 03/05/18.
//

#ifndef TPZCompElDiscScaled_h
#define TPZCompElDiscScaled_h

#include <stdio.h>
#include <iostream>
#include "TPZCompElDisc.h"
#include "pzgeoel.h"


class TPZCompElDiscScaled : public TPZCompElDisc {

protected:
    
    REAL fScale = 1.;
   
    /** @brief Compute the shape functions scale related to reference characteristic size*/
    void ComputeScale();
    
public:
    
    /** @brief Default constructor */

    TPZCompElDiscScaled() : TPZCompElDisc() {
        ComputeScale();
    }
    
    /** @brief Constructor of the discontinuous element associated with geometric element */
    TPZCompElDiscScaled(TPZCompMesh &mesh,TPZGeoEl *ref);

    /** @brief Copy constructor */
    TPZCompElDiscScaled(TPZCompMesh &mesh, const TPZCompElDiscScaled &copy) : TPZCompElDisc(mesh,copy), fScale(copy.fScale){
    }
    
    /** @brief Default destructor */
    virtual ~TPZCompElDiscScaled();
    
    /** @brief Return the shape functions scale */
    REAL Scale() const {
        return fScale;
    }

    /** @brief Set a constant value to shape functions scale */
    void SetScale(REAL c){
        fScale = c;
    }
    
    /** @brief Compute shape functions multiplied by scale value*/
    void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                                     TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                                     REAL &detjac, TPZFMatrix<REAL> &jacinv,
                              TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphix) override;    

    
    TPZCompEl *Clone(TPZCompMesh &mesh) const override {
            return new TPZCompElDiscScaled(mesh,*this);
    }
    
};

#endif /* TPZCompelDiscScaled_h */

