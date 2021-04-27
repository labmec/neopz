/**
 * @file
 * @brief Contains the TPBSpStructMatrix class which assembles on the pair equations.
 */

#ifndef TPBSPSTRUCTMATRIX_H
#define TPBSPSTRUCTMATRIX_H

#include "TPZSpStructMatrix.h"

#include "pzstrmatrixor.h"
/** 
 * @ingroup structural
 * @brief Assembles only the pair equations.
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZBSpStructMatrix : public TPZSpStructMatrix<TVar,TPar> {
public:    
    using TPZSpStructMatrix<TVar,TPar>::TPZSpStructMatrix;
    
    int ClassId() const override;

    TPZMatrix<TVar> * Create() override;
	
    TPZStructMatrix * Clone() override;
};

#endif //TPBSPSTRUCTMATRIX_H
