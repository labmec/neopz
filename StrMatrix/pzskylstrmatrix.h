/**
 * @file
 * @brief Contains the TPZSkylineStructMatrix class which implements SkyLine Structural Matrices.
 */

#ifndef TPZSKYLINESTRUCTMATRIX_H
#define TPZSKYLINESTRUCTMATRIX_H

#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"

class TPZCompMesh;

/**
 * @brief Implements a skyline structural matrix using TPZSkylMatrix as a storage format.
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSkylineStructMatrix : public TPZStructMatrixT<TVar>,
                               public TPar{
protected:
    
    /** @brief the equations which should actually be assembled */
    TPZVec<int64_t> fActiveEquations;
    
    /** @brief Equation destination */
    TPZVec<int64_t> fEquationDestination;
    
    /** Returns the skyline matrix object */
    virtual TPZMatrix<TVar> * ReallyCreate(int64_t neq, const TPZVec<int64_t> &skyline);
    
public:    
    using TPZStructMatrixT<TVar>::TPZStructMatrixT;
    
    TPZMatrix<TVar> * Create() override;
	
    TPZStructMatrix * Clone() override;
    //@{
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
private :
    friend TPZPersistenceManager;
};

#endif //TPZSKYLINESTRUCTMATRIX_H
