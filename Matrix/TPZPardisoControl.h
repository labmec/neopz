//
//  TPZPardisoControl.hpp
//  PZ
//
//  Created by Philippe Devloo on 5/5/16.
//
//

#ifndef TPZPardisoControl_hpp
#define TPZPardisoControl_hpp

#include "pz_config.h"

#ifdef USING_MKL

#include <stdio.h>

//#include "mkl_pardiso.h"
#include "pzmanvector.h"
#include "tpzautopointer.h"
#include "pzfmatrix.h"
template<class TVar>
class TPZSYsmpMatrix;

template<class TVar>
class TPZFYsmpMatrix;

/// class to control the pardiso solution process
// inspired by PardisoSolver by Armando Duarte
template<class TVar>
class TPZPardisoControl
{
public:
    enum MSystemType {ESymmetric, EHermitian, ENonSymmetric};
    
    enum MStructure {EStructureSymmetric, EStructureNonSymmetric};
    
    enum MProperty {EPositiveDefinite, EIndefinite};
    
    /// empty constructor (non symetric and LU decomposition
    TPZPardisoControl();
    
public:
    /// constructor indicating if the matrix is symmetric and what type of decomposition is expected
    TPZPardisoControl(MSystemType systemtype, MProperty prop);
    
    TPZPardisoControl(const TPZPardisoControl &copy);
    
    TPZPardisoControl &operator=(const TPZPardisoControl &copy);
    
    virtual ~TPZPardisoControl();
    
    /// change the matrix type
    // this method should only be called if the pardiso control is zero (non initialized)
    void SetMatrixType(MSystemType systemtype, MProperty prop);
    
    /// initialize the pointer to the nonsymmetric data structure
    void SetMatrix(TPZFYsmpMatrix<TVar> *matrix)
    {
        fNonSymmetricSystem = matrix;
    }
    
    /// initialize the pointer to the symmetric data structure
    void SetMatrix(TPZSYsmpMatrix<TVar> *matrix)
    {
        fSymmetricSystem = matrix;
    }
    
    /// decompose the matrix acording to the method determined by SetMatrixType
    void Decompose();
    
    /// Use the decomposed matrix to invert the system of equations
    void Solve(TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol) const;
    
protected:
    
    MSystemType fSystemType;
    
    MStructure fStructure;
    
    MProperty fProperty;
    
    // Solver internal data address pointers
    // 32-bit: int pt[64]; 64-bit: int64_t pt[64]
    // or void *pt[64] should be OK on both architectures
    // this datastructure should not be copied or duplicated, therefore the "autopointer" protection
    
    TPZAutoPointer<TPZManVector<int64_t, 64> > fPardisoControl;
    
    // adress of the first element of pt;
    int64_t *fHandle;
    //  ConcreteRigidArray1d<int64_t, 64> pt;
    
    // Array used to pass parameters to Pardiso
    TPZManVector<int64_t, 64> fParam;
    
    // Maximum number of factors we will pass to the solver
    int64_t fMax_num_factors;
    
    // Factor number we are using
    int64_t fMatrix_num;
    
    // Message level information
    int64_t fMessageLevel;
    
    // error flag from Pardiso
    int64_t fError;
    
    /// permutation vector computed by Pardiso
    TPZVec<int64_t> fPermutation;
    
    // matrix type, computed based on the structural information and TVar
    int64_t fMatrixType;
    
    /// Compute the matrix type
    int64_t MatrixType();
    
    /// pointer to the nonsymmetric system (where the data structures are stored
    TPZFYsmpMatrix<TVar> *fNonSymmetricSystem;
    
    /// pointer to the symmetric system (where the data structures are stored
    TPZSYsmpMatrix<TVar> *fSymmetricSystem;
    
};

#endif
#endif /* TPZPardisoControl_hpp */
