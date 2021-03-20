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


#include "pzmanvector.h"
#include "tpzautopointer.h"
#include "pzfmatrix.h"
template<class TVar>
class TPZSYsmpMatrix;

template<class TVar>
class TPZFYsmpMatrix;

/// class to control the pardiso solution process
// inspired by PardisoSolver by Armando Duarte
// it has no use unless USING_MKL=ON is set during
// CMake configuration of the NeoPZ library.
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
    // 32-bit: int pt[64]; 64-bit: long long pt[64]
    // or void *pt[64] should be OK on both architectures
    // this datastructure should not be copied or duplicated, therefore the "autopointer" protection
    
    TPZAutoPointer<TPZManVector<long long, 64> > fPardisoControl;
    
    // adress of the first element of pt;
    long long *fHandle;
    
    // Array used to pass parameters to Pardiso
    TPZManVector<long long, 64> fParam;
    
    // Maximum number of factors we will pass to the solver
    long long fMax_num_factors;
    
    // Factor number we are using
    long long fMatrix_num;
    
    // Message level information
    long long fMessageLevel;
    
    // error flag from Pardiso
    long long fError;
    
    /// permutation vector computed by Pardiso
    TPZVec<long long> fPermutation;
    
    // matrix type, computed based on the structural information and TVar
    long long fMatrixType;
    
    /// Compute the matrix type
    long long MatrixType();
    
    /// pointer to the nonsymmetric system (where the data structures are stored
    TPZFYsmpMatrix<TVar> *fNonSymmetricSystem;
    
    /// pointer to the symmetric system (where the data structures are stored
    TPZSYsmpMatrix<TVar> *fSymmetricSystem;
    
    /// Provides a explanation for the given error
    void Error_check(int error) const;
    
};

#endif /* TPZPardisoControl_hpp */
