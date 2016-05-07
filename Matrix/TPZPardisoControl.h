//
//  TPZPardisoControl.hpp
//  PZ
//
//  Created by Philippe Devloo on 5/5/16.
//
//

#ifndef TPZPardisoControl_hpp
#define TPZPardisoControl_hpp
#ifdef USING_MKL

#include <stdio.h>

#include "mkl_pardiso.h"
#include "pzmanvector.h"


/// class to control the pardiso solution process
// inspired by PardisoSolver by Armando Duarte
template<class TVar>
class TPZPardisoControl
{
public:
    enum MSystemType {ESymmetric, EHermitian, EnonSymmetric};
    
    enum MStructure {EStructureSymmetric, EStructureNonSymmetric};
    
    enum MProperty {EPositiveDefinite, EIndefinite};
    
    
protected:
    MSystemType fSystemType;
    
    MStructure fStructure;
    
    MProperty fProperty;
    
    // Solver internal data address pointers
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    TPZManVector<void*, 64> fPardisoControl;
    
    // adress of the first element of pt;
    _MKL_DSS_HANDLE_t fHandle;
    //  ConcreteRigidArray1d<long int, 64> pt;
    
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
    
public:
    
    TPZPardisoControl(MSystemType systemtype, MProperty prop);
    
    
    
};

#endif
#endif /* TPZPardisoControl_hpp */
