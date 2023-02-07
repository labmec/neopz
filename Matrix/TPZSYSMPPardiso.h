/**
 * @file
 * @brief Contains TPZSYsmpMatrix class which implements a nonsymmetric sparse matrix. \n
 * Purpose: Defines operations on symmetric sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 * Some of the functionalities of this class depends on the MKL library and thus needs the NeoPZ library
 * to be configured using USING_MKL=ON during the CMake process. Search on this header for MKL to see which functionalities are these.
 */

#ifndef SYSMPMATPARDISO_H
#define SYSMPMATPARDISO_H

#include "TPZSYSMPMatrix.h"

#ifdef USING_MKL

#include "TPZPardisoSolver.h"

 /**
  * @brief Implements a symmetric sparse matrix using Pardiso. \ref matrix "Matrix"
  * @ingroup matrix
  */
template<class TVar>
class TPZSYsmpMatrixPardiso : public TPZSYsmpMatrix<TVar>{
	
    friend class TPZPardisoSolver<TVar>;
    
public :

  TPZSYsmpMatrixPardiso();
  /** @brief Constructors from parent class*/
  using TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix;
	/** @brief Copy constructor */
  TPZSYsmpMatrixPardiso(const TPZSYsmpMatrixPardiso<TVar> &cp) = default;
  /** @brief Move constructor*/
  TPZSYsmpMatrixPardiso(TPZSYsmpMatrixPardiso<TVar> &&cp) = default;
  /** @brief Copy-assignment operator*/
  TPZSYsmpMatrixPardiso &operator=(const TPZSYsmpMatrixPardiso<TVar> &copy) = default;
  /** @brief Move-assignment operator*/
  TPZSYsmpMatrixPardiso &operator=(TPZSYsmpMatrixPardiso<TVar> &&copy) = default;


  /** @brief Copy constructor from generic sparse matrix*/
  TPZSYsmpMatrixPardiso(const TPZSYsmpMatrix<TVar> &cp)
    : TPZSYsmpMatrix<TVar>(cp) {}
  /** @brief Move constructor from generic sparse matrix*/
  TPZSYsmpMatrixPardiso(TPZSYsmpMatrix<TVar> &&rval)
    : TPZSYsmpMatrix<TVar>(rval) {}
  /** @brief Copy-assignment operator from generic sparse matrix*/
  TPZSYsmpMatrixPardiso &operator=(const TPZSYsmpMatrix<TVar> &cp)
  { TPZSYsmpMatrix<TVar>::operator=(cp); return *this;}
  /** @brief Move-assignment operator from generic sparse matrix*/
  TPZSYsmpMatrixPardiso &operator=(TPZSYsmpMatrix<TVar> &&rval)
  { TPZSYsmpMatrix<TVar>::operator=(rval); return *this;}
  
  inline TPZSYsmpMatrixPardiso<TVar>*NewMatrix() const override {return new TPZSYsmpMatrixPardiso<TVar>{};}
  CLONEDEF(TPZSYsmpMatrixPardiso)
	/** @brief Destructor */
	~TPZSYsmpMatrixPardiso() = default;

  /** @brief Creates a copy from another sparse matrix*/
  void CopyFrom(const TPZMatrix<TVar> *  mat) override;

  void SetIsDecomposed(DecomposeType val) override;
  virtual int Decompose(const DecomposeType dt) override;
  virtual int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) override;
  virtual int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const override;

  int ClassId() const override;

  //! Gets reference to TPZPardisoSolver instance for fine-tuning
  TPZPardisoSolver<TVar> & GetPardisoControl()
  {return fPardisoControl;}
private:
  TPZPardisoSolver<TVar> fPardisoControl;
};

#endif
#endif
