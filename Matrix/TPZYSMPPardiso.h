/**
 * @file
 * @brief Contains TPZSYsmpMatrix class which implements a nonsymmetric sparse matrix. \n
 * Purpose: Defines operations on symmetric sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 * Some of the functionalities of this class depends on the MKL library and thus needs the NeoPZ library
 * to be configured using USING_MKL=ON during the CMake process. Search on this header for MKL to see which functionalities are these.
 */

#ifndef YSMPMATPARDISO_H
#define YSMPMATPARDISO_H

#include "TPZYSMPMatrix.h"

#ifdef USING_MKL

#include "TPZPardisoSolver.h"

 /**
  * @brief Implements a symmetric sparse matrix using Pardiso. \ref matrix "Matrix"
  * @ingroup matrix
  */
template<class TVar>
class TPZFYsmpMatrixPardiso : public TPZFYsmpMatrix<TVar>{
	
    friend class TPZPardisoSolver<TVar>;
    
public :

  TPZFYsmpMatrixPardiso();
  /** @brief Constructors from parent class*/
  using TPZFYsmpMatrix<TVar>::TPZFYsmpMatrix;
	/** @brief Copy constructor */
  TPZFYsmpMatrixPardiso(const TPZFYsmpMatrixPardiso<TVar> &cp) = default;
  /** @brief Move constructor*/
  TPZFYsmpMatrixPardiso(TPZFYsmpMatrixPardiso<TVar> &&cp) = default;
  /** @brief Copy-assignment operator*/
  TPZFYsmpMatrixPardiso &operator=(const TPZFYsmpMatrixPardiso<TVar> &copy) = default;
  /** @brief Move-assignment operator*/
  TPZFYsmpMatrixPardiso &operator=(TPZFYsmpMatrixPardiso<TVar> &&copy) = default;


  /** @brief Copy constructor from generic sparse matrix*/
  TPZFYsmpMatrixPardiso(const TPZFYsmpMatrix<TVar> &cp)
    : TPZFYsmpMatrix<TVar>(cp) {}
  /** @brief Move constructor from generic sparse matrix*/
  TPZFYsmpMatrixPardiso(TPZFYsmpMatrix<TVar> &&rval)
    : TPZFYsmpMatrix<TVar>(rval) {}
  /** @brief Copy-assignment operator from generic sparse matrix*/
  TPZFYsmpMatrixPardiso &operator=(const TPZFYsmpMatrix<TVar> &cp)
  { TPZFYsmpMatrix<TVar>::operator=(cp); return *this;}
  /** @brief Move-assignment operator from generic sparse matrix*/
  TPZFYsmpMatrixPardiso &operator=(TPZFYsmpMatrix<TVar> &&rval)
  { TPZFYsmpMatrix<TVar>::operator=(rval); return *this;}
  
  inline TPZFYsmpMatrixPardiso<TVar>*NewMatrix() const override {return new TPZFYsmpMatrixPardiso<TVar>{};}
  CLONEDEF(TPZFYsmpMatrixPardiso)
	/** @brief Destructor */
	~TPZFYsmpMatrixPardiso() = default;

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
