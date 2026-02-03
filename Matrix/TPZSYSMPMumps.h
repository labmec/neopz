/**
 * @file
 * @brief Contains TPZSYsmpMatrix class which implements a nonsymmetric sparse matrix. \n
 * Purpose: Defines operations on symmetric sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 * Some of the functionalities of this class depends on the MKL library and thus needs the NeoPZ library
 * to be configured using USING_MKL=ON during the CMake process. Search on this header for MKL to see which functionalities are these.
 */

#ifndef SYSMPMATMUMPS_H
#define SYSMPMATMUMPS_H

#include "TPZSYSMPMatrix.h"
#include "TPZMumpsSolver.h"

 /**
  * @brief Implements a symmetric sparse matrix using Mumps. \ref matrix "Matrix"
  * @ingroup matrix
  */
template<class TVar>
class TPZSYsmpMatrixMumps : public TPZSYsmpMatrix<TVar>{
	
    friend class TPZMumpsSolver<TVar>;
    
public :

  TPZSYsmpMatrixMumps();
  /** @brief Constructors from parent class*/
  using TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix;
	
  /** @brief Copy constructor - creates new MUMPS instance */
  TPZSYsmpMatrixMumps(const TPZSYsmpMatrixMumps<TVar> &cp);
  
  /** @brief Move constructor*/
  TPZSYsmpMatrixMumps(TPZSYsmpMatrixMumps<TVar> &&cp);
  
  /** @brief Copy-assignment operator - creates new MUMPS instance*/
  TPZSYsmpMatrixMumps &operator=(const TPZSYsmpMatrixMumps<TVar> &copy);
  
  /** @brief Move-assignment operator*/
  TPZSYsmpMatrixMumps &operator=(TPZSYsmpMatrixMumps<TVar> &&copy);


  /** @brief Copy constructor from generic sparse matrix*/
  TPZSYsmpMatrixMumps(const TPZSYsmpMatrix<TVar> &cp)
    : TPZSYsmpMatrix<TVar>(cp) {}
  /** @brief Move constructor from generic sparse matrix*/
  TPZSYsmpMatrixMumps(TPZSYsmpMatrix<TVar> &&rval)
    : TPZSYsmpMatrix<TVar>(rval) {}
  /** @brief Copy-assignment operator from generic sparse matrix*/
  TPZSYsmpMatrixMumps &operator=(const TPZSYsmpMatrix<TVar> &cp)
  { TPZSYsmpMatrix<TVar>::operator=(cp); return *this;}
  /** @brief Move-assignment operator from generic sparse matrix*/
  TPZSYsmpMatrixMumps &operator=(TPZSYsmpMatrix<TVar> &&rval)
  { TPZSYsmpMatrix<TVar>::operator=(rval); return *this;}
  
  inline TPZSYsmpMatrixMumps<TVar>*NewMatrix() const override {return new TPZSYsmpMatrixMumps<TVar>{};}
  CLONEDEF(TPZSYsmpMatrixMumps)
	/** @brief Destructor */
	~TPZSYsmpMatrixMumps() = default;

  /** @brief Creates a copy from another sparse matrix*/
  void CopyFrom(const TPZMatrix<TVar> *  mat) override;

  void SetIsDecomposed(DecomposeType val) override;
  virtual int Decompose(const DecomposeType dt) override;
  virtual int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) override;
  virtual int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const override;

   void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const override;
  
  int ClassId() const override;

  //! Gets reference to TPZMumpsSolver instance for fine-tuning
  TPZMumpsSolver<TVar> & GetMumpsControl()
  {return fMumpsControl;}
  
  //! Update COO format from CSR (called by StructMatrix during assembly)
  void UpdateCOOFormat();
  
  //! Get COO format arrays (1-based indexing for MUMPS)
  void GetCOOFormat(TPZVec<MUMPS_INT> &irn, TPZVec<MUMPS_INT> &jcn) const;
  
private:
  
  TPZMumpsSolver<TVar> fMumpsControl;
  
  // COO format arrays (1-based indexing for MUMPS)
  mutable TPZManVector<MUMPS_INT> fIRN1Based;
  mutable TPZManVector<MUMPS_INT  > fJCN1Based;
  mutable bool fCOOValid{false};
};

#endif
