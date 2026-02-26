/**
 * @file
 * @brief Contains TPZSYsmpMatrixMumps class which implements a symmetric sparse matrix with MUMPS solver. \n
 * Purpose: Defines operations on symmetric sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 * @note Only real-valued types (double, float) are supported.
 * Attempting to use complex types will abort at runtime.
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
  
  void SetData(const TPZVec<int64_t> &IA, const TPZVec<int64_t> &JA, const TPZVec<TVar> &A) override;
  void SetData(TPZVec<int64_t> &&IA, TPZVec<int64_t> &&JA, TPZVec<TVar> &&A) override;
  
  //! Update COO format from CSR (called by StructMatrix during assembly)
  /**
    * @brief Converts the internal CSR (Compressed Sparse Row) representation
    *        to the COO (Coordinate) format required by MUMPS.
    *
    * NeoPZ stores sparse matrices in CSR format (fIA, fJA, fA), where fIA holds
    * row pointers, fJA holds column indices (0-based), and fA holds values.
    * MUMPS requires COO format with 1-based indices: explicit (row, col) pairs
    * for every non-zero entry.
    *
    * This method builds fIRN1Based (row indices) and fJCN1Based (column indices),
    * both 1-based, from the CSR structure. The values array fA is passed directly
    * to MUMPS as a pointer — no copy is made.
    *
    * The converted arrays are held only until MUMPS completes factorization and
    * are freed immediately after (see Decompose), since MUMPS uses its own
    * internal representation for the solve phase.
    *
    * @note fCOOValid tracks whether fIRN1Based/fJCN1Based are up-to-date.
    *       It is set to false by SetData (sparsity pattern changed) and after
    *       Decompose (arrays freed). UpdateCOOFormat is called lazily at the
    *       start of Decompose if fCOOValid is false.
    */
  void UpdateCOOFormat();
  
private:
  
  TPZMumpsSolver<TVar> fMumpsControl;
  
  // COO format arrays (1-based indexing for MUMPS)
  mutable TPZManVector<MUMPS_INT> fIRN1Based;
  mutable TPZManVector<MUMPS_INT  > fJCN1Based;
  mutable bool fCOOValid{false};
};

#endif
