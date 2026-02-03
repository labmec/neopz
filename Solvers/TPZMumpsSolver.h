#ifndef TPZMUMPSSOLVER_H
#define TPZMUMPSSOLVER_H

#include "TPZMatrixSolver.h"

#include "dmumps_c.h"
#include "pzmanvector.h"
#include "tpzautopointer.h"
#include <optional>

#define JOB_INIT -1
#define JOB_END -2
#define JOB_ANALYSIS 1
#define JOB_FACTORIZE 2
#define JOB_SOLVE 3

#define USE_COMM_WORLD -987654 // sequential run without MPI
#define MUMPS_HOST_PAR 1

template <class TVar>
class TPZFYsmpMatrixMumps;

template <class TVar>
class TPZSYsmpMatrixMumps;

template <typename TVar>
class TPZMumpsSolver : public TPZMatrixSolver<TVar> {

  // they need access to SetRawMatrix and ReallySolve
  friend class TPZFYsmpMatrixMumps<TVar>;
  friend class TPZSYsmpMatrixMumps<TVar>;

public:
  enum class MProperty { ENonInitialized = 0,
                         EPositiveDefinite,
                         EIndefinite };

  TPZMumpsSolver();

  // Delete copy operations - MUMPS data structures cannot be safely copied
  TPZMumpsSolver(const TPZMumpsSolver &copy) = delete;
  TPZMumpsSolver &operator=(const TPZMumpsSolver &copy) = delete;

  // Move operations are allowed
  TPZMumpsSolver(TPZMumpsSolver &&copy) noexcept;
  TPZMumpsSolver &operator=(TPZMumpsSolver &&copy) noexcept;

  virtual ~TPZMumpsSolver();

  void SetMatrix(TPZAutoPointer<TPZBaseMatrix> Refmat) override;

  void Decompose() override;

  void Solve(const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol, TPZFMatrix<TVar> *residual = nullptr) override;

  void SetMatrixType(SymProp symtype, MProperty prop);

  void SetMatrixType(long long type) { fMatrixType = type; }

  void SetMessageLevel(int lvl);

  int GetMessageLevel() const { return fMessageLevel; }

  TPZMumpsSolver<TVar> *Clone() const override;

  /**
   * Returns the ICNTL array used by MUMPS
   * @return Pointer to the ICNTL array
   */
  [[nodiscard]] inline MUMPS_INT *GetICNTL() { return fMumpsData.icntl; }

  /**
   * Sets a specific ICNTL parameter in MUMPS
   * @param index Index of the ICNTL parameter (1-based)
   * @param value Value to set for the ICNTL parameter
   */
  void SetICNTL(int index, MUMPS_INT value);

  void ResetICNTL();

  [[nodiscard]] bool HasCustomSettings() const { return fCustomSettings; }

  [[nodiscard]] inline DMUMPS_STRUC_C &GetMumpsData() { return fMumpsData; }

  [[nodiscard]] inline const DMUMPS_STRUC_C &GetMumpsData() const { return fMumpsData; }

  TPZVec<long long> &GetPermutationVec() { return fPermutation; }

  void FreeMumpsMemory();

protected:
  long long MatrixType();

  void Decompose(TPZMatrix<TVar> *mat);

  void Solve(const TPZMatrix<TVar> *mat, const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol) const;

  SymProp fSymmetry{SymProp::NonSym};

  MProperty fProperty{MProperty::ENonInitialized};

  mutable DMUMPS_STRUC_C fMumpsData;

  long long fMax_num_factors{1};

  long long fMatrix_num{1};

  long long fMessageLevel{0};

  long long fError{0};

  TPZVec<long long> fPermutation;

  long long fMatrixType{0};

  bool fDecomposed{false};

  bool fMumpsInitialized{false};

  bool fCustomSettings{false};

  std::optional<MUMPS_INT> fICNTL[60];
};

#endif /* TPZMUMPSSOLVER_H */
