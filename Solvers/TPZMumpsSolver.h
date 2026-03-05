#ifndef TPZMUMPSSOLVER_H
#define TPZMUMPSSOLVER_H

#include "TPZMatrixSolver.h"

#ifdef MUMPS_HAVE_SINGLE
#include "smumps_c.h"
#endif
#ifdef MUMPS_HAVE_DOUBLE
#include "dmumps_c.h"
#endif
#ifdef MUMPS_HAVE_COMPLEX
#include "cmumps_c.h"
#endif
#ifdef MUMPS_HAVE_COMPLEX16
#include "zmumps_c.h"
#endif

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

/**
 * @brief Maps TVar to the corresponding MUMPS struct type.
 *
 * Specializations are enabled only for the variants built into MUMPS
 * (controlled by MUMPS_HAVE_SINGLE/DOUBLE/COMPLEX/COMPLEX16 defines).
 * Using an unsupported TVar will produce a clear "incomplete type" error.
 */
template<class TVar> struct MumpsStrucSelector {
       static_assert(sizeof(TVar)==0,
           "No MUMPS support for this TVar. Check MUMPS_HAVE_* compile definitions.");
   };
   #ifdef MUMPS_HAVE_SINGLE
   template<> struct MumpsStrucSelector<float>                { using type = SMUMPS_STRUC_C; };
   #endif
   #ifdef MUMPS_HAVE_DOUBLE
   template<> struct MumpsStrucSelector<double>               { using type = DMUMPS_STRUC_C; };
   template<> struct MumpsStrucSelector<long double>          { using type = DMUMPS_STRUC_C; };
   #endif
   #ifdef MUMPS_HAVE_COMPLEX
   template<> struct MumpsStrucSelector<std::complex<float>>  { using type = CMUMPS_STRUC_C; };
   #endif
   #ifdef MUMPS_HAVE_COMPLEX16
   template<> struct MumpsStrucSelector<std::complex<double>> { using type = ZMUMPS_STRUC_C; };
   template<> struct MumpsStrucSelector<std::complex<long double>> { using type = ZMUMPS_STRUC_C; };
   #endif
   template<typename TVar> using MumpsStruc_t = typename MumpsStrucSelector<TVar>::type;

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

  [[nodiscard]] inline MumpsStruc_t<TVar> &GetMumpsData() { return fMumpsData; }

  [[nodiscard]] inline const MumpsStruc_t<TVar> &GetMumpsData() const { return fMumpsData; }

  TPZVec<long long> &GetPermutationVec() { return fPermutation; }

  void FreeMumpsMemory();

protected:
  long long MatrixType();

  void Decompose(TPZMatrix<TVar> *mat);

  void Solve(const TPZMatrix<TVar> *mat, const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol) const;

  void CallMumps() const;

  SymProp fSymmetry{SymProp::NonSym};

  MProperty fProperty{MProperty::ENonInitialized};

  mutable MumpsStruc_t<TVar> fMumpsData;

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
