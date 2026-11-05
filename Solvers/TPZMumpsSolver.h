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
 * Specializations exist for variants built into MUMPS
 * (controlled by MUMPS_HAVE_SINGLE/DOUBLE/COMPLEX/COMPLEX16 defines).
 * The primary template uses `char` as a placeholder so that
 * TPZMumpsSolver<T> and TPZFYsmpMatrixMumps<T>/TPZSYsmpMatrixMumps<T>
 * are well-formed for any T. Methods that call actual MUMPS routines are
 * guarded by MumpsTypeIsSupported_v<T> and will DebugStop() if called
 * for an unsupported type at runtime.
 */
template<class TVar> struct MumpsStrucSelector { using type = char; };
   #ifdef MUMPS_HAVE_SINGLE
   template<> struct MumpsStrucSelector<float>                { using type = SMUMPS_STRUC_C; };
   #endif
   #ifdef MUMPS_HAVE_DOUBLE
   template<> struct MumpsStrucSelector<double>               { using type = DMUMPS_STRUC_C; };
   #endif
   #ifdef MUMPS_HAVE_COMPLEX
   template<> struct MumpsStrucSelector<std::complex<float>>  { using type = CMUMPS_STRUC_C; };
   #endif
   #ifdef MUMPS_HAVE_COMPLEX16
   template<> struct MumpsStrucSelector<std::complex<double>> { using type = ZMUMPS_STRUC_C; };
   #endif
   template<typename TVar> using MumpsStruc_t = typename MumpsStrucSelector<TVar>::type;

/**
 * @brief True when TVar has an available MUMPS variant (i.e., MUMPS_HAVE_* is set).
 * Used to guard methods that access the MUMPS struct fields.
 */
template<class TVar> struct MumpsTypeIsSupported : std::false_type {};
#ifdef MUMPS_HAVE_SINGLE
template<> struct MumpsTypeIsSupported<float>                      : std::true_type {};
#endif
#ifdef MUMPS_HAVE_DOUBLE
template<> struct MumpsTypeIsSupported<double>                     : std::true_type {};
#endif
#ifdef MUMPS_HAVE_COMPLEX
template<> struct MumpsTypeIsSupported<std::complex<float>>        : std::true_type {};
#endif
#ifdef MUMPS_HAVE_COMPLEX16
template<> struct MumpsTypeIsSupported<std::complex<double>>       : std::true_type {};
#endif
template<class TVar>
inline constexpr bool MumpsTypeIsSupported_v = MumpsTypeIsSupported<TVar>::value;

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
  [[nodiscard]] inline MUMPS_INT *GetICNTL() {
    if constexpr (MumpsTypeIsSupported_v<TVar>) {
      return fMumpsData.icntl;
    } else {
      DebugStop();
      return nullptr;
    }
  }

  /**
   * Sets a specific ICNTL parameter in MUMPS
   * @param index Index of the ICNTL parameter (1-based)
   * @param value Value to set for the ICNTL parameter
   */
  void SetICNTL(int index, MUMPS_INT value);

  void ResetICNTL();

  [[nodiscard]] bool HasCustomSettings() const { return fCustomSettings; }

  /**
   * Locks the current settings so that TPZSYsmpMatrixMumps::Decompose will
   * not auto-detect the matrix type from the matrix's SymProp.
   * Call this after SetMatrixType() when you need to override the auto-detection
   * (e.g., to use SymProp::Sym for a complex matrix that reports SymProp::Herm).
   */
  void LockSettings() { fCustomSettings = true; }

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
