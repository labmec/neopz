//
//  TPZMumpsSolver.cpp
//  PZ
//
//
//
//

#ifdef USING_MUMPS
#include "TPZMumpsSolver.h"
#else
#define NOMUMPS
PZError << "The class TPZMumpsSolver should not be used ";
PZError << "if NeoPZ was configured with USING_MUMPS=OFF\n";
PZError << "Aborting..." << std::endl;
DebugStop();
#endif

#include "TPZSYSMPMumps.h"
#include "TPZYSMPMumps.h"
#include "pzlog.h"

/*auxiliary functions*/
template <typename TVar>
void Error_check(MUMPS_INT error, TVar info2);

template <class TVar>
int DataType([[maybe_unused]] TVar a);

#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.mumpscontrol");
#endif

/// empty constructor (non symetric and LU decomposition
template <class TVar>
TPZMumpsSolver<TVar>::TPZMumpsSolver() : TPZMatrixSolver<TVar>() {}

template <class TVar>
TPZMumpsSolver<TVar>::TPZMumpsSolver(TPZMumpsSolver &&copy) noexcept
    : TPZMatrixSolver<TVar>(std::move(copy)),
      fSymmetry(copy.fSymmetry),
      fProperty(copy.fProperty),
      fMumpsData(copy.fMumpsData),
      fMax_num_factors(copy.fMax_num_factors),
      fMatrix_num(copy.fMatrix_num),
      fMessageLevel(copy.fMessageLevel),
      fError(copy.fError),
      fPermutation(std::move(copy.fPermutation)),
      fMatrixType(copy.fMatrixType),
      fDecomposed(copy.fDecomposed),
      fMumpsInitialized(copy.fMumpsInitialized),
      fCustomSettings(copy.fCustomSettings) {
  // Mark source as uninitialized so it won't try to free MUMPS memory
  copy.fMumpsInitialized = false;
  copy.fDecomposed = false;
}

template <class TVar>
void TPZMumpsSolver<TVar>::CallMumps() const {
#ifdef MUMPS_HAVE_SINGLE
  if constexpr (std::is_same_v<TVar, float>) {
    smumps_c(&fMumpsData);
    return;
  }
#endif
#ifdef MUMPS_HAVE_DOUBLE
  if constexpr (std::is_same_v<TVar, double> || std::is_same_v<TVar, long double>) {
    dmumps_c(&fMumpsData);
    return;
  }
#endif
#ifdef MUMPS_HAVE_COMPLEX
  if constexpr (std::is_same_v<TVar, std::complex<float>>) {
    cmumps_c(&fMumpsData);
    return;
  }
#endif
#ifdef MUMPS_HAVE_COMPLEX16
  if constexpr (std::is_same_v<TVar, std::complex<double>> || std::is_same_v<TVar, std::complex<long double>>) {
    zmumps_c(&fMumpsData);
    return;
  }
#endif
  // If we reach here, no MUMPS variant matched TVar. This means the
  // MumpsStrucSelector specialisation exists but CallMumps() was not updated
  // to handle the new type, or the matching MUMPS_HAVE_* macro is missing.
  PZError << __PRETTY_FUNCTION__
          << ": no MUMPS variant found for the current TVar."
             " Check MUMPS_HAVE_* compile definitions.\n";
  DebugStop();
}

template <class TVar>
TPZMumpsSolver<TVar> &TPZMumpsSolver<TVar>::operator=(TPZMumpsSolver &&copy) noexcept {
  if (this != &copy) {
    // Free our own MUMPS memory first
    FreeMumpsMemory();

    // Move base class
    TPZMatrixSolver<TVar>::operator=(std::move(copy));

    // Move members
    fSymmetry = copy.fSymmetry;
    fProperty = copy.fProperty;
    fMumpsData = copy.fMumpsData;
    fMax_num_factors = copy.fMax_num_factors;
    fMatrix_num = copy.fMatrix_num;
    fMessageLevel = copy.fMessageLevel;
    fError = copy.fError;
    fPermutation = std::move(copy.fPermutation);
    fMatrixType = copy.fMatrixType;
    fDecomposed = copy.fDecomposed;
    fMumpsInitialized = copy.fMumpsInitialized;
    fCustomSettings = copy.fCustomSettings;

    // Mark source as uninitialized so it won't try to free MUMPS memory
    copy.fMumpsInitialized = false;
    copy.fDecomposed = false;
  }
  return *this;
}

template <class TVar>
void TPZMumpsSolver<TVar>::FreeMumpsMemory() {
#ifdef USING_MUMPS
  if (fMumpsInitialized) {
    // MUMPS memory release phase
    fMumpsData.job = JOB_END;

    // Call MUMPS to release the internal memory
    CallMumps();

    // Check for errors
    if (fMumpsData.info[0] < 0) {
      std::cerr << "MUMPS Error during memory deallocation: "
                << "INFO(1) = " << fMumpsData.info[0]
                << ", INFO(2) = " << fMumpsData.info[1] << std::endl;
      Error_check(MUMPS_INT(fMumpsData.info[0]), fMumpsData.info[1]);
      DebugStop();
    }

    fMumpsInitialized = false;
  }
#else
  NOMUMPS
#endif
}

template <class TVar>
TPZMumpsSolver<TVar>::~TPZMumpsSolver() {
  FreeMumpsMemory();
  // we should NOT delete fSymmetricSystem and fNonSymmetricSystem
}

template <class TVar>
void TPZMumpsSolver<TVar>::SetMatrix(TPZAutoPointer<TPZBaseMatrix> refmat) {
  auto *symSystem =
      dynamic_cast<TPZMatrix<TVar> *>(refmat.operator->());
  auto *nSymSystem =
      dynamic_cast<TPZMatrix<TVar> *>(refmat.operator->());
#ifdef PZDEBUG
  if (!symSystem && !nSymSystem) {
    PZError << __PRETTY_FUNCTION__;
    PZError << "This solver is only compatible with sparse matrices.\nAborting...\n";
    DebugStop();
  }
#endif

  fDecomposed = refmat->IsDecomposed();
  const MProperty prop = refmat->IsDefPositive()
                             ? MProperty::EPositiveDefinite
                             : MProperty::EIndefinite;
  SetMatrixType(refmat->GetSymmetry(), prop);
  TPZMatrixSolver<TVar>::SetMatrix(refmat);
}

template <class TVar>
void TPZMumpsSolver<TVar>::Solve(const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol, TPZFMatrix<TVar> *res) {

  if (!this->Matrix()) {
    PZError << __PRETTY_FUNCTION__;
    PZError << " called without a matrix pointer\n";
    DebugStop();
  }
  if (!fDecomposed) Decompose();
  sol = rhs;
  Solve(this->Matrix().operator->(), rhs, sol);

  if (res) res->Redim(rhs.Rows(), rhs.Cols());
}

template <class TVar>
void TPZMumpsSolver<TVar>::Decompose() {
  // Cast directly to the concrete MUMPS matrix types, skipping the
  // intermediate parent-class cast that was previously done here.
  auto *symSystem =
      dynamic_cast<TPZSYsmpMatrixMumps<TVar> *>(this->Matrix().operator->());
  auto *nSymSystem =
      dynamic_cast<TPZFYsmpMatrixMumps<TVar> *>(this->Matrix().operator->());
  if (symSystem) {
    Decompose(symSystem);
  } else if (nSymSystem) {
    Decompose(nSymSystem);
  } else {
    PZError << __PRETTY_FUNCTION__;
    PZError << "This solver is only compatible with TPZSYsmpMatrixMumps or "
               "TPZFYsmpMatrixMumps.\nAborting...\n";
    DebugStop();
  }
}

template <class TVar>
void TPZMumpsSolver<TVar>::Decompose(TPZMatrix<TVar> *mat) {
#ifndef USING_MUMPS
  NOMUMPS
#else
  auto *symSystem = dynamic_cast<TPZSYsmpMatrixMumps<TVar> *>(mat);
  auto *nSymSystem = dynamic_cast<TPZFYsmpMatrixMumps<TVar> *>(mat);

  MUMPS_INT n = 0;    // number of equations
  TVar *a;            // array of non zero values
  MUMPS_INT *ia, *ja; // row and column indices

  if (symSystem) { // symmetric matrix
    if (symSystem->Rows() == 0) {
      return;
    }
    a = &(symSystem->fA[0]);
    ia = (MUMPS_INT *)&(symSystem->fIA[0]);
    ja = (MUMPS_INT *)&(symSystem->fJA[0]);
    n = symSystem->Rows();
  }
  if (nSymSystem) { // non symmetric matrix
    a = &(nSymSystem->fA[0]);
    ia = (MUMPS_INT *)&(nSymSystem->fIA[0]);
    ja = (MUMPS_INT *)&(nSymSystem->fJA[0]);
    n = nSymSystem->Rows();
  }

  // Initialize MUMPS control parameters
  if (!fMumpsInitialized) {
    fMumpsData.comm_fortran = USE_COMM_WORLD; // USE_COMM_WORLD in C
    fMumpsData.par = MUMPS_HOST_PAR;          // host participates in the factorization -
                                              // This allows MUMPS to run on a single
                                              // processor and prevents the host processor
                                              // being idle during the factorization and
                                              // solve phases. MUMPS generally recommends par=1
    fMumpsData.sym = fMatrixType;             // 0: unsymmetric
                                              // 1: symmetric positive definite
                                              // 2: symmetric indefinite
    fMumpsData.job = JOB_INIT;                // Initialize MUMPS

    // Call MUMPS to initialize the internal data structures
    CallMumps();

    if (fMumpsData.info[0] < 0) {
      std::cerr << "MUMPS Error during initialization: "
                << "INFO(1) = " << fMumpsData.info[0]
                << ", INFO(2) = " << fMumpsData.info[1] << std::endl;
      Error_check(MUMPS_INT(fMumpsData.info[0]), fMumpsData.info[1]);
      DebugStop();
    }
    fMumpsInitialized = true;
  }

  MUMPS_INT nnz = (symSystem) ? symSystem->fA.size() : nSymSystem->fA.size();

  // Get COO format from matrix (conversion done once in matrix class)
  TPZVec<MUMPS_INT> irn, jcn;
  if (symSystem) {
    fMumpsData.irn = symSystem->fIRN1Based.begin();
    fMumpsData.jcn = symSystem->fJCN1Based.begin();
  } else {
    fMumpsData.irn = nSymSystem->fIRN1Based.begin();
    fMumpsData.jcn = nSymSystem->fJCN1Based.begin();
  }

  fMumpsData.n = n;
  fMumpsData.nz = nnz;
  fMumpsData.a = reinterpret_cast<decltype(fMumpsData.a)>(a);

  /**
   * INCTL(1): is the output stream for error messages
   *
   * Possible values:
   * ≤ 0: these messages will be suppressed.
   * > 0 : is the output stream.
   * Default value: 6 (standard output stream)
   */
  fMumpsData.icntl[0] = fICNTL[0].has_value() ? fICNTL[0].value() : (fMessageLevel ? 6 : -1);

  /**
   * INCTL(2): is the output stream for diagnostic printing and statistics local to each MPI process.
   *
   * Possible values:
   * ≤ 0: these messages will be suppressed.
   * > 0 : is the output stream.
   * Default value: 0
   * Remarks: If ICNTL(2) > 0 and ICNTL(4) ≥ 2, then information on advancement (flops done)
   * is also printed.
   */
  fMumpsData.icntl[1] = fICNTL[1].has_value() ? fICNTL[1].value() : (fMessageLevel ? 6 : -1);

  /**
   * ICNTL(3) is the output stream for global information, collected on the host.
   *
   * Possible values:
   * ≤ 0: these messages will be suppressed.
   * > 0 : is the output stream.
   * Default value: 6 (standard output stream)
   */
  fMumpsData.icntl[2] = fICNTL[2].has_value() ? fICNTL[2].value() : (fMessageLevel > 0 ? 6 : -1);

  /**
   * ICNTL(4) is the level of printing for error, warning, and diagnostic messages.
   *
   * Possible values:
   * ≤ 0: No messages output.
   * 1 : Only error messages printed.
   * 2 : Errors, warnings, and main statistics printed.
   * 3 : Errors and warnings and terse diagnostics (only first ten entries of arrays) printed.
   * ≥ 4 : Errors, warnings and information on input, output parameters printed.
   * Default value: 2 (errors, warnings and main statistics printed)
   */
  fMumpsData.icntl[3] = fICNTL[3].has_value() ? fICNTL[3].value() : (fMessageLevel);

  /**
   * ICNTL(5) controls the matrix input format
   *
   * Phase: accessed by the host and only during the analysis phase
   * Possible variables/arrays involved: N, NNZ (or NZ for backward compatibility), IRN, JCN,
   *    NNZ loc (or NZ loc for backward compatibility), IRN loc, JCN loc, A loc, NELT, ELTPTR,
   *    ELTVAR, and A ELT
   *
   * Possible values:
   * 0: assembled format
   * 1: elemental format. The matrix must be input in the structure components N, NELT, ELTPTR, ELTVAR, and A ELT
   */
  fMumpsData.icntl[4] = fICNTL[4].has_value() ? fICNTL[4].value() : (0);

  /**
   * ICNTL(6) computes a permutation to permute the matrix to a zero-free diagonal and/or computes a matrix scaling
   * Phase: accessed by the host and only during sequential analysis
   * Possible values:
   * 0 : No column permutation is computed.
   * 1 : The permuted matrix has as many entries on its diagonal as possible. The values on the diagonal are of arbitrary size.
   * 2 : The permutation is such that the smallest value on the diagonal of the permuted matrix is maximized. The numerical values of the original matrix, (mumps par%A), must be provided by the user during the analysis phase.
   * 3 : Variant of option 2 with different performance. The numerical values of the original matrix (mumps par%A) must be provided by the user during the analysis phase.
   * 4 : The sum of the diagonal entries of the permuted matrix is maximized. The numerical values of the original matrix (mumps par%A) must be provided by the user during the analysis phase.
   * 5 : The product of the diagonal entries of the permuted matrix is maximized.
   * 6 : Similar to 5 but with a more costly (time and memory footprint) algorithm. The numerical values of the original matrix, mumps par%A, must be provided by the user during the analysis phase.
   * 7 : Based on the structural symmetry of the input matrix and on the availability of the numerical values, the value of ICNTL(6) is automatically chosen by the software.
   *
   * Default value: 7 (automatic choice done by the package)
   */
  fMumpsData.icntl[5] = fICNTL[5].has_value() ? fICNTL[5].value() : (7);

  /**
   * ICNTL(7) computes a symmetric permutation (ordering) to determine the pivot order to be used for the factorization
   *
   * Possible values:
   * 0: Approximate Minimum Degree (AMD) is used
   * 1: The pivot order should be set by the user in PERM IN, on the host processor
   * 2: Approximate Minimum Fill (AMF) is used
   * 3: SCOTCH is used if installed by user, otherwise treats as 7
   * 4: PORD is used if installed by user, otherwise treats as 7
   * 5: METIS is used if installed by user, otherwise treats as 7
   * 6: Approximate Minimum Degree with automatic quasi-dense row detection (QAMD) is used
   * 7: Automatic choice done by the package
   */
  fMumpsData.icntl[6] = fICNTL[6].has_value() ? fICNTL[6].value() : (7);

  /**
   * ICNTL(13) controls the parallelism of the root node (enabling or not the use of ScaLAPACK) and also its splitting.
   *
   * Phase: accessed by the host during the analysis phase.
   * Possible values :
   * < -1 : treated as 0.
   * -1 : force splitting of the root node in all cases (even sequentially)
   * 0 : parallel factorization of the root node based on ScaLAPACK.
   * > 0 : ScaLAPACK is not used (recommended value is 1 to partly recover parallelism of the root node). It forces a sequential factorization of the root node (ScaLAPACK will not be used).
   *
   * Default value: 0 (parallel factorization on the root node)
   */
  fMumpsData.icntl[12] = fICNTL[12].has_value() ? fICNTL[12].value() : (0);

  /**
   * ICNTL(14) controls the percentage increase in the estimated working space
   *
   * Phase: accessed by the host both during the analysis and the factorization phases.
   * Default value: between 20 and 35 (which corresponds to at most 35 % increase) and depends on
   * the number of MPI processes. It is set to 5 % with SYM=1 and one MPI process.
   */
  fMumpsData.icntl[13] = fICNTL[13].has_value() ? fICNTL[13].value() : (35);

  /**
   * ICNTL(18) defines the strategy for the distributed input matrix (only for assembled matrix).
   *
   * Phase: accessed by the host during the analysis phase.
   * Possible values :
   * 0 : the input matrix is centralized on the host (see Subsection 5.4.2.1).
   * 1 : the user provides the structure of the matrix on the host at analysis, MUMPS returns a mapping and the user should then provide the matrix entries distributed according to the mapping on entry to the numerical factorization phase.
   * 2 : the user provides the structure of the matrix on the host at analysis, and the distributed matrix entries on all slave processors at factorization. Any distribution is allowed
   * 3 : user directly provides the distributed matrix, pattern and entries, input both for analysis and factorization.
   *
   * Other values are treated as 0.
   * Default value: 0 (input matrix centralized on the host)
   */
  fMumpsData.icntl[17] = fICNTL[17].has_value() ? fICNTL[17].value() : (0);

  /**
   * ICNTL(28) determines whether a sequential or parallel computation of the ordering is performed
   *
   * Phase: accessed by the host process during the analysis phase.
   * Possible values:
   * 0: automatic choice.
   * 1: sequential computation. In this case the ordering method is set by ICNTL(7) and the ICNTL(29) parameter is meaningless (choice of the parallel ordering tool).
   * 2: parallel computation. A parallel ordering and parallel symbolic factorization is requested by the user.
   *
   * Any other values will be treated as 0.
   * Default value: 0 (automatic choice)
   */
  fMumpsData.icntl[27] = fICNTL[27].has_value() ? fICNTL[27].value() : (0);

  if (!fCustomSettings) {
    /**
     * CNTL(1) is the relative threshold for numerical pivoting.
     *
     * Phase: accessed by the host during the factorization phase.
     * Possible values :
     * < 0.0: Automatic choice
     * = 0.0: no numerical pivoting performed and the subroutine will fail if a zero pivot is encountered.
     * > 0.0: numerical pivoting performed.
     * For unsymmetric matrices values greater than 1.0 are treated as 1.0
     * For symmetric matrices values greater than 0.5 are treated as 0.5
     *
     * Default value: -1.0 (automatic choice):
     * 0.1: in case of rank-revealing (ICNTL(56)= 1)
     * 0.01: for unsymmetric or general symmetric matrices
     * 0.0: for symmetric positive definite matrices
     */
    fMumpsData.cntl[0] = -1.0;
  }

  /**
   * ICNTL(24) controls the detection of “null pivot rows”.
   * Phase: accessed by the host during the factorization phase
   * Possible variables/arrays involved: PIVNUL LIST
   * Possible values :
   * 0: Nothing done. A null pivot row will result in error INFO(1)=-10.
   * 1: Null pivot row detection.
   * Other values are treated as 0.
   * Default value: 0 (no null pivot row detection)
   */
  if (fProperty == MProperty::EIndefinite) {
    fMumpsData.icntl[23] = fICNTL[23].has_value() ? fICNTL[23].value() : (1);
  }

  // --- Analysis phase
  fMumpsData.job = JOB_ANALYSIS;
  CallMumps();

  if (fMumpsData.info[0] < 0) {
    std::cerr << "MUMPS analysis error: "
              << "INFO(1) = " << fMumpsData.info[0]
              << ", INFO(2) = " << fMumpsData.info[1] << std::endl;
    Error_check(MUMPS_INT(fMumpsData.info[0]), fMumpsData.info[1]);
    DebugStop();
  }

  // --- Factorization phase
  fMumpsData.job = JOB_FACTORIZE;
  CallMumps();

  if (fMumpsData.info[0] < 0) {
    std::cerr << "MUMPS factorization error: "
              << "INFO(1) = " << fMumpsData.info[0]
              << ", INFO(2) = " << fMumpsData.info[1] << std::endl;
    Error_check(MUMPS_INT(fMumpsData.info[0]), fMumpsData.info[1]);
    DebugStop();
  }

  fDecomposed = true;
#endif
}

/// Use the decomposed matrix to invert the system of equations
template <class TVar>
void TPZMumpsSolver<TVar>::Solve(const TPZMatrix<TVar> *mat, const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol) const {
#ifndef USING_MUMPS
  NOMUMPS
#else

#ifdef PZDEBUG
  if (!fDecomposed) {
    PZError << __PRETTY_FUNCTION__;
    PZError << "\nError: Matrix has not been decomposed.\nAborting..." << std::endl;
    DebugStop();
  }
#endif

  auto *symSystem = dynamic_cast<const TPZSYsmpMatrixMumps<TVar> *>(mat);
  auto *nSymSystem = dynamic_cast<const TPZFYsmpMatrixMumps<TVar> *>(mat);

  long long n = 0;
  TVar *a;
  long long *ia, *ja;

  if (symSystem) {
    if (symSystem->Rows() == 0) {
      return;
    }
    a = &(symSystem->fA[0]);
    ia = (long long *)&(symSystem->fIA[0]);
    ja = (long long *)&(symSystem->fJA[0]);
    n = symSystem->Rows();
  }
  if (nSymSystem) {
    a = &(nSymSystem->fA[0]);
    ia = (long long *)&(nSymSystem->fIA[0]);
    ja = (long long *)&(nSymSystem->fJA[0]);
    n = nSymSystem->Rows();
  }

#ifdef PZ_LOG
  if (logger.isDebugEnabled()) {
    std::stringstream sout;
    sout << "MUMPS control parameters:\n";
    for (int i = 0; i < 40; i++) {
      if (fMumpsData.icntl[i] != 0) {
        sout << "ICNTL(" << i + 1 << ") = " << fMumpsData.icntl[i] << "\n";
      }
    }
    LOGPZ_DEBUG(logger, sout.str())
  }
#endif

  long long nrhs = rhs.Cols();
  n = rhs.Rows();

  auto normb = Norm(rhs);
  if (IsZero(normb)) {
    sol.Zero();
    return;
  }

  // Prepare RHS - MUMPS needs a non-const pointer
  TPZFMatrix<TVar> rhs_copy = rhs;
  sol.Redim(n, nrhs);

  // Configure RHS and solution in MUMPS data structure
  fMumpsData.nrhs = nrhs;
  fMumpsData.lrhs = n;
  fMumpsData.rhs = reinterpret_cast<decltype(fMumpsData.rhs)>(&rhs_copy.g(0, 0));

  // MUMPS overwrites RHS with the solution, so copy to sol afterwards

  // Check if MUMPS was properly initialized and decomposed
  if (!fMumpsInitialized) {
    std::cerr << "MUMPS Error: MUMPS was not initialized before solving.\n";
    std::cerr << "fMumpsInitialized = " << fMumpsInitialized << std::endl;
    DebugStop();
  }

  if (!fDecomposed) {
    std::cerr << "MUMPS Error: Matrix was not decomposed before solving.\n";
    std::cerr << "fDecomposed = " << fDecomposed << std::endl;
    DebugStop();
  }

  // Solve phase
  fMumpsData.job = JOB_SOLVE;

  CallMumps();

  // Check for errors
  if (fMumpsData.info[0] < 0) {
    std::cerr << "MUMPS solve error: "
              << "INFO(1) = " << fMumpsData.info[0]
              << ", INFO(2) = " << fMumpsData.info[1] << std::endl;
    Error_check(MUMPS_INT(fMumpsData.info[0]), fMumpsData.info[1]);

    // Try refactorization if there was an error
    if (fMumpsData.info[0] == -10 || fMumpsData.info[0] == -13) {
      std::cout << "MUMPS: Calling numerical factorization due to error...\n";

      // Redo analysis and factorization
      fMumpsData.job = JOB_ANALYSIS; // Analysis
      CallMumps();

      if (fMumpsData.info[0] < 0) {
        std::cerr << "MUMPS re-analysis error: "
                  << "INFO(1) = " << fMumpsData.info[0] << std::endl;
        Error_check(MUMPS_INT(fMumpsData.info[0]), fMumpsData.info[1]);
        DebugStop();
      }

      fMumpsData.job = JOB_FACTORIZE; // Factorization
      CallMumps();

      if (fMumpsData.info[0] < 0) {
        std::cerr << "MUMPS re-factorization error: "
                  << "INFO(1) = " << fMumpsData.info[0] << std::endl;
        Error_check(MUMPS_INT(fMumpsData.info[0]), fMumpsData.info[1]);
        DebugStop();
      }

      // Try solving again
      fMumpsData.rhs = reinterpret_cast<decltype(fMumpsData.rhs)>(&rhs_copy.g(0, 0));
      fMumpsData.job = JOB_SOLVE;
      CallMumps();

      if (fMumpsData.info[0] < 0) {
        std::cerr << "MUMPS solve error after refactorization: "
                  << "INFO(1) = " << fMumpsData.info[0] << std::endl;
        Error_check(MUMPS_INT(fMumpsData.info[0]), fMumpsData.info[1]);
        DebugStop();
      }
    } else {
      DebugStop();
    }
  }

  // Copy solution (MUMPS overwrites RHS with the solution)
  for (long long i = 0; i < n; i++) {
    for (long long j = 0; j < nrhs; j++) {
      sol(i, j) = rhs_copy(i, j);
    }
  }

  // Report useful information
  if (fMessageLevel > 0) {
    // INFO(3): estimated number of operations for elimination (millions)
    if (fMumpsData.info[2] > 0) {
      std::cout << "MUMPS: Estimated operations: "
                << fMumpsData.info[2] << " million flops\n";
    }

    // INFO(9): real size of workspace used
    if (fMumpsData.info[8] > 0) {
      std::cout << "MUMPS: Real workspace used: "
                << fMumpsData.info[8] << " MB\n";
    }

    // RINFOG(3): number of operations during elimination
    if (fMumpsData.rinfog[2] > 0) {
      std::cout << "MUMPS: Actual operations: "
                << fMumpsData.rinfog[2] << " flops\n";
    }
  }

  // Warnings about solution quality
  if (fMumpsData.info[0] > 0) {
    switch (fMumpsData.info[0]) {
    case 1:
      std::cout << "MUMPS warning: Index out of range.\n";
      break;
    case 2:
      std::cout << "MUMPS warning: Allocated array too small.\n";
      break;
    case 4:
      std::cout << "MUMPS warning: Matrix is singular in structure.\n";
      break;
    case 8:
      std::cout << "MUMPS warning: Main internal integer workarray too small.\n";
      break;
    default:
      if (fMessageLevel > 1) {
        std::cout << "MUMPS warning code: " << fMumpsData.info[0] << "\n";
      }
      break;
    }
  }

  // #ifdef PZDEBUG
  //   std::cout << "MUMPS: linear solve complete.\n";
  // #endif
#endif // USING_MUMPS
}

template <class TVar>
void TPZMumpsSolver<TVar>::SetMessageLevel(int lvl) {
  fMessageLevel = lvl;
}

template <class TVar>
TPZMumpsSolver<TVar> *TPZMumpsSolver<TVar>::Clone() const {
  // Cannot clone a MUMPS solver - each instance needs its own MUMPS initialization
  // Return a new uninitialized solver instead
  auto *newSolver = new TPZMumpsSolver<TVar>();
  newSolver->fSymmetry = fSymmetry;
  newSolver->fProperty = fProperty;
  newSolver->fMatrixType = fMatrixType;
  newSolver->fMessageLevel = fMessageLevel;
  newSolver->fCustomSettings = fCustomSettings;
  // Note: fMumpsData, fDecomposed, and fMumpsInitialized are NOT copied
  // The cloned solver will need to be initialized and decomposed separately
  return newSolver;
}

/**
 * Sets a specific ICNTL parameter to a custom value.
 * Note: MUMPS uses 1-based indexing for ICNTL parameters.
 * @param index Index of the ICNTL parameter (1 to 60)
 * @param value Value to set for the ICNTL parameter
 */
template <class TVar>
void TPZMumpsSolver<TVar>::SetICNTL(int index, MUMPS_INT value) {
  const int max_size = 60; // ICNTL array size in MUMPS
  if (index < 1 || index > max_size) {
    PZError << __PRETTY_FUNCTION__
            << "\nInvalid ICNTL index!"
            << "\nExpected index between 1 and " << max_size
            << " but got " << index
            << "\nAborting..." << std::endl;
    DebugStop();
  }
  // MUMPS uses 1-based indexing, C++ uses 0-based
  fICNTL[index - 1] = value;
  fCustomSettings = true;
}

template <class TVar>
void TPZMumpsSolver<TVar>::ResetICNTL() {
  // will setup param with info regarding matrix' structure
  this->MatrixType();
  fCustomSettings = false;
}

template <class TVar>
void TPZMumpsSolver<TVar>::SetMatrixType(SymProp symtype, MProperty prop) {

  fSymmetry = symtype;
  fProperty = prop;
  fMatrixType = MatrixType();
}

template <class TVar>
long long TPZMumpsSolver<TVar>::MatrixType() {
  // Determine matrix symmetry type for MUMPS
  // MUMPS SYM parameter:
  // 0 = unsymmetric
  // 1 = symmetric positive definite
  // 2 = symmetric general (indefinite)

  if constexpr (is_complex<TVar>::value) {
    switch (fSymmetry) {
    case SymProp::NonSym:
      fMatrixType = 0; // Unsymmetric
      break;
    case SymProp::Sym:
      fMatrixType = 2; // Symmetric complex (general)
      break;
    case SymProp::Herm:
      if (fProperty == MProperty::EPositiveDefinite) {
        fMatrixType = 1; // Hermitian positive definite
      } else {
        fMatrixType = 2; // Hermitian general
      }
      break;
    }
  } else {
    // Real matrices (double)
    switch (fSymmetry) {
    case SymProp::NonSym:
      fMatrixType = 0; // Unsymmetric
      break;
    case SymProp::Sym:
    case SymProp::Herm: // Same as Sym for real matrices
      if (fProperty == MProperty::EPositiveDefinite) {
        fMatrixType = 1; // Symmetric positive definite
      } else {
        fMatrixType = 2; // Symmetric general (indefinite)
      }
      break;
    }
  }

  // ICNTL(1-3): Output streams
  fMumpsData.icntl[0] = fICNTL[0].has_value() ? fICNTL[0].value() : (fMessageLevel > 0 ? 6 : -1); // Error messages
  fMumpsData.icntl[1] = fICNTL[1].has_value() ? fICNTL[1].value() : (fMessageLevel > 1 ? 6 : -1); // Diagnostics
  fMumpsData.icntl[2] = fICNTL[2].has_value() ? fICNTL[2].value() : (fMessageLevel > 0 ? 6 : -1); // Global info

  // ICNTL(4): Print level
  fMumpsData.icntl[3] = fICNTL[3].has_value() ? fICNTL[3].value() : (fMessageLevel);

  // ICNTL(5): Matrix format (0 = assembled/elemental)
  fMumpsData.icntl[4] = fICNTL[4].has_value() ? fICNTL[4].value() : (0);

  // ICNTL(6): Permutation for input matrix (0 = automatic)
  fMumpsData.icntl[5] = fICNTL[5].has_value() ? fICNTL[5].value() : (0);

  // ICNTL(7): Ordering strategy
  // 7 = automatic, 5 = METIS, 0 = AMD
  fMumpsData.icntl[6] = fICNTL[6].has_value() ? fICNTL[6].value() : (7);

  // ICNTL(8): Scaling strategy
  // 77 = automatic, -1 = none
  fMumpsData.icntl[7] = fICNTL[7].has_value() ? fICNTL[7].value() : (77);

  // ICNTL(10): Max steps of iterative refinement
  fMumpsData.icntl[9] = fICNTL[9].has_value() ? fICNTL[9].value() : (0); // Default: no iterative refinement

  // ICNTL(11): Error analysis
  // 0 = no statistics, 1 = compute statistics
  fMumpsData.icntl[10] = fICNTL[10].has_value() ? fICNTL[10].value() : (1);

  // ICNTL(12): Ordering for symmetric matrices
  // 0 = automatic, 1 = use A+A^T
  fMumpsData.icntl[11] = fICNTL[11].has_value() ? fICNTL[11].value() : (0);

  // ICNTL(13): ScaLAPACK (sequential: always 0)
  fMumpsData.icntl[12] = fICNTL[12].has_value() ? fICNTL[12].value() : (0);

  // ICNTL(14): Percentage increase in estimated working space
  fMumpsData.icntl[13] = fICNTL[13].has_value() ? fICNTL[13].value() : (20);

  // ICNTL(18): Distributed matrix strategy
  // 0 = centralized on host (sequential)
  fMumpsData.icntl[17] = fICNTL[17].has_value() ? fICNTL[17].value() : (0);

  // ICNTL(19): Schur complement (0 = no Schur)
  fMumpsData.icntl[18] = fICNTL[18].has_value() ? fICNTL[18].value() : (0);

  // ICNTL(20): RHS format (0 = dense)
  fMumpsData.icntl[19] = fICNTL[19].has_value() ? fICNTL[19].value() : (0);

  // ICNTL(21): Solution format (0 = centralized)
  fMumpsData.icntl[20] = fICNTL[20].has_value() ? fICNTL[20].value() : (0);

  // ICNTL(22): Out-of-core factorization (0 = in-core)
  fMumpsData.icntl[21] = fICNTL[21].has_value() ? fICNTL[21].value() : (0);

  // ICNTL(23): Maximum memory for working space (0 = automatic)
  fMumpsData.icntl[22] = fICNTL[22].has_value() ? fICNTL[22].value() : (0);

  // ICNTL(24): Null space detection (0 = no)
  fMumpsData.icntl[23] = fICNTL[23].has_value() ? fICNTL[23].value() : (0);

  // ICNTL(28): Parallel ordering (0 = sequential, 1 = PT-SCOTCH, 2 = ParMetis)
  // For sequential MUMPS, always use 0
  fMumpsData.icntl[27] = fICNTL[27].has_value() ? fICNTL[27].value() : (0);

  // ICNTL(29): Parallel ordering method details
  fMumpsData.icntl[28] = fICNTL[28].has_value() ? fICNTL[28].value() : (0);

  if (!fCustomSettings) {
    // CNTL(1): Relative threshold for numerical pivoting
    // Default 0.01 for unsymmetric, 0.0 for symmetric positive definite
    if (fMatrixType == 0) { // Unsymmetric
      fMumpsData.cntl[0] = 0.01;
    } else if (fMatrixType == 1) { // Symmetric positive definite
      fMumpsData.cntl[0] = 0.0;
    } else { // Symmetric indefinite
      fMumpsData.cntl[0] = 0.01;
    }

    // CNTL(2): Stopping criterion for iterative refinement
    // (only if ICNTL(10) > 0)
    fMumpsData.cntl[1] = 1.0e-8;

    // CNTL(3): Threshold for null pivot detection
    fMumpsData.cntl[2] = 0.0;

    // CNTL(4): Threshold for static pivoting (advanced)
    fMumpsData.cntl[3] = -1.0; // Disabled by default

    // CNTL(5): Fixation for null pivots
    fMumpsData.cntl[4] = 0.0;
  }

  return fMatrixType;
}

template <class TVar>
int DataType([[maybe_unused]] TVar a) {
  DebugStop();
  return 0;
}

template <>
inline int DataType([[maybe_unused]] double a) {
  return 0;
}

template <>
inline int DataType([[maybe_unused]] float a) {
  return 1;
}

template <>
inline int DataType([[maybe_unused]] std::complex<double> a) {
  return 0;
}

template <>
inline int DataType([[maybe_unused]] std::complex<float> a) {
  return 1;
}

/**
 * @brief Check MUMPS error codes and print corresponding messages.
 * @param INFO1 The primary MUMPS error code (INFO(1)).
 * @param INFO2 Secondary information (INFO(2)).
 */
template <typename TVar>
void Error_check(MUMPS_INT INFO1, TVar INFO2) {
  // MUMPS error codes are stored in INFO(1)
  // Negative values indicate errors
  // Positive values indicate warnings
  // INFO(2) provides additional context for some errors

  switch (INFO1) {
  case -1:
    std::cout << "MUMPS error -1: An error occurred on processor " << INFO2 << "." << std::endl;
    break;
  case -2:
    std::cout << "MUMPS error -2: NZ is out of range. NZ = " << INFO2 << " (should be >= 1)." << std::endl;
    break;
  case -3:
    std::cout << "MUMPS error -3: MUMPS was called with an invalid value for JOB. JOB = " << INFO2 << "." << std::endl;
    break;
  case -4:
    std::cout << "MUMPS error -4: Error in user-provided permutation array PERM_IN at position " << INFO2 << "." << std::endl;
    break;
  case -5:
    std::cout << "MUMPS error -5: Problem with real workspace allocation during analysis. Required size = " << INFO2 << " MB." << std::endl;
    break;
  case -6:
    std::cout << "MUMPS error -6: Matrix is singular in structure. " << INFO2 << " entries are missing from diagonal." << std::endl;
    break;
  case -7:
    std::cout << "MUMPS error -7: Problem with integer workspace allocation during analysis. Required size = " << INFO2 << " MB." << std::endl;
    break;
  case -8:
    std::cout << "MUMPS error -8: Main internal integer workarray is too small. Minimum required = " << INFO2 << "." << std::endl;
    break;
  case -9:
    std::cout << "MUMPS error -9: Main internal real/complex workarray is too small. Minimum required = " << INFO2 << " MB." << std::endl;
    break;
  case -10:
    std::cout << "MUMPS error -10: Numerically singular matrix. Pivot " << INFO2 << " is zero or too small." << std::endl;
    break;
  case -11:
    std::cout << "MUMPS error -11: Internal real/complex workarray is too small for solution phase. Minimum required = " << INFO2 << " MB." << std::endl;
    break;
  case -12:
    std::cout << "MUMPS error -12: Internal real/complex workarray is too small for iterative refinement/error analysis. Minimum required = " << INFO2 << " MB." << std::endl;
    break;
  case -13:
    std::cout << "MUMPS error -13: Problem with memory allocation during factorization/solve. Required additional memory = " << INFO2 << " MB." << std::endl;
    break;
  case -14:
    std::cout << "MUMPS error -14: Integer workarray is too small for solution phase. Minimum required = " << INFO2 << "." << std::endl;
    break;
  case -15:
    std::cout << "MUMPS error -15: Integer workarray is too small for iterative refinement/error analysis. Minimum required = " << INFO2 << "." << std::endl;
    break;
  case -16:
    std::cout << "MUMPS error -16: N is out of range. N = " << INFO2 << " (should be >= 1)." << std::endl;
    break;
  case -17:
    std::cout << "MUMPS error -17: Internal send buffer is too small (parallel version only). Required size = " << INFO2 << "." << std::endl;
    break;
  case -18:
    std::cout << "MUMPS error -18: Blocking size for multiple RHS is too large. ICNTL(27) = " << INFO2 << "." << std::endl;
    break;
  case -19:
    std::cout << "MUMPS error -19: Maximum allowed size of working memory is too small. Required = " << INFO2 << " MB." << std::endl;
    break;
  case -20:
    std::cout << "MUMPS error -20: Recv buffer is too small (parallel version only). Required size = " << INFO2 << "." << std::endl;
    break;
  case -21:
    std::cout << "MUMPS error -21: Value of PAR parameter is not valid. PAR = " << INFO2 << "." << std::endl;
    break;
  case -22:
    std::cout << "MUMPS error -22: A pointer array is provided but constraints are violated. Problem at position " << INFO2 << "." << std::endl;
    break;
  case -23:
    std::cout << "MUMPS error -23: MPI was not initialized." << std::endl;
    break;
  case -24:
    std::cout << "MUMPS error -24: NELT is out of range. NELT = " << INFO2 << " (should be >= 1)." << std::endl;
    break;
  case -25:
    std::cout << "MUMPS error -25: Problem with Schur complement ordering. Invalid entry in LISTVAR_SCHUR at position " << INFO2 << "." << std::endl;
    break;
  case -26:
    std::cout << "MUMPS error -26: Null pivot detected during Schur complement factorization. Pivot row = " << INFO2 << "." << std::endl;
    break;
  case -27:
    std::cout << "MUMPS error -27: Problem with Schur complement I/O. Error during operation " << INFO2 << "." << std::endl;
    break;
  case -28:
    std::cout << "MUMPS error -28: NRHS is out of range. NRHS = " << INFO2 << " (should be >= 1)." << std::endl;
    break;
  case -29:
    std::cout << "MUMPS error -29: LRHS is out of range. LRHS = " << INFO2 << ", should be >= max(1, N) = " << std::max(MUMPS_INT(1), static_cast<MUMPS_INT>(INFO2)) << "." << std::endl;
    break;
  case -30:
    std::cout << "MUMPS error -30: NRHS or LRHS incompatible with JOB. JOB = " << INFO2 << "." << std::endl;
    break;
  case -31:
    std::cout << "MUMPS error -31: NZ_RHS is out of range. NZ_RHS = " << INFO2 << " (should be >= 1)." << std::endl;
    break;
  case -32:
    std::cout << "MUMPS error -32: JOB=8 forward elimination not supported with certain settings. ICNTL(" << INFO2 << ") incompatible." << std::endl;
    break;
  case -33:
    std::cout << "MUMPS error -33: ICNTL(26) incompatible with ICNTL(20) and/or ICNTL(21). Conflicting parameter = ICNTL(" << INFO2 << ")." << std::endl;
    break;
  case -34:
    std::cout << "MUMPS error -34: LREDRHS is out of range. LREDRHS = " << INFO2 << "." << std::endl;
    break;
  case -35:
    std::cout << "MUMPS error -35: Problem with Schur matrix I/O. Error code = " << INFO2 << "." << std::endl;
    break;
  case -36:
    std::cout << "MUMPS error -36: Incompatible values of ICNTL(25) and INFOG(28). ICNTL(25) = " << INFO2 << "." << std::endl;
    break;
  case -37:
    std::cout << "MUMPS error -37: Value of ICNTL(25) incompatible with other parameters. ICNTL(25) = " << INFO2 << "." << std::endl;
    break;
  case -38:
    std::cout << "MUMPS error -38: Parallel analysis error (distributed input matrix). Failure on processor " << INFO2 << "." << std::endl;
    break;
  case -39:
    std::cout << "MUMPS error -39: Problem with Schur complement. SIZE_SCHUR = " << INFO2 << " is incompatible." << std::endl;
    break;
  case -40:
    std::cout << "MUMPS error -40: Incompatible values between ICNTL(28) and ICNTL(5) and/or ICNTL(19). ICNTL(28) = " << INFO2 << "." << std::endl;
    break;
  case -44:
    std::cout << "MUMPS error -44: Error in solve phase (internal). Sub-error code = " << INFO2 << "." << std::endl;
    break;
  case -45:
    std::cout << "MUMPS error -45: NRHS incompatible with option for sparse RHS. NRHS = " << INFO2 << "." << std::endl;
    break;
  case -46:
    std::cout << "MUMPS error -46: NZ_RHS is larger than ICNTL(27)*N. NZ_RHS = " << INFO2 << "." << std::endl;
    break;
  case -47:
    std::cout << "MUMPS error -47: Problem with distributed solution. Processor " << INFO2 << " has invalid ISOL_loc or LSOL_loc." << std::endl;
    break;
  case -48:
    std::cout << "MUMPS error -48: A-inverse request error. Element (" << INFO2 << ") is invalid." << std::endl;
    break;
  case -49:
    std::cout << "MUMPS error -49: Problem with SIZE_SCHUR. SIZE_SCHUR = " << INFO2 << " is invalid." << std::endl;
    break;
  case -51:
    std::cout << "MUMPS error -51: Out-of-range value in permutation array on entry. Position " << INFO2 << " is invalid." << std::endl;
    break;
  case -52:
    std::cout << "MUMPS error -52: Inconsistency in input permutation array. Entry " << INFO2 << " appears more than once." << std::endl;
    break;
  case -53:
    std::cout << "MUMPS error -53: Problem with ICNTL(23) and out-of-core mode. ICNTL(23) = " << INFO2 << " incompatible." << std::endl;
    break;
  case -54:
    std::cout << "MUMPS error -54: Problem with element size. Element " << INFO2 << " has invalid size." << std::endl;
    break;
  case -55:
    std::cout << "MUMPS error -55: Problem with variable block low-rank parameters. Invalid parameter at position " << INFO2 << "." << std::endl;
    break;
  case -70:
    std::cout << "MUMPS error -70: Internal error in BLR (Block Low-Rank). Sub-error = " << INFO2 << "." << std::endl;
    break;
  case -71:
    std::cout << "MUMPS error -71: Problem with ICNTL(35) (BLR CB_SIZE parameter). ICNTL(35) = " << INFO2 << " is invalid." << std::endl;
    break;
  case -72:
    std::cout << "MUMPS error -72: Problem allocating BLR structures. Required memory = " << INFO2 << " MB." << std::endl;
    break;
  case -73:
    std::cout << "MUMPS error -73: NSLAVES is out of range (distributed input on master). NSLAVES = " << INFO2 << "." << std::endl;
    break;
  case -90:
    std::cout << "MUMPS error -90: Error in out-of-core management. Error code = " << INFO2 << "." << std::endl;
    break;

  // Positive values are warnings
  case 1:
    std::cout << "MUMPS warning +1: Index (in IRN or JCN) is out of range. " << INFO2 << " such entries detected." << std::endl;
    break;
  case 2:
    std::cout << "MUMPS warning +2: During error analysis, some integer statistical information is computed in double precision." << std::endl;
    break;
  case 4:
    std::cout << "MUMPS warning +4: Matrix is singular in structure. Rank = " << INFO2 << "." << std::endl;
    break;
  case 8:
    std::cout << "MUMPS warning +8: More than ICNTL(14) pivot adjustments were performed. Total adjustments = " << INFO2 << "." << std::endl;
    break;

  default:
    if (INFO1 < 0) {
      std::cout << "MUMPS error " << INFO1 << ": Unrecognized error code. INFO(2) = " << INFO2 << ". Check MUMPS documentation." << std::endl;
    } else if (INFO1 > 0) {
      std::cout << "MUMPS warning " << INFO1 << ": Unrecognized warning code. INFO(2) = " << INFO2 << ". Check MUMPS documentation." << std::endl;
    } else {
      std::cout << "MUMPS: No error (INFO(1) = 0)." << std::endl;
    }
    break;
  }
}

#ifdef MUMPS_HAVE_SINGLE
template class TPZMumpsSolver<float>;
#endif
#ifdef MUMPS_HAVE_DOUBLE
template class TPZMumpsSolver<double>;
template class TPZMumpsSolver<long double>;
#endif
#ifdef MUMPS_HAVE_COMPLEX
template class TPZMumpsSolver<std::complex<float>>;
#endif
#ifdef MUMPS_HAVE_COMPLEX16
template class TPZMumpsSolver<std::complex<double>>;
template class TPZMumpsSolver<std::complex<long double>>;
#endif