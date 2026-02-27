/**
 * @file TestSolverComparison.cpp
 * @brief Unit test comparing MUMPS and Pardiso sparse direct solvers.
 *
 * Solves the same 3D Darcy flow problem with both solvers and checks
 * that their solutions agree to within a tight tolerance.
 *
 * Based on the setup used in main_mumps_3d.cpp.
 *
 * --- Notes on the solver API ---
 *
 * TPZStepSolver::SetDirect(decomp) is required to configure the analysis to
 * use a direct solver. The decomp argument (ECholesky, ELU, …) is only used
 * as bookkeeping inside NeoPZ; it does NOT select the factorization algorithm
 * used by MUMPS or Pardiso. For those external solvers the factorization mode
 * is determined entirely by two other sources:
 *
 *   1. The matrix class itself:
 *        TPZSYsmpMatrixMumps / TPZSYsmpMatrixPardiso  → symmetric storage
 *        TPZFYsmpMatrixMumps / TPZFYsmpMatrixPardiso  → non-symmetric storage
 *
 *   2. Matrix::SetDefPositive(true/false):
 *        For symmetric matrices only — switches MUMPS from sym=2 (LDL^T) to
 *        sym=1 (Cholesky) and Pardiso from mtype=-2 to mtype=2. For
 *        non-symmetric matrices this call has no effect (sym=0 always).
 */

#if defined(PZ_USING_MKL) && defined(PZ_USING_MUMPS)

#include "DarcyFlow/TPZDarcyFlow.h"
#include "TPZGenGrid3D.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSSpStructMatrixMumps.h"
#include "TPZSpStructMatrix.h"
#include "TPZSpStructMatrixMumps.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzstepsolver.h"
#include "tpzautopointer.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// -----------------------------------------------------------------------
// Mesh helpers (identical to those in main_mumps_3d.cpp)
// -----------------------------------------------------------------------

namespace solvercomparison {

enum MatIds { EMatId = 1 };

static TPZAutoPointer<TPZGeoMesh> CreateGeoMesh3D(const int nDiv) {
    TPZVec<REAL> minX(3, 0.);
    TPZVec<REAL> maxX(3, 1.);
    TPZVec<int> nelDiv = {nDiv, nDiv, nDiv};
    MMeshType elType = MMeshType::ETetrahedral;
    TPZGenGrid3D gen3d(minX, maxX, nelDiv, elType);
    gen3d.BuildVolumetricElements(EMatId);
    // face material ids follow the convention used in main_mumps_3d.cpp
    TPZGeoMesh *gmesh = gen3d.BuildBoundaryElements(-1, -3, -1, -2, -1, -1);
    return TPZAutoPointer<TPZGeoMesh>(gmesh);
}

static TPZAutoPointer<TPZCompMesh> CreateCompMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                                                  const int pord) {
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(gmesh->Dimension());
    cmesh->SetDefaultOrder(pord);
    cmesh->SetAllCreateFunctionsContinuous();

    TPZDarcyFlow *mat = new TPZDarcyFlow(EMatId, gmesh->Dimension());
    mat->SetConstantPermeability(1.0);
    cmesh->InsertMaterialObject(mat);

    TPZManVector<REAL, 1> val2(1, 0.);
    TPZFMatrix<REAL> val1(1, 1, 0.);
    const int neumanntype = 1, diritype = 0;

    val2[0] = 0.;
    cmesh->InsertMaterialObject(mat->CreateBC(mat, -1, neumanntype, val1, val2));
    val2[0] = 1.;
    cmesh->InsertMaterialObject(mat->CreateBC(mat, -2, diritype, val1, val2));
    val2[0] = 3.;
    cmesh->InsertMaterialObject(mat->CreateBC(mat, -3, diritype, val1, val2));

    cmesh->AutoBuild();
    return TPZAutoPointer<TPZCompMesh>(cmesh);
}

// -----------------------------------------------------------------------
// Solver helpers
// -----------------------------------------------------------------------

/**
 * Symmetric storage + Pardiso.
 * SetDirect(ECholesky): configures the analysis for a direct solve.
 *   ECholesky is passed as convention but the argument has no effect on
 *   the actual Pardiso algorithm.
 * SetDefPositive(true): switches Pardiso mtype from -2 (general symmetric)
 *   to 2 (SPD), enabling the Cholesky factorization path. This IS what
 *   controls the factorization mode for symmetric external solvers.
 */
static TPZFMatrix<STATE> SolveWithPardisoSym(TPZAutoPointer<TPZCompMesh> cmesh,
                                             int nthreads) {
    TPZLinearAnalysis an(cmesh, RenumType::ENone);
    TPZSSpStructMatrix<STATE> matsp(cmesh);
    matsp.SetNumThreads(nthreads);
    an.SetStructuralMatrix(matsp);

    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Assemble();

    // Must be called after Assemble (matrix object exists) and before Solve
    // (Decompose reads this flag to set Pardiso mtype=2 vs mtype=-2).
    an.MatrixSolver<STATE>().Matrix()->SetDefPositive(true);

    an.Solve();
    return an.Solution();
}

/**
 * Symmetric storage + MUMPS.
 * Same reasoning as SolveWithPardisoSym: SetDirect triggers direct solve,
 * SetDefPositive(true) switches MUMPS sym from 2 (LDL^T) to 1 (Cholesky).
 */
static TPZFMatrix<STATE> SolveWithMumpsSym(TPZAutoPointer<TPZCompMesh> cmesh,
                                           int nthreads) {
    TPZLinearAnalysis an(cmesh, RenumType::ENone);
    TPZSSpStructMatrixMumps<STATE> matsp(cmesh);
    matsp.SetNumThreads(nthreads);
    an.SetStructuralMatrix(matsp);

    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Assemble();

    an.MatrixSolver<STATE>().Matrix()->SetDefPositive(true);

    an.Solve();
    return an.Solution();
}

/**
 * Non-symmetric storage + Pardiso.
 * SetDirect(ELU): configures direct solve; the ELU argument has no effect
 *   on Pardiso's algorithm (non-symmetric matrices always use mtype=11).
 * SetDefPositive: omitted because it has no effect on non-symmetric matrices
 *   (sym=0 / mtype=11 is always used regardless).
 */
static TPZFMatrix<STATE> SolveWithPardisoNonSym(TPZAutoPointer<TPZCompMesh> cmesh,
                                                int nthreads) {
    TPZLinearAnalysis an(cmesh, RenumType::ENone);
    TPZSpStructMatrix<STATE> matsp(cmesh);
    matsp.SetNumThreads(nthreads);
    an.SetStructuralMatrix(matsp);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an.SetSolver(step);
    an.Assemble();

    an.Solve();
    return an.Solution();
}

/**
 * Non-symmetric storage + MUMPS.
 * Same reasoning as SolveWithPardisoNonSym: non-symmetric matrices always
 * use MUMPS sym=0; SetDefPositive has no effect here.
 */
static TPZFMatrix<STATE> SolveWithMumpsNonSym(TPZAutoPointer<TPZCompMesh> cmesh,
                                              int nthreads) {
    TPZLinearAnalysis an(cmesh, RenumType::ENone);
    TPZSpStructMatrixMumps<STATE> matsp(cmesh);
    matsp.SetNumThreads(nthreads);
    an.SetStructuralMatrix(matsp);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an.SetSolver(step);
    an.Assemble();

    an.Solve();
    return an.Solution();
}

// -----------------------------------------------------------------------
// Helper: relative L2 error between two solution vectors
// -----------------------------------------------------------------------

/// Returns ||a - b||_2 / ||a||_2  (or ||a - b||_2 when ||a||_2 ~= 0).
static STATE RelativeL2Error(const TPZFMatrix<STATE> &ref,
                             const TPZFMatrix<STATE> &other) {
    const int64_t n = ref.Rows();
    STATE normDiff = 0., normRef = 0.;
    for (int64_t i = 0; i < n; i++) {
        const STATE d = ref.GetVal(i, 0) - other.GetVal(i, 0);
        normDiff += d * d;
        normRef  += ref.GetVal(i, 0) * ref.GetVal(i, 0);
    }
    normDiff = std::sqrt(normDiff);
    normRef  = std::sqrt(normRef);
    constexpr STATE smallVal = 1e-10;
    return (normRef > smallVal) ? normDiff / normRef : normDiff;
}

} // namespace solvercomparison

// -----------------------------------------------------------------------
// Macro helper: builds a fresh independent mesh pair and runs a solver
// -----------------------------------------------------------------------

#define SOLVE_WITH(SolverFn, nDivVal, pordVal, nthreadsVal)                 \
    [&]() -> TPZFMatrix<STATE> {                                            \
        using namespace solvercomparison;                                   \
        TPZAutoPointer<TPZGeoMesh> gmesh = CreateGeoMesh3D(nDivVal);        \
        TPZAutoPointer<TPZCompMesh> cmesh = CreateCompMesh(gmesh, pordVal); \
        return SolverFn(cmesh, nthreadsVal);                                \
    }()

// -----------------------------------------------------------------------
// Test cases
// -----------------------------------------------------------------------

TEST_CASE("Symmetric solvers: Mumps vs Pardiso - 3D Darcy flow",
          "[solver_comparison][symmetric]") {
    using namespace solvercomparison;

    const int nDiv    = GENERATE(1, 2, 3);
    const int pord    = GENERATE(1, 2);
    constexpr int nthreads = 4;
    constexpr STATE tol    = 1e-10;

    CAPTURE(nDiv, pord);

    const auto solPardiso = SOLVE_WITH(SolveWithPardisoSym, nDiv, pord, nthreads);
    const auto solMumps   = SOLVE_WITH(SolveWithMumpsSym,   nDiv, pord, nthreads);

    REQUIRE(solPardiso.Rows() == solMumps.Rows());

    const STATE err = RelativeL2Error(solPardiso, solMumps);
    CAPTURE(err);
    REQUIRE(err < tol);
}

TEST_CASE("Non-symmetric solvers: Mumps vs Pardiso - 3D Darcy flow",
          "[solver_comparison][nonsymmetric]") {
    using namespace solvercomparison;

    const int nDiv    = GENERATE(1, 2, 3);
    const int pord    = GENERATE(1, 2);
    constexpr int nthreads = 4;
    constexpr STATE tol    = 1e-10;

    CAPTURE(nDiv, pord);

    const auto solPardiso = SOLVE_WITH(SolveWithPardisoNonSym, nDiv, pord, nthreads);
    const auto solMumps   = SOLVE_WITH(SolveWithMumpsNonSym,   nDiv, pord, nthreads);

    REQUIRE(solPardiso.Rows() == solMumps.Rows());

    const STATE err = RelativeL2Error(solPardiso, solMumps);
    CAPTURE(err);
    REQUIRE(err < tol);
}

TEST_CASE("Symmetric vs non-symmetric storage: same solution - 3D Darcy flow",
          "[solver_comparison][cross]") {
    using namespace solvercomparison;

    const int nDiv    = GENERATE(1, 2, 3);
    const int pord    = GENERATE(1, 2);
    constexpr int nthreads = 4;
    constexpr STATE tol    = 1e-10;

    CAPTURE(nDiv, pord);

    // Pardiso symmetric (reference) vs MUMPS non-symmetric (cross-check)
    const auto solRef   = SOLVE_WITH(SolveWithPardisoSym,   nDiv, pord, nthreads);
    const auto solCross = SOLVE_WITH(SolveWithMumpsNonSym,  nDiv, pord, nthreads);

    REQUIRE(solRef.Rows() == solCross.Rows());

    const STATE err = RelativeL2Error(solRef, solCross);
    CAPTURE(err);
    REQUIRE(err < tol);
}

#undef SOLVE_WITH

#endif // PZ_USING_MKL && PZ_USING_MUMPS
