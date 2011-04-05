#ifndef NONLINANALYSISH
#define NONLINANALYSISH

#include "pzanalysis.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <iostream>
class TPZCompMesh;
class TPZFMatrix;

/** Jorge
 * @brief Derived class from TPZAnalysis implements non linear analysis (Newton's method)
 * @ingroup analysis
 */
class TPZNonLinearAnalysis : public TPZAnalysis {

public:

TPZNonLinearAnalysis(TPZCompMesh *mesh,std::ostream &out);

TPZNonLinearAnalysis();

virtual ~TPZNonLinearAnalysis();

/**
 * It process a Newton's method to solve the non-linear problem.
 * It has the possibility of line search with parameter linesearch = true.
 */
virtual void IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch = false, bool checkconv = false);

/** Implements a golden section line search.
 * Parameter DeltaW must be a copy. Please do not put a &
 * It is because usually here and in derived classes fSolution was passed
 * as DeltaW. But fSolution changes in the linesearch procedure when LoadSolution
 * is called before AssembleResidual.
 */
REAL LineSearch(const TPZFMatrix &Wn, TPZFMatrix DeltaW, TPZFMatrix &NextW, REAL tol, int niter);

REAL SolutionNorm();

virtual void ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &coefs, int icase);

int NumCases();

virtual void Residual(TPZFMatrix &residual, int icase);

virtual void LoadSolution();

virtual void LoadSolution(const TPZFMatrix &state);

void LoadState(TPZFMatrix &state);

};

#endif

