#ifndef NONLINANALYSISH
#define NONLINANALYSISH

#include "pzanalysis.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <iostream>
class TPZCompMesh;
class TPZFMatrix;

class TPZNonLinearAnalysis : public TPZAnalysis {

public:

TPZNonLinearAnalysis(TPZCompMesh *mesh,std::ostream &out);

TPZNonLinearAnalysis();

virtual ~TPZNonLinearAnalysis(void);

void IterativeProcess(std::ostream &out,REAL tol,int numiter);

REAL SolutionNorm();

void ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &coefs, int icase);

int NumCases();

void Residual(TPZFMatrix &residual, int icase);

void LoadSolution(TPZFMatrix &state);

void LoadState(TPZFMatrix &state);

};

#endif

