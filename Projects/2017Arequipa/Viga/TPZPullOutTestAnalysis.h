//---------------------------------------------------------------------------

#ifndef TPZPullOutTestAnalysisH
#define TPZPullOutTestAnalysisH

#include "pznonlinanalysis.h"
#include "map"
#include "pzfmatrix.h"
#include "pzbndcond.h"

class TPZPullOutTestAnalysis : public TPZNonLinearAnalysis {

private:

  std::map<TPZBndCond*, TPZFMatrix<STATE> > fFullBCData;

  void FillBCData();

  void SetBCValue(REAL subStep);

  std::pair<int,REAL> SolveOneStep(int TimeStep, std::ostream &out,REAL tol,int numiter, bool linesearch);

  void PostProcess(int istep, REAL subStep, const std::string &filename);
  void PostProcessAverageValues();

  void SupportReaction(int istep, int supportMatId, int direction, std::ostream &myfile);

  void NodalDisplacement(int istep, int nodeMatId, std::ostream &myfile);

  void UpdatePlasticDeformation();

  void UpdateDamage(int istep);

  void ForgetDamageStep();

  REAL RhsNorm(const TPZFMatrix<STATE> &origrhs, bool removeBigNumbers, bool infinity) const;

public:

TPZPullOutTestAnalysis();

TPZPullOutTestAnalysis(TPZCompMesh *mesh);

virtual ~TPZPullOutTestAnalysis();

/** Implements a golden section line search.
 * Parameter DeltaW must be a copy. Please do not put a &
 * It is because usually here and in derived classes fSolution was passed
 * as DeltaW. But fSolution changes in the linesearch procedure when LoadSolution
 * is called before AssembleResidual.
 */
virtual REAL LineSearch(const TPZFMatrix<STATE> &Wn, const TPZFMatrix<STATE> DeltaW, REAL tol, int niter, TPZFMatrix<STATE> &NextW);

REAL LineSearch(const TPZFMatrix<STATE> &Wn, const TPZFMatrix<STATE> DeltaW, TPZFMatrix<STATE> &NextW);

void TimeSteps(std::ostream &out,REAL tol,int numiter, bool linesearch, int nsubsteps,
               const std::string &filename);

};
#endif
