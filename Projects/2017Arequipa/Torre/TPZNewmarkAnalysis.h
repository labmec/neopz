#ifndef TPZNewmarkAnalysisH
#define TPZNewmarkAnalysisH

#include "pznonlinanalysis.h"
#include "map"
#include "pzfmatrix.h"
#include "pzbndcond.h"
#include "TPZEulerBernoulliBeamData.h"

class TPZNewmarkAnalysis : public TPZAnalysis {

  private:

  ///Structure to hold solution vectors: u, v, and a
  struct SSolution{
    TPZFMatrix<STATE> fDisplacement, fVelocity, fAcceleration;
    bool HasData() const{
      return fDisplacement.Rows();
    }
  };

  //Pointer to global data
  TPZAutoPointer< TPZEulerBernoulliBeamData > fPropertyData;

  ///Solution at time step n
  SSolution fLastSol;

  ///Solution at time step n+1
  SSolution fNextSol;

  ///Mass and viscous matrices
  TPZAutoPointer< TPZMatrix<STATE> > fMassMatrix, fViscousMatrix, fStiffnessMatrix;

  ///Coefficients of viscous matrix = Mass fViscousMassCoeff + Stiffness fViscousStiffCoef
  REAL fViscousMassCoeff, fViscousStiffCoef;

  ///Newmark's parameters
  REAL fNewmarkBeta, fNewmarkGamma;

  ///Time step
  REAL fDeltaT;

  ///One time step solver
  std::pair<int,REAL> SolveOneStep(int TimeStep, std::ostream &out,REAL tol,int numiter, bool linesearch);

  void PostProcess(int istep, REAL time, const std::string &filename);

  void SupportReaction(int istep, REAL time, int supportId, std::ostream &myfile);

  void NodalDisplacement(int istep, REAL time, int nodeId, std::ostream &myfile);

  //Returns norm in terms of force and moment equilibrium
  std::pair<REAL,REAL> RhsNorm(const TPZFMatrix<STATE> &origrhs) const;

  ///Constructs mass matrix
  void AssembleMassMatrix(const TPZFMatrix<STATE> & u);

  ///Update displacement, velocity and acceleration vectors
  void UpdateNextSol(const TPZFMatrix<STATE> &nextAcceleration);

  //fRhs = F - K u
  void AssembleStiffnessMatrix(const TPZFMatrix<STATE> & u);

  //fRhs = F - K u
  void AssembleElasticForces(const TPZFMatrix<STATE> & u);

  //Build fviscousMatrix. It requires Mass and Stiffness have already been computed.
  void BuildViscousMatrix();

  /** Computes residual and its jacobian. Jacobian is the derivative of R with respect to the acceleration
   */
  void ResidualAndJacobian(bool updateJacobian);

  /** Implements a golden section line search.
   * Parameter DeltaAccel must be a copy. Please do not put a &
   * It is because usually here and in derived classes fSolution was passed
   * as DeltaAccel. But fSolution changes in the linesearch procedure when UpdateNextSol
   * is called before ResidualAndJacobian.
   */
  virtual REAL LineSearch(const TPZFMatrix<STATE> &Acceleration, const TPZFMatrix<STATE> DeltaAccel, REAL tol, int niter, TPZFMatrix<STATE> &NextAcc);

  public:

  ///Default constructor
  TPZNewmarkAnalysis();

  ///Class constructor
  TPZNewmarkAnalysis(TPZCompMesh *mesh, TPZAutoPointer< TPZEulerBernoulliBeamData > PropertyData);

  ///Destructor
  virtual ~TPZNewmarkAnalysis();

  ///Defines initial solution: u, v, a
  void SetInitialSolution( const TPZFMatrix<STATE> & Displacement, const TPZFMatrix<STATE> & Velocity, const TPZFMatrix<STATE> & Acceleration );

  ///Coefficients of viscous matrix = Mass fViscousMassCoeff + Stiffness fViscousStiffCoef
  void SetViscousMatrixCoefficients(REAL ViscousMassCoeff = 0* 1./100., REAL ViscousStiffCoef = 0* 1./100.);

  ///Newmark's parameters
  void SetNewmarkParameters(REAL NewmarkBeta = 0.25, REAL NewmarkGamma = 0.5);

  ///Performs time steps
  void TimeSteps(std::ostream &out,
                 REAL DeltaT,
                 REAL SteadyStateTol, int MaxTimeSteps,
                 REAL NewtonTol, int NewtonMaxIter,
                 const std::string &filename);


};

#endif

