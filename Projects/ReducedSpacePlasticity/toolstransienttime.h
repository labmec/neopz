//
//  toolstransienttime.h
//  PZ
//
//  Created by Agnaldo Farias on 9/5/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef PZ_toolstransienttime_h
#define PZ_toolstransienttime_h

#include "pzcmesh.h"
#include "pzanalysis.h"
#include "tpzcompmeshreferred.h"
#include "TPZPlasticFrac2D.h"
#include "TPZSandlerDimaggio.h"
#include "TPZCohesiveBC.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"


class TPZElastoPlasticAnalysis;
class TPZNonLinearAnalysis;


class ToolsTransient
{
public:
  
  ///Default Constructor - Should never be called
  ToolsTransient();
  
  /// Constructor with polinomial order of approximation
  ToolsTransient(int pOrder);
  
  /// Destructor
  ~ToolsTransient();
  
  /// Main Methods
  //---------------------------------------------------------------
  
  /// Main of the program. Method to be called to run everything
  void Run();
  
  /// Initialization of the problem methods
  //---------------------------------------------------------------
  
  /// Creates all the uncoupled parametes and values. Cmesh for elastoplastic with ReducedSpace and cmesh for pressure with H1
  void InitializeUncoupledMeshesAttributes();

  /// Creates the geometric mesh for all cmeshes, and create the solution vectors of the reduced space
  TPZCompMesh * ElastCMeshReferenceProcessed();
  
  /// Creates the geometric mesh. Called in ElastCMeshReferenceProcessed()
  void Mesh2D();

  /// Creates the elastic computational mesh for the reduced space. It uses nloadcases
  TPZCompMesh * CMeshElastic();

  /// Creates the Cmesh using reduces space for coupled simulation
  TPZCompMeshReferred * CMeshReduced(TPZCompMesh * cmeshref);

  /// Creates the computational mesh for pressure. 1 dimensional
  TPZCompMesh * CMeshPressure();
  
  /// Creates the multphysic computational mesh for coupled simulation. fmphysics
  void CMeshMultiphysics();
  
  /// Used to construct cmesh for hat function used in the reduces space. Deprecated! Now using same CmeshElastic() with nloadcases
  TPZCompMesh * CMeshHat(int &dirid, int &porder);
  
  /// Set Methods
  //---------------------------------------------------------------
  
  /// Sets the parameters of the MohrCoulomb model
  void SetMohrCoulombParameters(REAL poisson, REAL elast, REAL cohesion, REAL Phi, REAL Psi);
  
  /// Sets The Sigma in determined stripe of the compmesh
  void SetSigmaNStripeNum(TPZCompMesh * cmeshref, int actStripe);
  
  /// Run Methods
  //---------------------------------------------------------------
  
  /// Solves an elasticity problem given an analysis and a cmesh
  void SolveInitialElasticity(TPZAnalysis &an, TPZCompMesh *Cmesh);
  
  /// Solves the transient coupled problem until it propagates or reach the final time
  bool SolveSistTransient(TPZAnalysis *an, bool initialElasticKickIsNeeded);
  
  /// Assemble the last step
  void MassMatrix(TPZFMatrix<REAL> & Un);
  
  /// Creates the jacobian and residuum for one step of newton method
  void StiffMatrixLoadVec(TPZAnalysis *an,
                          TPZAutoPointer< TPZMatrix<REAL> > & matK1, TPZFMatrix<REAL> &fvec);

  /// Apply equation filter on first solution if using prestress. NOT USED ANYMORE
  void ApplyEquationFilter(TPZAnalysis * an);
  
  
  /// Post Process Methods
  //---------------------------------------------------------------
  
  /// Transfer Solution between meshes. Elastoplastic will change????
  void TransferSolutions(TPZCompMesh * lastMPhysicsCMesh, TPZCompMesh * lastElastReferredCMesh);
  
  /// Transfer Solution between meshes. Elastoplastic will change????
  void TransferElasticSolution(TPZCompMesh * cmeshFrom);
  
  /// Integrates the solution
  REAL IntegrateSolution(TPZCompMesh * cmesh, int variable);//0 = meshvec[0] ; 1 = meshvec[1]
  
  /// Transfer the leakoff between meshes
  void TransferLeakoff(TPZCompMesh * oldMphysicsCMesh);
  
  /// Post Process pressure
  void PostprocessPressure();

  /// Post Process acumulated volume in oppening
  void PostProcessAcumVolW();
  
  /// Post Process leak off volume
  void PostProcessVolLeakoff();

  /// Print Methods
  //---------------------------------------------------------------
  
  /// Method to show the displacement x sigmay of the first cohesive elements
	void ShowDisplacementSigmaYCohesive(TPZCompMesh *cmesh);

  /// Method to show the displacement x sigmay of all the bottom el
	void ShowDisplacementSigmaYBottom(TPZCompMesh *cmesh);

  /// Plot all the hat functions separately using cmeshhat. Deprecated!
  void PlotAllHatsVTK();
  
  /// Tests
  //---------------------------------------------------------------

  /// Apply equation filter on pressure. Used for tests
	void ApplyEquationFilterOnPressure(TPZAnalysis * an);
  
  /// ChechConv of the coupled problem. Uses Equation filter and prestress. Depecrated, now I uses net pressure, so no EqFilter
  void CheckConv(bool OptimizeBandwidth);
  
  /// Main to test plasticity uncoupled. Just run a problem with constante pressure
  void RunPlasticity();
  
  /// Creates elastoplastic mesh for simple uncoupled problem. Used for early tests
  TPZCompMesh * CMeshElastoPlastic(TPZGeoMesh *gmesh, REAL SigmaN);
  
  /// Runs a simulation with uncoupled elastoplasticity given an elastoplasticanalysis and a cmesh
  void SolveInitialElastoPlasticity(TPZElastoPlasticAnalysis &analysis, TPZCompMesh *Cmesh);
  
  /// Creates an H1 mesh for uncoupled simulation. Used to test Cohesive material
	TPZCompMesh* CMeshCohesiveH1(REAL pressure);
  
  /// Solves an NonLinear elasticity problem given an NonLinearAnalysis and a cmesh. Used to test Cohesive
	void SolveNLElasticity(TPZCompMesh *cmesh, TPZNonLinearAnalysis &an);

  
  /// Atributes
  //---------------------------------------------------------------
  
  /// Polinomyal order of approximation
  int fpOrder;
  
  /// bool that indicates when the max time of simulation is reached
  bool fMustStop;
  
  /// Sandler Dimaggio material for coupled problem (Good idea to exchange to pzsandlerextend!)
  //TPZPlasticFrac2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> > * fCouplingMaterial1;
  
  /// MohrCoulomb material for coupled problem
  TPZPlasticFrac2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem> * fCouplingMaterial1;
  TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> fPlasticStepPV;

  /// Cohesive Material
	TPZCohesiveBC * fCohesiveMaterial;
  TPZCohesiveBC * fCohesiveMaterialFirst;

	/// Geometric mesh
  TPZGeoMesh * fgmesh;
 
  /** fmeshvec[0] = Malha computacional elastica do tipo referred */
  /** fmeshvec[1] = Malha computacional de fluxo 1D */
  TPZVec<TPZCompMesh *> fmeshvec;
  TPZCompMesh * fmphysics;
};



template<class TVar>
class TElastSolFunction : public TPZFunction<TVar>
{
public:
  
  /**
   * Class constructor
   */
  TElastSolFunction();
  
  TElastSolFunction(TPZCompMesh * cmesh);
  
  /**
   * Class destructor
   */
  ~TElastSolFunction();
  
  virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f);
  virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df);
  
  /** Returns number of functions.
   */
  virtual int NFunctions();
  
  /** Polynomial order of this function. In case of non-polynomial
   * function it can be a reasonable approximation order.
   */
  virtual int PolynomialOrder();
  
  TPZCompMesh * fcmesh;
  long fIniElIndex;
  
};

template<class TVar>
class TLeakoffFunction : public TPZFunction<TVar>
{
public:
  
  /**
   * Class constructor
   */
  TLeakoffFunction();
  
  TLeakoffFunction(TPZCompMesh * cmesh);
  
  /**
   * Class destructor
   */
  ~TLeakoffFunction();
  
  virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f);
  virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df);
  
  /** Returns number of functions.
   */
  virtual int NFunctions();
  
  /** Polynomial order of this function. In case of non-polynomial
   * function it can be a reasonable approximation order.
   */
  virtual int PolynomialOrder();
  
  TPZCompMesh * fcmesh;
  long fIniElIndex;
};

#endif
