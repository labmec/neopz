#ifndef TPZFRACDATA_H
#define TPZFRACDATA_H

#include "pzvec.h"
#include "pzfmatrix.h"

/**
 * @author Omar Duran and Nathan Shauer
 * @since 19/08/2014
 * @brief Material to store data of fracture propagation simulation
 */
class TPZFracData {
  
public:
  
  /** @brief Default Constructor */
  TPZFracData();
  
  /** @brief Destructor */
  ~TPZFracData();
  
  /** @brief copy constructor */
  TPZFracData(const TPZFracData &copy);
  
  /** @brief operator equal */
  TPZFracData &operator=(const TPZFracData &copy);
  
private:
  
  /** @brief State: Stiffness or Mass Matrix Calculations */
  enum nState { n = 0, nplusone = 1 };
  nState State;
  
  /** @brief Fluid Viscosity - Pa.s */
  REAL fmu;
  
  /** @brief Fluid Density - kg/m3 */
  REAL fRho;
  
  /** @brief Intrinsic absolute permeability - m2 */
  TPZFMatrix<STATE> fKab;
  
  /** @brief Rock porosity - fraction */
  REAL fPhi;
  
  /** @brief Simulation time step */
  REAL fDeltaT;
  
  /** @brief Simulation current time */
  REAL fTime;
  
  /** @brief Simulation Total time */
  REAL fTtot;
  
  /** @brief Simulation temporal scheme (theta = 1, means full implicit) */
  REAL fTheta;
  
  /** @brief Fracture length */
  REAL fLfrac;
  
  /** @brief Fracture height */
  REAL fHf;
  
  /** @brief Reservoir young modulus */
  REAL fE;
  
  /** @brief Reservoir poisson */
  REAL fnu;
  
  /** @brief Flow BC */
  REAL fQ;
  
  /** @brief Confinement Stress */
  REAL fSigmaConf;
  
  /** @brief P order of pressure (p) analysis for fracturing simulation */
  int fPorderPressure;
  
  /** @brief P order of flow (Q) analysis for fracturing simulation */
  int fPorderFlow;
  
  /** @brief Number of elements in the fracture */
  int fnelFrac;
  
  /** @brief Name of the postprocess file for fracture vtk */
  std::string fpostProcessFileName;
  
public:
  
  /** @brief Sets next time of simulation */
  void SetNextTime();
  
  /** @brief Sets the postprocess file for fracture vtk */
  void SetPostProcessFileName(std::string &postProcessFileName);
  
  /** @brief Returns the postprocess file for fracture vtk */
  std::string PostProcessFileName();
  
  /** @brief Set fluid viscosity. */
  void SetViscosity(REAL mu){this->fmu = mu;}
  
  /** @brief Returns fluid viscosity. */
  REAL Viscosity(){return this->fmu;}
  
  // Nathan: Our indentation setting are different
  
  /** @brief Set fluid density. */
  void SetDensity(REAL rho){this->fRho = rho;}
  
  /** @brief Returns fluid density. */
  REAL Density(){return this->fRho;}
  
  /** @brief Set absolute permeability. */
  void SetK(TPZFMatrix<STATE> &Kab){ this->fKab = Kab;}
  
  /** @brief Returns absolute permeability. */
  TPZFMatrix<STATE> K();
  
  /** @brief Returns absolute permeability inverse. */
  TPZFMatrix<STATE>  Kinv();
  
  /** @brief Set rock porosity. */
  void SetPorosity(REAL phi){this->fPhi = phi;}
  
  /** @brief Returns rock porosity. */
  REAL Porosity(){return this->fPhi;}
  
  /**
   * @brief \f$ Rock porosity function. \f$ Phi = Phi( p ) \f$
   * @param p pressure
   */
  void Porosity(REAL p, REAL &porosity, REAL &dPorosityDp);
  
  /**
   * @brief \f$ Fluid density function. \f$ RhoFluid = RhoFluid( p ) \f$
   * @param p pressure
   */
  void Density(REAL p, REAL &RhoFluid, REAL &dRhoDp);
  
  /**
   * @brief Fluid viscosity function. \f$ FluidViscosity = Visc( p ) \f$
   * @param p pressure
   */
  void Viscosity(REAL p, REAL &FluidViscosity, REAL &dFluidViscosityDp);
  
  
  /** @brief Set fluid viscosity. */
  void SetTotalTime(REAL Ttot){this->fTtot = Ttot;}
  
  /** @brief Returns fluid viscosity. */
  REAL TotalTime(){return this->fTtot;}
  
  /** @brief Defines simulation time step. */
  void SetTimeStep(REAL timestep){ this->fDeltaT = timestep;}
  
  /** @brief Returns simulation time step. */
  REAL TimeStep(){return this->fDeltaT;}
  
  /** @brief Defines simulation time */
  void SetTime(REAL time){ this->fTime = time;}
  
  /** @brief Returns simulation time */
  REAL Time(){return this->fTime;}
  
  /** @brief Defines simulation temporal scheme  */
  void SetTheta(REAL theta){ this->fTheta = theta;}
  
  /** @brief Returns simulation temporal scheme  */
  REAL Theta(){return this->fTheta;}
  
  /** @brief Defines simulation fracture length */
  void SetLfrac(REAL Lfrac){ this->fLfrac = Lfrac;}
  
  /** @brief Returns simulation fracture length */
  REAL Lfrac(){return this->fLfrac;}
  
  /** @brief Defines simulation fracture height */
  void SetHf(REAL Hf){ this->fHf = Hf;}
  
  /** @brief Returns simulation fracture height */
  REAL Hf(){return this->fHf;}
  
  /** @brief Defines simulation elasticity */
  void SetE(REAL E){ this->fE = E;}
  
  /** @brief Returns simulation elasticity */
  REAL E(){return this->fE;}
  
  /** @brief Defines simulation poisson */
  void SetPoisson(REAL nu){ this->fnu = nu;}
  
  /** @brief Returns simulation poisson */
  REAL Poisson(){return this->fnu;}
  
  /** @brief Defines simulation Flow BC */
  void SetQ(REAL Q){ this->fQ = Q;}
  
  /** @brief Returns simulation Flow BC */
  REAL Q(){return this->fQ;}
  
  /** @brief Defines simulation Confinement stress */
  void SetSigmaConf(REAL SigmaConf){ this->fSigmaConf = SigmaConf;}
  
  /** @brief Returns simulation confinement stress */
  REAL SigmaConf(){return this->fSigmaConf;}
  
  /** @brief Returns simulation G */
  REAL G();
  
  /** @brief Defines p order of the pressure in L2 space */
  void SetPorderPressure(int PorderPressure){ this->fPorderPressure = PorderPressure;}
  
  /** @brief Returns p order of the pressure */
  int PorderPressure(){return this->fPorderPressure;}
  
  /** @brief Defines p order of the flow in H1 space */
  void SetPorderFlow(int PorderFlow){ this->fPorderFlow = PorderFlow;}
  
  /** @brief Returns p order of the flow */
  int PorderFlow(){return this->fPorderFlow;}
  
  /** @brief Defines p order of the flow in H1 space */
  void SetNelFrac(int NelFrac){ this->fnelFrac = NelFrac;}
  
  /** @brief Returns p order of the flow */
  int NelFrac(){return this->fnelFrac;}
  
  /** @brief Set evaluating step n */
  void SetLastState(){ State = n;}
  
  /** @brief Set evaluating step n + 1 */
  void SetCurrentState(){ State = nplusone;}
  
  bool IsLastState() { return State == n;}
  
};

#endif