#ifndef TPZFRACDATA_H
#define TPZFRACDATA_H

#include "pzvec.h"
#include "pzfmatrix.h"

/**
 * @author Omar Duran and Nathan Shauer
 * @since 19/08/2014
 * @brief class to store data of fracture propagation simulation
 */
class TPZFracData {
  
public:
  
  
  
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
  
  // --------- Leak off ---------
  /** @brief Carter coefficient */
  REAL fCl;
  
  /** @brief Static pressure */
  REAL fPe;
  
  /** @brief Reference pressure */
  REAL fPref;
  
  /** @brief Spurt loss */
  REAL fvsp;
  
  /** @brief Vl for new element after propagation  */
  REAL fAccumVl;
  
  /** @brief P order of pressure (p) analysis for fracturing simulation */
  int fPorderPressure;
  
  /** @brief P order of flow (Q) analysis for fracturing simulation */
  int fPorderFlow;
  
  /** @brief P order of pressure (p) analysis for Darcy simulation */
  int fPorderDarcyPressure;
  
  /** @brief P order of flow (Q) analysis for Darcy simulation */
  int fPorderDarcyFlow;
  
  /** @brief Number of threads used in assemble */
  int fNThreadsForAssemble;
  
  /** @brief Number propagations allowed */
  int fnpropag;
  
  /** @brief Number of refinements in y */
  int fhyref;
  
  /** @brief size of darcy domain in y */
  REAL fLy;
  
  /** @brief Number of elements in the fracture */
  int fnelFrac;

  /** @brief bool to determine if we are coupling fracture with darcy */
  bool fIsCoupled;

  /** @brief bool to determine if we want to plot VTK. It is way faster without it */
  bool fplotVTK;
    
  /** @brief Refinement for postprocessing */
  int fhrefPost;
  
  /** @brief Size of the elements in the fracture */
  REAL felSize;
  
  /** @brief Derivative of opening with relation to the pressure */
  REAL fdwdp;
  
  /** @brief Last q tip */
  REAL fLastQtip;
  
  /** @brief Name of the postprocess file for fracture vtk */
  std::string fpostProcessFileName;
  
  /** @brief Map used ONLY for debbuging */
  std::map<REAL,REAL> fDebugMap, fDebugMap2;
  
public:
    
  /** @brief Enum for materials Ids */
  enum MadIds {EMatDarcy = 1, EBCBottom = 2, EBCRight = 3, EBCTop = 4, EBCLeft = 5, EMatFrac = 6, EMatInterFrac = 20, EBCInlet = 7, EBCOutlet = 8 , EBCAuxBottom = 10};
  
  /** @brief Returns w based on pfrac */
  REAL GetW(REAL pfrac) const;
  
  /** @brief Returns dwdp based on pfrac */
  REAL GetDwDp() const;
  
  /** @brief Sets dwdp which is constant for the simulation. Can only be called after defining some other atributes */
  void SetDwDp();
  
  /** @brief Sets next time of simulation */
  void SetNextTime();
  
  /** @brief Sets the postprocess file for fracture vtk */
  void SetPostProcessFileName(std::string &postProcessFileName);
  
  /** @brief Returns the postprocess file for fracture vtk */
  std::string PostProcessFileName();
  
  /** @brief Set fluid viscosity. */
  void SetViscosity(REAL mu){this->fmu = mu;}
  
  /** @brief Returns fluid viscosity. */
  REAL Viscosity() const {return this->fmu;}
  
  /** @brief Set fluid density. */
  void SetDensity(REAL rho){this->fRho = rho;}
  
  /** @brief Returns fluid density. */
  REAL Density() const {return this->fRho;}
  
  /** @brief Set absolute permeability. */
  void SetK(TPZFMatrix<STATE> &Kab){ this->fKab = Kab;}
  
  /** @brief Returns absolute permeability. */
  TPZFMatrix<STATE> K() const;
  
  /** @brief Returns absolute permeability inverse. */
  TPZFMatrix<STATE>  Kinv() const;
  
  /** @brief Set rock porosity. */
  void SetPorosity(REAL phi){this->fPhi = phi;}
  
  /** @brief Returns rock porosity. */
  REAL Porosity() const {return this->fPhi;}
  
  /**
   * @brief \f$ Rock porosity function. \f$ Phi = Phi( p ) \f$
   * @param p pressure
   */
  void Porosity(REAL p, REAL &porosity, REAL &dPorosityDp) const;
  
  /**
   * @brief \f$ Fluid density function. \f$ RhoFluid = RhoFluid( p ) \f$
   * @param p pressure
   */
  void Density(REAL p, REAL &RhoFluid, REAL &dRhoDp) const;
  
  /**
   * @brief Fluid viscosity function. \f$ FluidViscosity = Visc( p ) \f$
   * @param p pressure
   */
  void Viscosity(REAL p, REAL &FluidViscosity, REAL &dFluidViscosityDp) const;
  
  
  /** @brief Set fluid viscosity. */
  void SetTotalTime(REAL Ttot){this->fTtot = Ttot;}
  
  /** @brief Returns fluid viscosity. */
  REAL TotalTime() const {return this->fTtot;}
  
  /** @brief Defines simulation time step. */
  void SetTimeStep(REAL timestep){ this->fDeltaT = timestep;}
  
  /** @brief Returns simulation time step. */
  REAL TimeStep() const {return this->fDeltaT;}
  
  /** @brief Defines simulation time */
  void SetTime(REAL time){ this->fTime = time;}
  
  /** @brief Returns simulation time */
  REAL Time() const {return this->fTime;}
  
  /** @brief Defines simulation temporal scheme  */
  void SetTheta(REAL theta){ this->fTheta = theta;}
  
  /** @brief Returns simulation temporal scheme  */
  REAL Theta() const {return this->fTheta;}
  
  /** @brief Defines simulation fracture length */
  void SetLfrac(REAL Lfrac){ this->fLfrac = Lfrac;}
  
  /** @brief Returns simulation fracture length */
  REAL Lfrac() const {return this->fLfrac;}
  
  /** @brief Defines simulation fracture height */
  void SetHf(REAL Hf){ this->fHf = Hf;}
  
  /** @brief Returns simulation fracture height */
  REAL Hf() const {return this->fHf;}
  
  /** @brief Defines simulation elasticity */
  void SetE(REAL E){ this->fE = E;}
  
  /** @brief Returns simulation elasticity */
  REAL E() const {return this->fE;}
  
  /** @brief Defines simulation poisson */
  void SetPoisson(REAL nu){ this->fnu = nu;}
  
  /** @brief Returns simulation poisson */
  REAL Poisson() const {return this->fnu;}
  
  /** @brief Defines simulation Flow BC */
  void SetQ(REAL Q){ this->fQ = Q;}
  
  /** @brief Returns simulation Flow BC */
  REAL Q() const {return this->fQ;}
  
  /** @brief Defines simulation Confinement stress */
  void SetSigmaConf(REAL SigmaConf){ this->fSigmaConf = SigmaConf;}
  
  /** @brief Returns simulation confinement stress */
  REAL SigmaConf() const {return this->fSigmaConf;}
  
  /** @brief Returns simulation G */
  REAL G() const;

  /** @brief Defines if coupled simulation */
  void SetIsCoupled(bool isCoupled){ this->fIsCoupled = isCoupled;}
  
  /** @brief Returns if using coupled simulation */
  bool IsCoupled() const {return this->fIsCoupled;}

  /** @brief Defines if plot vtk */
  void SetIsPlotVTK(bool plotVTK){ this->fplotVTK = plotVTK;}

  /** @brief Returns if ploting vtk */
  bool IsPlotVTK() const {return this->fplotVTK;}
    
  /** @brief Defines if plot vtk */
  void SetHrefPostPro(int hrefpost){ this->fhrefPost = hrefpost;}
    
  /** @brief Returns if ploting vtk */
  int HrefPostPro() const {return this->fhrefPost;}
    
  /** @brief Defines p order of the pressure in L2 space */
  void SetPorderPressure(int PorderPressure){ this->fPorderPressure = PorderPressure;}
  
  /** @brief Returns p order of the pressure */
  int PorderPressure() const {return this->fPorderPressure;}
  
  /** @brief Defines p order of the flow in H1 space */
  void SetPorderFlow(int PorderFlow){ this->fPorderFlow = PorderFlow;}
  
  /** @brief Returns p order of the flow */
  int PorderFlow() const {return this->fPorderFlow;}
  
  
  /** @brief Defines p order of the flow in hdiv space */
  void SetPorderDarcyFlow(int PorderDarcyFlow){ this->fPorderDarcyFlow = PorderDarcyFlow;}
  
  /** @brief Returns p order of the Darcy flow */
  int PorderDarcyFlow() const {return this->fPorderDarcyFlow;}
  
  /** @brief Defines p order of the Darcy pressure in L2 space */
  void SetPorderDarcyPressure(int PorderDarcyPressure){ this->fPorderDarcyPressure = PorderDarcyPressure;}
  
  /** @brief Returns p order of the Darcy pressure */
  int PorderDarcyPressure(){return this->fPorderDarcyPressure;}
  
  
  /** @brief Defines Number of elements in the fracture */
  void SetNelFrac(int NelFrac){ this->fnelFrac = NelFrac;}
  
  /** @brief Returns number of elements in the fracture */
  int NelFrac() const {return this->fnelFrac;}
  
  /** @brief Defines p order of the flow in H1 space */
  void SetElSize(REAL elSize){ this->felSize = elSize;}
  
  /** @brief Returns p order of the flow */
  REAL ElSize() const {return this->felSize;}
  
  /** @brief Set evaluating step n */
  void SetLastState(){ State = n;}
  
  /** @brief Set evaluating step n + 1 */
  void SetCurrentState(){ State = nplusone;}
  
  /** @brief Returns true if is last state */
  bool IsLastState() const { return State == n;}
  
  /** @brief Sets Carter coefficient */
  void SetCl(REAL Cl) {fCl = Cl;}
  
  /** @brief Return Carter coefficient */
  REAL Cl() const {return fCl;}
  
  /** @brief Sets static pressure */
  void SetPe(REAL Pe) {fPe = Pe;}
  
  /** @brief Returns static pressure */
  REAL Pe() const {return fPe;}
  
  /** @brief Sets reference pressure for leak off */
  void SetPref(REAL Pref){fPref = Pref;}
  
  /** @brief Returns reference pressure for leak off */
  REAL Pref() const {return fPref;}
  
  /** @brief Sets spurt loss */
  void SetVsp(REAL vsp){fvsp = vsp;}
  
  /** @brief Returns spurt loss */
  REAL Vsp() const {return fvsp;}
  
  /** @brief Sets vl of the element to be propagated */
  void SetAccumVl(REAL AccumVl){fAccumVl = AccumVl;}
  
  /** @brief Returns vl of the element to be propagated */
  REAL AccumVl() const {return fAccumVl;}
  
  /** @brief Sets number of threads used for assemble */
  void SetNThreadsForAssemble(int nthreads) { fNThreadsForAssemble = nthreads; } const
  
  /** @brief Returns number of threads used for assemble */
  int NThreadsForAssemble() { return fNThreadsForAssemble; } const

  /** @brief Sets maximum number of propagations allowerd */
  void SetNPropagations(int npropag) { fnpropag = npropag; } const
  
  /** @brief Returns number of threads used for assemble */
  int NPropagations() { return fnpropag; } const

  /** @brief Sets number of refinements in y */
  void SetHy(int Hy) { fhyref = Hy; } const
  
  /** @brief Returns number of refinements in y */
  int Hy() { return fhyref; } const
  
  /** @brief Sets last q tip */
  void SetLastQtip(REAL lastQtip){fLastQtip = lastQtip;}
  
  /** @brief Returns last q tip */
  REAL LastQtip() const {return fLastQtip;}

  /** @brief Sets size of domain in y */
  void SetLy(REAL Ly) {fLy = Ly;}
  
  /** @brief Returns size of domain in y */
  REAL Ly() const {return fLy;}
  
  /** @brief Return vl based on exposition time and pressure */
  REAL VlFtau(REAL pfrac, REAL tau, REAL pe = -1.) const;
  
  /** @brief Return the "ficticious time" based on acumulated volume and pressure */
  REAL FictitiousTime(REAL VlAcum, REAL pfrac, REAL pe = -1.) const;

  /** @brief Corrects Cl based on the pe and pfrac */
  REAL ClCorrected(REAL pfrac, REAL pe) const;
  
  /** @brief Return the seapage velocity ql (mm/s) */
  REAL QlFVl(REAL VlAcum, REAL pfrac, REAL pe = -1.) const;
  
  /** @brief Return the seapage velocity derivative */
  REAL dQlFVl(REAL VlAcum, REAL pfrac, REAL pe = -1.) const;
  
  /** @brief Return seapage velocity derivative for pore pressure */
  REAL dQlFVlPoros(REAL VlAcum, REAL pfrac, REAL pe = -1.) const;
  
  /** @brief Return debug map */
  std::map<REAL,REAL> & DebugMap() {return fDebugMap;}
  
  /** @brief Prints debug map in Mathematica style */
  void PrintDebugMapForMathematica(std::string filename);
  
  /** @brief Return debug map */
  std::map<REAL,REAL> & DebugMap2() {return fDebugMap2;}
  
};

#endif