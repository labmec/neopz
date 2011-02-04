//$Id: pzblackoil2p3d.h,v 1.4 2011-02-04 08:53:02 fortiago Exp $

#ifndef PZBLACKOIL2P3D_H
#define PZBLACKOIL2P3D_H

#include "pzmaterial.h"
#include "pzdiscgal.h"

#ifdef _AUTODIFF

#include "fad.h"


/**
 * Implements a 3D two-phase (oil-water) black-oil flow.
 * @since November 11, 2008
*/
class TPZBlackOil2P3D : public TPZDiscontinuousGalerkin{

public:

typedef TFad<4, REAL> BFadREAL;

protected:

  /** 
   * Interpolacao linear
   */
  void Interpolate(std::map<REAL,REAL> &dados, double x, double &y, double &dy);
  void Interpolate(std::map<REAL,REAL> &dados, BFadREAL x, BFadREAL &y);

  /** Simulation time step
   */
  double fDeltaT;

  /** Estados: un ou un+1
   */
  enum EState { ELastState = 0, ECurrentState = 1 };

  static EState gState;
  
  void testedoBo();
  
  void testeKrw();

public:

  static void SetLastState(){ gState = ELastState; }
  static void SetCurrentState(){ gState = ECurrentState; }

  /** Class constructor 
   * @param id material id
   * @param dim problem dimension
   * @param nstate number of state variables
   * @param sol constant solution vector
   */
  TPZBlackOil2P3D(int id, double deltaT);

  /** Class destructor
   */
  ~TPZBlackOil2P3D();

  /** Copy constructor
   */
  TPZBlackOil2P3D(const TPZBlackOil2P3D &cp);

  /** Defines simulation time step.
   */
  void SetTimeStep(double timestep){
    this->fDeltaT = timestep;
  }

  /** Creates another material of the same type 
   */
  virtual TPZAutoPointer<TPZMaterial> NewMaterial();

  /** Returns problem dimension 
   */
  virtual int Dimension(){ return 3; }

  /** Returns number of state variables: oil pressure and oil saturation
   */
  virtual int NStateVariables(){ return 2; }

  /** Contribute method
   */
  virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

   /**
   * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
   * @param data[in] stores all input data
   * @param weight[in] is the weight of the integration rule
   * @param ek[out] is the stiffness matrix
   * @param ef[out] is the load vector
   * @param bc[in] is the boundary condition material
   */
  virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc);

  /** To satisfy base class interface.
   */
  virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

  /** To satisfy base class interface.
   */
  virtual void ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);


  ///Solution methods

  /** Solution indices of post-processing
    */
  enum ESolutionVars { ENone = 0, EWaterPressure = 1, EOilPressure, EWaterSaturation, EOilSaturation, EDarcyVelocity };

  /** It returns the variable index associated with the name 
   */
  virtual int VariableIndex(const std::string &name);

  /** It returns the number of variables associated with the variable
   *  indexed by var.  
   * var is obtained by calling VariableIndex
   */
  virtual int NSolutionVariables(int var);

  /** It returns the solution associated with the var index based on
   * the finite element approximation
   */
  virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
                        TPZFMatrix &axes, int var, TPZVec<REAL> &Solout);


  /** Fill material data parameter with necessary requirements for the
    * Contribute method. Here, in base class, all requirements are considered
    * as necessary. Each derived class may optimize performance by selecting
    * only the necessary data.
    */
  virtual void FillDataRequirements(TPZMaterialData &data){
    data.SetAllRequirements(true);
    data.fNeedsNeighborSol = false;
    data.fNeedsNeighborCenter = false;
  }

  /** Fill material data parameter with necessary requirements for the
    * ContributeInterface method. Here, in base class, all requirements are considered
    * as necessary. Each derived class may optimize performance by selecting
    * only the necessary data.
    */
  virtual void FillDataRequirementsInterface(TPZMaterialData &data){
    data.SetAllRequirements(true);
    data.fNeedsSol = false;
  }

  ///Dados

  /**
   * Permeabilidade relativa do oleo
   * Kro = Kro( Sw )
   */
  void Kro(double So, double &Kro, double &dKroSo);
  void Kro(BFadREAL So, BFadREAL &Kro);

  /**
   * Permeabilidade relativa da agua
   * Krw = Krw( Sw )
   */
  void Krw(double So, double &Krw, double &dKrwSo);
  void Krw(BFadREAL So, BFadREAL &Krw);

  /** Bo = Bo( po )
   */
  void Bo(double po, double &Bo, double &dBoDpo);
  void Bo(BFadREAL po, BFadREAL &Bo);

  /** Viscosidade do oleo ViscOleo = ViscOleo( po )
   */
  void ViscOleo(double po, double &ViscOleo, double &dViscOleoDpo);
  void ViscOleo(BFadREAL po, BFadREAL &ViscOleo);

  /** Pressao capilar
   * pc = pc( Sw )
   */
  void PressaoCapilar(double So, double &pc, double &DpcDSo);
  void PressaoCapilar(BFadREAL So, BFadREAL &pc);

  /** Porosidade
   * Phi = Phi( pw ) - fizemos como Phi ( po )
  */
  void Porosidade(double po, double &poros, double &dPorosDpo);
  void Porosidade(BFadREAL po, BFadREAL &poros);

  ///Dados constantes

  /** Densidade do oleo em condicoes padroes - kg/m3
   */
  double RhoOleoSC();

  /** Densidade da agua em condicoes padroes - kg/m3
   */
  double RhoAguaSC();

  /** Aceleracao da gravidade
   */
  double g();

  /** Bw = constante
   */
  double Bw();

  /** Viscosidade da agua (constante)
  */
  double ViscAgua();

  /** Permeabilidade absoluta 
   */
  void K(TPZFMatrix &K);

};

#endif

#endif
