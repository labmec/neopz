#ifndef ELASAXIMATHPP
#define ELASAXIMATHPP

#include <iostream>

#include "pzmaterial.h"
#include "pzdiscgal.h"

#include <vector>
#include <math.h>

const REAL Pi = 4.*atan(1);

/// This class implements a two dimensional elastic material in plane stress or strain
class TPZElasticityAxiMaterial : public TPZDiscontinuousGalerkin {

public :

  TPZElasticityAxiMaterial();
  /**Creates an elastic material with:
            elasticity modulus  =   E
           poisson coefficient  =   nu
            forcing function -x =   fx
            forcing function -y =   fy
	    fplainstress = 1 indicates use of plainstress
  */
  TPZElasticityAxiMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy);

//------------------- FEITO POR AGNALDO : 05/02/10 . Só para teste ----------
/**Creates an elastic material test with:
	elasticity modulus  =   E
	poisson coefficient  =   nu
	forcing function -x =   fx
	forcing function -y =   fy
	symmetrical method = coefTheta (-1-> symmetric, 1->nonsymmetric)
	penalty term = coefAlpha
 */
TPZElasticityAxiMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy, REAL coefTheta,  REAL coefAlpha);

    /**Copies the data of one TPZElasticityAxiMaterial object to
     another*/
    TPZElasticityAxiMaterial(const TPZElasticityAxiMaterial &copy);
    
    /**Creates a new material from the current object   ??*/
    virtual TPZAutoPointer<TPZMaterial> NewMaterial() { return new TPZElasticityAxiMaterial(*this);}
    
    /**Destructor*/
    virtual ~TPZElasticityAxiMaterial();
    

//-----------------------------------------------------------------------------------
	
 /** Set the origin of Revolution Axis (Z),
     the direction of Revolution Axis (Z),
     and the Radius vector (orthogonal with respect of Z axis)
 */
 void SetOrigin(TPZManVector<REAL> &Orig, TPZManVector<REAL> &AxisZ, TPZManVector<REAL> &AxisR);

 REAL ComputeR(TPZVec<REAL> &x);

 void SetMohrCoulomb(double c, double phi)
 {
    f_c = c;
    f_phi = phi;
 }
    
    void SetThermalExpansionCoefficient(REAL alpha)
    {
        fAlpha = alpha;
    }
    
    void SetTemperature(REAL delt)
    {
        fDelTemperature = delt;
    }

  /**Returns the dimension*/
  int Dimension() { return 2;}

  /**Returns the number of state variables associated
     with the material*/
 virtual  int NStateVariables();

  /**Print the material data*/
  virtual void Print(std::ostream & out = std::cout);

  /**Returns the material name*/
  std::string Name() { return "TPZElasticityAxiMaterial"; }

  /**Return the number of components which form the flux function*/
  virtual short NumberOfFluxes(){return 3;}

  /**Return the number of components which form the flux function*/
  virtual int NFluxes(){ return 3;}

  /**Cria as condicoes de contorno*/
  //virtual TPZBndCond *CreateBc(long num, int typ, TPZFMatrix &val1, TPZFMatrix &val2);

  /**Calculates the element stiffness matrix*/
  virtual void Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix &ek,TPZFMatrix &ef);

  //*Applies the element boundary conditions*/
  virtual void ContributeBC(TPZMaterialData &data,REAL weight,
			    TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

  virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef);

  virtual void ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
	  return;
    PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
  }

  /**Returns the variable index associated with the name*/
  virtual int VariableIndex(const std::string &name);

  /**Returns the number of variables associated with the variable
	   indexed by var. var is obtained by calling VariableIndex*/
  virtual int NSolutionVariables(int var);

  /**returns the solution associated with the var index based on
   * the finite element approximation*/
  virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);

  /**Compute the value of the flux function to be used
     by ZZ error estimator*/
  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);

  /**Compute the error due to the difference between
     the interpolated flux and the flux computed based
     on the derivative of the solution*/
  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				  TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
		      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);//Cedric

  //virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix &dsol, TPZFMatrix &axes, TPZVec<REAL> &flux,
  //		      TPZVec<REAL> &uexact,TPZFMatrix &duexact,TPZVec<REAL> &val){}

  /**Returns the elasticity modulus E*/
  REAL E() {return fE;}

  /**Returns the poison coefficient modulus E*/
  REAL Nu() {return fnu;}

  /**SetPresStress Tensor*/
  void SetPreStress(REAL Sigxx, REAL Sigyy, REAL Sigxy);
    
    void SetTemperatureFunction(void (*func)(const TPZVec<REAL> &rz, REAL &temperature))
    {
        fTemperatureFunction = func;
    }

  virtual int ClassId() const;

  virtual void Read(TPZStream &buf, void *context);

  virtual void Write(TPZStream &buf, int withclassid);

  TPZManVector<REAL> GetAxisR();
  TPZManVector<REAL> GetAxisZ();
  TPZManVector<REAL> GetOrigin();

	REAL fIntegral;

private:

  /// Mohr Coulomb Plasticity Criteria Data 
  REAL f_phi;
    /// Mohr Coulomb Plasticity Criteria Data 
  REAL f_c;

  /// Elasticity modulus
  REAL fE;

  /// Poison coeficient
  REAL fnu;

    /// Thermal expansion coeficient
    REAL fAlpha;
    
  /// Forcing vector
  REAL ff[3];

    /// Temperature difference
    REAL fDelTemperature;
    
  /**G = E/2(1-nu)*/
  REAL fEover21PlusNu;

  /** E/(1-nu)*/
  REAL fEover1MinNu2;

  /**Direction of Surface*/
  TPZManVector<REAL> f_AxisR;

  /**Revolution Axis*/
  TPZManVector<REAL> f_AxisZ;

  /**Origin of AxisR and AxisZ*/
  TPZManVector<REAL> f_Origin;
	
 /**symmetric*/ 
 REAL fSymmetric;

 /**penalty term*/
 REAL fPenalty;
    
    ///Function which defines the temperature
    void (*fTemperatureFunction)(const TPZVec<REAL> &rz, REAL &temperature);
	
};

#endif
