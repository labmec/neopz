#ifndef ELASMATHPP
#define ELASMATHPP

#include <iostream>

#include "pzmaterial.h"

/// This class implements a two dimensional elastic material in plane stress or strain
class TPZElasticityMaterial : public TPZMaterial {

public :

  TPZElasticityMaterial();
  /**Creates an elastic material with:
            elasticity modulus  =   E
           poisson coefficient  =   nu
            forcing function -x =   fx
            forcing function -y =   fy
	    fplainstress = 1 indicates use of plainstress
  */
  TPZElasticityMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy, int plainstress = 1);
  
  /**Copies the data of one TPZElasticityMaterial object to 
     another*/
  TPZElasticityMaterial(const TPZElasticityMaterial &copy);

  /**Creates a new material from the current object   ??*/
  virtual TPZAutoPointer<TPZMaterial> NewMaterial() { return new TPZElasticityMaterial(*this);}

  /**Destructor*/
  virtual ~TPZElasticityMaterial();

  /**Returns the dimension*/
  int Dimension() { return 2;}

  /**Returns the number of state variables associated
     with the material*/
 virtual  int NStateVariables();

  /**Print the material data*/
  virtual void Print(std::ostream & out = std::cout);

  /**Returns the material name*/
  char *Name() { return "TPZElasticityMaterial"; }

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

  /**Returns the variable index associated with the name*/
  virtual int VariableIndex(char *name);

  /**Returns the number of variables associated with the variable
	   indexed by var. var is obtained by calling VariableIndex*/
  virtual int NSolutionVariables(int var);

  /**Returns the solution associated with the var index based
	   on the finite element approximation*/
  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);

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

  virtual int ClassId() const;

  virtual void Read(TPZStream &buf, void *context);

  virtual void Write(TPZStream &buf, int withclassid);



private:
  /**Elasticity modulus*/
  REAL fE;

  /**Poison coeficient*/
	REAL fnu;

  /**Forcing vector*/
  REAL ff[3];

  /**G = E/2(1-nu)*/
  REAL fEover21PlusNu;

  /** E/(1-nu)*/
  REAL fEover1MinNu2;  

  /**Pre Stress Tensor - Sigma XX*/
  REAL fPreStressXX;

  /**Pre Stress Tensor - Sigma YY*/
  REAL fPreStressYY;

  /**Pre Stress Tensor - Sigma XY*/
  REAL fPreStressXY;

  /** Uses plain stress*/
  int fPlaneStress;
};

#endif
