#ifndef PZELAST3D
#define PZELAST3D

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <vector>

/** This class implements an isotropic elasticity material.
 *  Problem dimension is 1, 2 or 3.
 *  @since Aug 31, 2005.
 */
class TPZElasticity3D : public TPZMaterial {

public :

/** Class constructor.
 * @param nummat - material ID.
 * @param E - Young's modulus.
 * @param poisson - poisson's ratio 
 * @param force - external forces
 */ 
TPZElasticity3D(int nummat, REAL E, REAL poisson, TPZVec<REAL> &force);

/** Class destructor.
 */
virtual ~TPZElasticity3D();

/** Returns material dimension.
 */
int Dimension() { return 3;}

/** Number of state variables.
 */
int NStateVariables(){ return 3;}

/** Print material report.
 */
virtual void Print(std::ostream & out);

/** Material name.
 */
char *Name() { return "TPZElasticity3D"; }

/** Contribute to stiff matrix and load vector.
 *  See base class to more informations.
 */
virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
                        TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);

/** Implements Dirichlet and Neumann boundary conditions.
 *  See base class to more informations.
 */
virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
                          TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

/** Returns index of post-processing variable.
 */
virtual int VariableIndex(char *name);

/** Number of data of variable var.
 */
virtual int NSolutionVariables(int var);

/** Post-processing method. Based on solution Sol and its derivatives DSol, it computes the post-processed variable var.
 */
virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);

/** Return the number of components which form the flux function
 * Method not implemented.
 */
virtual int NFluxes() { 
  PZError << "\nTPZElasticity3D::NFluxes() - Method not implemented\n";
  return 0;
}

/** Compute the value of the flux function to be used by ZZ error estimator.
 * Method not implemented.
 */
virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux){
  PZError << "\nTPZElasticity3D::Flux - Method not implemented\n";
}

/** Evaluate error between approximate (FEM) and exact solutions.
 */
virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u, TPZFMatrix &dudx, 
                    TPZFMatrix &axes, TPZVec<REAL> &flux,
                    TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
/**
 * Returns the number of norm errors: 3 (Semi H1, L2 and H1).
 */
virtual int NEvalErrors() {return 3;} 

private :

/** Young's modulus.
 */  
REAL fE;
  
/** Poisson's ratio.
 */  
REAL fPoisson;
  
/** External forces.
 */  
TPZVec<REAL> fForce;

};

#endif
