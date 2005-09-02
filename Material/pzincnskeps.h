// -*- c++ -*-

/** @file pzincnskeps.h
 *
 * Header file for class TPZIncNSKEps
 */

#ifndef PZINCNSKEPS
#define PZINCNSKEPS

#include <iostream>
#include <string.h>

#include "pzmaterial.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

class TPZBndCond;

/** This class implements an imcompressible Navier-Stokes formulation with modified KEpsilon turbulence model.
 * Variables are: {K, Eps, Pressure, Vx, Vy, Vz}.
 * This class is homework:
 * @author Professor Paulo Vatavuk
 * @author Professor Philippe Devloo
 * @author Edimar Cesar Rylo
 * @author Luis Fernando
 * @author Roberto H. Heiderich 
 * @author Tiago Forti
 * @since June 29, 2005
 */
class  TPZIncNavierStokesKEps : public TPZMaterial {

  private:
  
    int fDimension;
    
    REAL fMU, fRHO, fCmu, fSigmaK, fSigmaEps, fCepsilon1, fCepsilon2;
    
    TPZVec<REAL> fBodyForce; //fc = {0,0,-fc}
    
    /** Dot for matrices with same dimensions. Tr[A B]. No consistence test is made. */
    REAL Dot(TPZFMatrix &A, TPZFMatrix &B);
    
    /** Dot of vector A with row BRow of matrix B. */
    REAL Dot(TPZVec<REAL> &A, TPZFMatrix &B, int BRow);

    /** Dot for vectors with same dimensions. No consistence test is made. */        
    REAL Dot(TPZVec<REAL> &A, TPZVec<REAL> &B);

  public:
  
    enum VARIABLES {ENone = -1, EK = 0, EEpsilon, EPressure, EVx, EVy, EVz, EVvector};

    TPZIncNavierStokesKEps(int id, int dimension);
    
    void SetParameters(REAL MU, REAL RHO, REAL Cmu, REAL SigmaK, REAL SigmaEps, REAL Cepsilon1, REAL Cepsilon2, TPZVec<REAL> &BodyForce );
    
    void GetParameters(REAL &MU, REAL &RHO, REAL &Cmu, REAL &SigmaK, REAL &SigmaEps, REAL &Cepsilon1, REAL &Cepsilon2, TPZVec<REAL> &BodyForce );    

    virtual ~TPZIncNavierStokesKEps();

    virtual int Dimension();

    /** returns the number of state variables associated with the material*/
    virtual int NStateVariables();

    /** print out the data associated with the material*/
    virtual void Print(std::ostream &out = std::cout);

    /** returns the number of variables associated with the variable
     *  indexed by var.  var is obtained by calling VariableIndex*/
    virtual int NSolutionVariables(int var);

    /**returns the solution associated with the var index based on
     * the finite element approximation */
    virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
                          TPZFMatrix &axes, int var, TPZVec<REAL> &Solout);

    /**Compute contribution to the tangent matrix and residual
     * at an integration point*/
    virtual void Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,
                            TPZVec<REAL> &sol, TPZFMatrix &dsol,
                            REAL weight,TPZFMatrix &axes,
                            TPZFMatrix &phi, TPZFMatrix &dphi,
                            TPZFMatrix &ek, TPZFMatrix &ef);
                            
    /**Compute contribution to the residual at an integration point*/
    virtual void Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,
                              TPZVec<REAL> &sol, TPZFMatrix &dsol, REAL weight,
                              TPZFMatrix &axes, TPZFMatrix &phi,
                              TPZFMatrix &dphi, TPZFMatrix &ef);                            


    /** Compute contribution to the stiffness matrix and right hand
      * side at the integration point of a boundary*/
    virtual void ContributeBC(TPZVec<REAL> &x, TPZVec<REAL> &sol,
                              REAL weight, TPZFMatrix &axes,
                              TPZFMatrix &phi, TPZFMatrix &ek,
                              TPZFMatrix &ef, TPZBndCond &bc);


    /**Compute the error due to the difference between the
     * interpolated flux and the flux computed based on the
     * derivative of the solution
     */
    virtual void Errors(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix &dsol,
                        TPZFMatrix &axes, TPZVec<REAL> &flux,
                        TPZVec<REAL> &uexact, TPZFMatrix &duexact,
                        TPZVec<REAL> &val){
        PZError << __PRETTY_FUNCTION__ << std::endl;
        PZError << "Method not implemented! Error comparison not available. Please, implement it." << std::endl;
    }

};


#endif

