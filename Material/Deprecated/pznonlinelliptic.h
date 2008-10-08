//
// C++ Interface: pznonlinelliptic
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZNONLINELLIPTIC_H
#define PZNONLINELLIPTIC_H

#include <pzmaterial.h>

/**
 * Implements a non linear elliptic equation.
 * Laplac( u ) + div(V u) + Sigma u = u^2
 * @author Roberto Heiderich.
 * @since Feb 23, 2005
 */
class TPZNonLinElliptic : public TPZMaterial
{
public:

    /** Constructor.
     * @param id = material identification
     * @param dimension = set problem dimension
     */
    TPZNonLinElliptic(int id, int dimension);

    /** Defaul destructor.
     */
    ~TPZNonLinElliptic();
    
    /**returns the integrable dimension of the material*/
    virtual int Dimension(){ return fDim;}
    
    /** returns the number of state variables associated with the material*/
    virtual int NStateVariables(){return 1;}
    
    /**Compute contribution to the stiffness matrix and right hand
     * side at an integration point*/
    virtual void Contribute(TPZMaterialData &data,
                              REAL weight,
                              TPZFMatrix &ek,
                              TPZFMatrix &ef);

    virtual void Contribute(TPZMaterialData &data,
                              REAL weight,
                              TPZFMatrix &ef);

    /** Compute contribution to the stiffness matrix and right hand
     * side at the integration point of a boundary*/
    virtual void ContributeBC(TPZMaterialData &data,
                                REAL weight, 
                                TPZFMatrix &ek,
                                TPZFMatrix &ef,
                                TPZBndCond &bc);

    void SetParameters(REAL D, TPZVec<REAL> &V, REAL Sigma, REAL LambdaDivK, REAL F);
    
    void GetParameters(REAL &D, TPZVec<REAL> &V, REAL &Sigma, REAL &LambdaDivK, REAL &F);    
   
  virtual int VariableIndex(const std::string &name);

  virtual int NSolutionVariables(int var);

  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &/*axes*/,int var,TPZVec<REAL> &Solout);
  
  virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
			       TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
			       TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);

private:

    /** Problem dimension.
     */
    int fDim;
    
    /** Coeficiente de dispersao.
     * Dispersao populacional.
     * Multiplies laplacian operator.
     */
    REAL fCoeffD;
    
    /** Vetor de conveccao.
     * 
     * Inside divergent operator.
     */
    TPZVec<REAL> fConvDir;

    /** Coeficiente da reacao.
     * Hostilidade do meio - taxa de crescimento.
     */    
    REAL fSigma;
    
    /** Taxa de crescimento dividido pela capacidade de suporte.
     */
    REAL fLambdaDivK;
    
    /** Termo fonte constante.
     */
    REAL fSource;   
    
};

#endif
