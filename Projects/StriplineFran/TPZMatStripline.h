/**
 * @file TPZMatStripline.h
 * @brief Header file for class TPZMatStripline.\n
 * It implements the weak statement of the model problem from Oden's book, Chapter 1 within the PZ environment.
 */

#ifndef TPZMatStripline_H
#define TPZMatStripline_H

#include "TPZVecL2.h"
#include "pzaxestools.h"

const REAL M_C  (3*1e8); //velocidade da luz no vacuo
const REAL M_UZERO  (4*M_PI*1e-7);//permeabilidade do meio livre
const REAL M_EZERO  (8.854*1e-12);//permissividade do meio livre
const STATE imaginary(0.,1.);//unidade imaginaria
/**
 * @ingroup material
 * @brief This class implements the weak statement of the model problem from Oden's book, Chapter 1, within the PZ environment
 */
class  TPZMatStripline : public TPZVecL2
{
    
protected:
  
  //COM CERTEZA
  STATE (& fUr)( TPZVec<REAL>);
  STATE (& fEr)( TPZVec<REAL>);
  REAL fFreq;//frequencia da onda
  STATE fW;
   
	
public:
    TPZMatStripline(int id, REAL freq, STATE (& ur)( TPZVec<REAL>),STATE (& er)( TPZVec<REAL>));
  
    TPZMatStripline(int id);
  
    /** @brief Default constructor */
    TPZMatStripline();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatStripline(const TPZMatStripline &mat);
    /** @brief Default destructor */
    virtual ~TPZMatStripline();
	
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMatStripline"; }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const {return 2;}
  
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() { return 1;}
    
public:
    
  
  /**
   * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
   * @param data [in] stores all input data
   * @param weight [in] is the weight of the integration rule
   * @param ek [out] is the stiffness matrix
   * @param ef [out] is the load vector
   * @since April 16, 2007
   */
  virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
  
  
  /**
   * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
   * @param datavec [in] stores all input data
   * @param weight [in] is the weight of the integration rule
   * @param ek [out] is the stiffness matrix
   * @param ef [out] is the load vector
   */
  virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
  
  /**
   * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
   * @param datavec [in] stores all input data
   * @param weight [in] is the weight of the integration rule
   * @param ek [out] is the stiffness matrix
   * @param ef [out] is the load vector
   */
  virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
  
  /**
   * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
   * @param data [in] stores all input data
   * @param weight [in] is the weight of the integration rule
   * @param ek [out] is the stiffness matrix
   * @param ef [out] is the load vector
   * @param bc [in] is the boundary condition material
   * @since October 07, 2011
   */
  virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
  
  /**
   * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
   * to multiphysics simulation.
   * @param datavec [in]  stores all input data
   * @param weight [in] is the weight of the integration rule
   * @param ek [out] is the stiffness matrix
   * @param ef [out] is the load vector
   * @param bc [in] is the boundary condition material
   * @since October 18, 2011
   */
  virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
  
  /**
   * @brief It computes a contribution to the residual vector at one integration point.
   * @param data [in] stores all input data
   * @param weight [in] is the weight of the integration rule
   * @param ef [out] is the residual vector
   * @since April 16, 2007
   */
  virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
  
  /**
   * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
   * @param data [in] stores all input data
   * @param weight [in] is the weight of the integration rule
   * @param ef [out] is the load vector
   * @param bc [in] is the boundary condition material
   * @since April 16, 2007
   */
  virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
  
  virtual int VariableIndex(const std::string &name);
  
  /**
   * @brief Returns the number of variables associated with the variable indexed by var.
   * @param var Index variable into the solution, is obtained by calling VariableIndex
   */
  virtual int NSolutionVariables(int var);
  
  /** @brief Returns the solution associated with the var index based on the finite element approximation */
  virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
};

STATE urDefault( TPZVec<REAL>x );//default material has ur=1
STATE erDefault( TPZVec<REAL>x );//default material has er=1
#endif

