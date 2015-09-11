/**
 * @file TPZMatValidacaoHCurlFran2.h
 * @brief Header file for class TPZMatValidacaoHCurlFran2.\n
 */

#ifndef TPZMatValidacaoHCurlFran2_H
#define TPZMatValidacaoHCurlFran2_H

#include "TPZVecL2.h"
#include "pzaxestools.h"
#include "pzvec_extras.h"

const REAL M_C  (3*1e8); //velocidade da luz no vacuo
const REAL M_UZERO  (4*M_PI*1e-7);//permeabilidade do meio livre
const REAL M_EZERO  (8.854*1e-12);//permissividade do meio livre
#ifdef STATE_COMPLEX
const STATE imaginary(0.,1.);//unidade imaginaria
#endif
/**
 * @ingroup material
 * @brief This class implements the weak statement of the model problem from Oden's book, Chapter 1, within the PZ environment
 */
class  TPZMatValidacaoHCurlFran2 : public TPZVecL2
{
    
protected:
  
  //COM CERTEZA
  STATE (& fUr)( TPZVec<REAL>);
  STATE (& fEr)( TPZVec<REAL>);
  REAL fFreq;//frequencia da onda
  REAL fW;
  REAL fTheta;
   
	
public:
    TPZMatValidacaoHCurlFran2(int id, REAL freq, STATE (& ur)( TPZVec<REAL>),STATE (& er)( TPZVec<REAL>), REAL t);
  
    TPZMatValidacaoHCurlFran2(int id);
  
    /** @brief Default constructor */
    TPZMatValidacaoHCurlFran2();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatValidacaoHCurlFran2(const TPZMatValidacaoHCurlFran2 &mat);
    /** @brief Default destructor */
    virtual ~TPZMatValidacaoHCurlFran2();
	
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMatValidacaoHCurlFran2"; }
    
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
  
  /**
   * @brief This method defines which parameters need to be initialized in order to compute the contribution of an element
   * @param datavec [out] vector of TPZMaterialData, each position will specifie the requirements for its correspondent state variable
   */
  virtual void FillDataRequirements(TPZMaterialData &data)
  {
    data.SetAllRequirements(false);
    data.fNeedsNormal = true;
  }
  
  
  
  /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
  virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
  {
    data.SetAllRequirements(false);
    data.fNeedsNormal = true;
  }
  
  /**
   * @brief This method defines which parameters need to be initialized in order to compute the contribution of an element
   * @param datavec [out] vector of TPZMaterialData, each position will specifie the requirements for its correspondent state variable
   */
  virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
  {
    int nref = datavec.size();
    for(int iref = 0; iref<nref; iref++){
      datavec[iref].SetAllRequirements(false);
      datavec[iref].fNeedsNormal = true;
    }
  }
  
  
  
  /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
  virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
  {
    // default is no specific data requirements
    int nref = datavec.size();
    for(int iref = 0; iref<nref; iref++){
      datavec[iref].SetAllRequirements(false);
      datavec[iref].fNeedsNormal = true;
    }
  }
  
  virtual int VariableIndex(const std::string &name);
  
  /**
   * @brief Returns the number of variables associated with the variable indexed by var.
   * @param var Index variable into the solution, is obtained by calling VariableIndex
   */
  virtual int NSolutionVariables(int var);
  
  /** @brief Returns the solution associated with the var index based on the finite element approximation */
  virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
  
  void SetTheta(REAL t) {fTheta=t;}
};



STATE urDefault( TPZVec<REAL>x );//default material has ur=1
STATE erDefault( TPZVec<REAL>x );//default material has er=1
#endif

