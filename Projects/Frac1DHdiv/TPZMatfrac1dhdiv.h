#ifndef TPZMatfrac1dhdiv_H
#define TPZMatfrac1dhdiv_H

#include "pzmaterial.h"
#include "tpzautopointer.h"
#include "TPZFracData.h"
/**
 * @ingroup material
 * @author Omar Duran and Nathan Shauer
 * @since 19/08/2014
 * @brief Material to solve a 1d mixed formulation for fracture opening
 * @brief Here is used H1 for flux and L2 for pressure
 */
class TPZMatfrac1dhdiv : public TPZMaterial {
  
protected:
  
  /** @brief Problem dimension */
  int fDim;
  
public:
  
  /** @brief Default Constructor */
  TPZMatfrac1dhdiv();
  
  /** @brief Constructor with matid */
  TPZMatfrac1dhdiv(int matid);

  /** @brief Destructor */  
  virtual ~TPZMatfrac1dhdiv();
  
  /** @brief copy constructor */
  TPZMatfrac1dhdiv(const TPZMatfrac1dhdiv &copy);

  /** @brief operator equal */
  TPZMatfrac1dhdiv &operator=(const TPZMatfrac1dhdiv &copy);
  
  /** @brief Print Method */
  virtual void Print(std::ostream & out);
  
  /** @brief Name of the material */
  virtual std::string Name() { return "TPZMatfrac1dhdiv"; }
  
  /** @brief Returns the integrable dimension */
  virtual int Dimension() const;
  
  /** @brief Return the number of state variables */
  virtual int NStateVariables();
  
  /** @brief Contribute method for not multiphysics materials */
  virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
  {
    DebugStop();
  }
  
  /** @brief Contribute method beeing used */
  virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
  
  /** @brief ContributeBC method for not multiphysics materials */
  virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
  {
    DebugStop();
  }

  /** @brief ContributeBC method beeing used */
  virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
  
  /** @brief Fill material data parameter with necessary requirements for the Contribute method*/
  virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);

  /** @brief Fill material data parameter with necessary requirements for the ContributeBC method*/
  virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
  
  /** @brief Returns the variable index associated with the name */
  virtual int VariableIndex(const std::string &name);
  
  /** @brief Returns the number of variables associated with the variable indexed by var */
  virtual int NSolutionVariables(int var);
  
  /** @brief Calculates a solution given datavec*/
  virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
  
private:

  /** @brief Data of the simulation */
  TPZAutoPointer<TPZFracData> fData;
  
public:
  
  /** @brief Sets data of the simulation */
  void SetSimulationData(TPZAutoPointer<TPZFracData> Data) { fData = Data;}
  
  /** @brief Return w based on p of the material data */
  REAL Getw(REAL p);
  
  /** @brief Return dwdp based on p of the material data */
  REAL Getdwdp();
  
};

#endif