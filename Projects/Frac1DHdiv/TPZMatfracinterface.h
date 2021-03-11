#ifndef TPZMatfrac1dhdiv_H
#define TPZMatfrac1dhdiv_H

#include "TPZMaterial.h"
#include "tpzautopointer.h"
#include "TPZFracData.h"


/**
 * @ingroup material
 * @author Omar Duran and Nathan Shauer
 * @since 19/08/2014
 * @brief Material to couple fracture simulation with darcy
 * @brief It only implements ContributeInterface
 * @brief DOCUMENTATION OF WEAK FORMULATION IN LYX LOCATED AT THE SVN REPOSITORY
 */
class TPZMatfracinterface : public TPZMaterial {
  
protected:
  
  /** @brief Problem dimension */
  int fDim;
  
public:
  
  /** @brief Default Constructor */
  TPZMatfracinterface();
  
  /** @brief Constructor with matid */
  TPZMatfracinterface(int matid);

  /** @brief Destructor */  
  virtual ~TPZMatfracinterface();
  
  /** @brief copy constructor */
  TPZMatfracinterface(const TPZMatfracinterface &copy);

  /** @brief operator equal */
  TPZMatfracinterface &operator=(const TPZMatfracinterface &copy);
  
  /** @brief Print Method */
  virtual void Print(std::ostream & out);
  
  /** @brief Name of the material */
  virtual std::string Name() { return "TPZMatfracinterface"; }
  
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
  virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
  {
    DebugStop();
  }
  
  /** @brief ContributeBC method for not multiphysics materials */
  virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
  {
    DebugStop();
  }

  /** @brief ContributeBC method beeing used */
  virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
  {
    DebugStop();
  }
  
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
  {
    DebugStop();
  }
  
  virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
  
private:

  /** @brief Data of the simulation */
  TPZAutoPointer<TPZFracData> fData;
  
public:
  
  /** @brief Sets data of the simulation */
  void SetSimulationData(TPZAutoPointer<TPZFracData> Data) { fData = Data;}
    
};

#endif
