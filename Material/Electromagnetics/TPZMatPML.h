/**
 * @file TPZMatPML.h
 * @brief Header file for class TPZMatPML.\n
 */
#ifndef _TPZMATPML_H_
#define _TPZMATPML_H_

#include "pzreal.h"
/**
 * @ingroup material
 * @brief This class implements the interface for an anisotropic PML for a given existing material.
 * The material should access its permittivity and permeability through
 * methods called GetPermeability and GetPermittivity
 * with the same signature as shown here.
 */
template<class TMAT>
class  TPZMatPML : public TMAT
{
protected:
  //! For single space materials
  int MyIntegrationRuleOrder(const int elPMaxOrder) const{
    return 2 + TMAT::IntegrationRuleOrder(elPMaxOrder);
  }
  //! For combined spaces materials
  int MyIntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const{
    return 2 + TMAT::IntegrationRuleOrder(elPMaxOrder);
  }
  TPZMatPML() = default;
public:
  //! Creates PML based on another domain region
  TPZMatPML(const int id, const TMAT &mat) : TMAT(mat) {this->SetId(id);};
  //! Gets the permeability of the material
  void GetPermeability(const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &ur) const override = 0;
  //! Gets the permittivity of the material
  void GetPermittivity(const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &er) const override = 0;
    
  std::string Name() const override { return "TPZMatPML";}
};

#endif /* _TPZMATPML_H_ */
