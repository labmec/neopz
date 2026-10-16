/**
 * @file TPZCylindricalPML.h
 * @brief Header file for class TPZCylindricalPML.\n
 */
#ifndef _TPZCYLINDRICALPML_H_
#define _TPZCYLINDRICALPML_H_

#include "TPZMatPML.h"

template<class T>
class TPZVec;
template<class T>
class TPZFMatrix;

/**
 * @ingroup material
 * @brief This class implements a cylindrical PML for a given existing material.
 * The material should access its permittivity and permeability through
 * methods called GetPermeability and GetPermittivity
 * with the same signature as shown here.
 */

template<class TMAT>
class  TPZCylindricalPML : public TPZMatPML<TMAT>
{
protected:
  bool fAttR{false};
  REAL fPmlBeginR{-1};
  bool fAttZ{false};
  REAL fPmlBeginZ{-1};
  
  CSTATE fAlphaMaxR{0};
  CSTATE fAlphaMaxZ{0};
  REAL fDR{-1};
  REAL fDZ{-1};
  TPZCylindricalPML() = default;

  void ComputeTransformMat(TPZFMatrix<CSTATE> &mat,
                           const REAL rho,
                           const REAL phi,
                           const REAL z) const;
public:
  //! Creates PML based on another domain region
  TPZCylindricalPML(const int id, const TMAT &mat) : TPZMatPML<TMAT>(id,mat) {};
  //! Sets information regarding the attenuation of the PML in the r-direction
  void SetAttR(const REAL pmlBegin, const CSTATE alpha, const REAL d);
  //! Sets information regarding the attenuation of the PML in the z-direction
  void SetAttZ(const REAL pmlBegin, const CSTATE alpha, const REAL d);
  //! Gets information regarding the attenuation of the PML in the r-direction
  void GetAttR(REAL& pmlBegin, CSTATE& alpha, REAL& d) const
  {pmlBegin=fPmlBeginR;alpha=fAlphaMaxR;d=fDR;}
  void GetAttZ(REAL& pmlBegin, CSTATE& alpha, REAL& d) const
  {pmlBegin=fPmlBeginZ;alpha=fAlphaMaxZ;d=fDZ;}
  //! Gets the permeability of the material
  void GetPermeability(const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &ur) const override;
  //! Gets the permittivity of the material
  void GetPermittivity(const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &er) const override;
    
  TPZCylindricalPML * NewMaterial() const override;
    
  std::string Name() const override { return "TPZCylindricalPML";}
  
  int ClassId() const override;
};

template<class T>
class TPZSingleSpaceCylindricalPML : public virtual TPZCylindricalPML<T>{
public:
  using TPZCylindricalPML<T>::TPZCylindricalPML;
  [[nodiscard]] int IntegrationRuleOrder(const int elPMaxOrder) const override{
    return this->MyIntegrationRuleOrder(elPMaxOrder);
  }
};

template<class T>
class TPZCombinedSpacesCylindricalPML : public virtual TPZCylindricalPML<T>{
public:
  using TPZCylindricalPML<T>::TPZCylindricalPML;
  [[nodiscard]] int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override{
    return this->MyIntegrationRuleOrder(elPMaxOrder);
  }
};


#include "TPZCylindricalPML_impl.h"
#endif /* _TPZCYLINDRICALPML_H_ */
