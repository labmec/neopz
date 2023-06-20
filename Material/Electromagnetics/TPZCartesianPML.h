/**
 * @file TPZCartesianPML.h
 * @brief Header file for class TPZCartesianPML.\n
 */
#ifndef _TPZCARTESIANPML_H_
#define _TPZCARTESIANPML_H_

#include "TPZMatPML.h"
/**
 * @ingroup material
 * @brief This class implements a cartesian PML for a given existing material.
 * The material should access its permittivity and permeability through
 * methods called GetPermeability and GetPermittivity
 * with the same signature as shown here.
 */

template<class TMAT>
class  TPZCartesianPML : public TPZMatPML<TMAT>
{
protected:
  bool fAttX{false};
  REAL fPmlBeginX{-1};
  bool fAttY{false};
  REAL fPmlBeginY{-1};
  bool fAttZ{false};
  REAL fPmlBeginZ{-1};
  
  STATE fAlphaMaxX{-1};
  STATE fAlphaMaxY{-1};
  STATE fAlphaMaxZ{-1};
  STATE fDX{-1};
  STATE fDY{-1};
  STATE fDZ{-1};
  TPZCartesianPML() = default;

  void ComputeSParameters(const TPZVec<REAL> &x,
                          CSTATE&sx, CSTATE& sy, CSTATE& sz) const;
public:
  //! Creates PML based on another domain region
  TPZCartesianPML(const int id, const TMAT &mat) : TPZMatPML<TMAT>(id,mat) {};
  //! Sets information regarding the attenuation of the PML in the x-direction
  void SetAttX(const REAL pmlBegin, const STATE alpha, const REAL d);
  //! Sets information regarding the attenuation of the PML in the y-direction
  void SetAttY(const REAL pmlBegin, const STATE alpha, const REAL d);
  //! Sets information regarding the attenuation of the PML in the z-direction
  void SetAttZ(const REAL pmlBegin, const STATE alpha, const REAL d);
  //! Gets the permeability of the material
  void GetPermeability(const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &ur) const override;
  //! Gets the permittivity of the material
  void GetPermittivity(const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &er) const override;
    
  TPZCartesianPML * NewMaterial() const override;
    
  std::string Name() const override { return "TPZCartesianPML";}
  
  int ClassId() const override;
};

template<class T>
class TPZSingleSpaceCartesianPML : public virtual TPZCartesianPML<T>{
public:
  using TPZCartesianPML<T>::TPZCartesianPML;
  [[nodiscard]] int IntegrationRuleOrder(const int elPMaxOrder) const override{
    return this->MyIntegrationRuleOrder(elPMaxOrder);
  }
};

template<class T>
class TPZCombinedSpacesCartesianPML : public virtual TPZCartesianPML<T>{
public:
  using TPZCartesianPML<T>::TPZCartesianPML;
  [[nodiscard]] int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override{
    return this->MyIntegrationRuleOrder(elPMaxOrder);
  }
};


#include "TPZCartesianPML_impl.h"
#endif /* _TPZCARTESIANPML_H_ */
