/**
 * @file TPZMatPML.h
 * @brief Header file for class TPZMatPML.\n
 */
#ifndef _TPZMATPML_H_
#define _TPZMATPML_H_

#include "pzreal.h"
/**
 * @ingroup material
 * @brief This class implements a rectangular UPML for a given existing material.
 * The material should access its permittivity and permeability through
 * methods called GetPermeability and GetPermittivity
 * with the same signature as shown here.
 */
template<class TMAT>
class  TPZMatPML : public TMAT
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
  TPZMatPML() = default;

  void ComputeSParameters(const TPZVec<REAL> &x,
                          CSTATE&sx, CSTATE& sy, CSTATE& sz) const;
public:
  //! Creates PML based on another domain region
  TPZMatPML(const int id, const TMAT &mat);
  //! Sets information regarding the attenuation of the PML in the x-direction
  void SetAttX(const REAL pmlBegin, const STATE alpha, const REAL d);
  //! Sets information regarding the attenuation of the PML in the y-direction
  void SetAttY(const REAL pmlBegin, const STATE alpha, const REAL d);
  //! Sets information regarding the attenuation of the PML in the z-direction
  void SetAttZ(const REAL pmlBegin, const STATE alpha, const REAL d);
  //! Gets the permeability of the material
  void GetPermeability(const TPZVec<REAL> &x,TPZVec<CSTATE> &ur) const override;
  //! Gets the permittivity of the material
  void GetPermittivity(const TPZVec<REAL> &x,TPZVec<CSTATE> &er) const override;
    
  TPZMatPML * NewMaterial() const override;
    
  std::string Name() const override { return "TPZMatPML";}
  
  int ClassId() const override;
};

template<class T>
class TPZSingleSpacePML : public virtual TPZMatPML<T>{
public:
  using TPZMatPML<T>::TPZMatPML;
  [[nodiscard]] int IntegrationRuleOrder(const int elPMaxOrder) const override;
};

template<class T>
class TPZCombinedSpacesPML : public virtual TPZMatPML<T>{
public:
  using TPZMatPML<T>::TPZMatPML;
  [[nodiscard]] int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override;
};


#include "TPZMatPML_impl.h"
#endif /* _TPZMATPML_H_ */
