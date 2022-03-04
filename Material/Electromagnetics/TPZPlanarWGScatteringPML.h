/**

 * @file TPZPlanarWGScattering.h
 * @brief Header file for class TPZPlanarWGScattering.\n
 */

#ifndef TPZPLANARSCATTERINGPML_H
#define TPZPLANARSCATTERINGPML_H


#include "TPZPlanarWGScattering.h"

class  TPZPlanarWGScatteringPML :
  public TPZPlanarWGScattering
{
protected:
  bool fAttX{false};
  REAL fPmlBeginX{-1};
  bool fAttY{false};
  REAL fPmlBeginY{-1};
    
  STATE fAlphaMaxX{-1};
  STATE fAlphaMaxY{-1};
  STATE fDX{-1};
  STATE fDY{-1};
  TPZPlanarWGScatteringPML() = default;

  void ComputeSParameters(const TPZVec<REAL> &x, CSTATE&sx, CSTATE&sy, CSTATE &sz) const;
public:
  
  //! Creates PML based on another domain region
  TPZPlanarWGScatteringPML(const int id,
                               const TPZPlanarWGScattering &mat);
  //! Sets information regarding the attenuation of the PML in the x-direction
  void SetAttX(const REAL pmlBegin, const STATE alpha, const REAL d);
  //! Sets information regarding the attenuation of the PML in the y-direction
  void SetAttY(const REAL pmlBegin, const STATE alpha, const REAL d);

  //! Gets the permeability of the material
  TPZVec<CSTATE> GetPermeability(const TPZVec<REAL> &x) const override;
  //! Gets the permittivity of the material
  TPZVec<CSTATE> GetPermittivity(const TPZVec<REAL> &x) const override;
    
  TPZPlanarWGScatteringPML * NewMaterial() const override;
    
  std::string Name() const override { return "TPZPlanarWGScatteringPML";}

  int IntegrationRuleOrder(const int elPMaxOrder) const override;
};

#endif
