/**
 * @file TPZPlanarWgScatt.h
 * @brief Header file for class TPZPlanarWgScatt.\n
 */

#ifndef TPZPLANARSCATTERING_H
#define TPZPLANARSCATTERING_H


#include "TPZScalarField.h"
/**
 * @ingroup material
 * @brief This class implements the weak statement for the scattering analysis of planar waveguides using H1 elements.
 * It can either approximate TE or TM modes of a planar waveguide, periodic or not.
 * This source is implemented through the material TPZPlanarWgScattSrc.
 * @note Formulation taken from: 
 Yasuhide Tsuji and Masanori Koshiba, "Finite Element Method Using Port Truncation by Perfectly Matched Layer Boundary Conditions for Optical Waveguide Discontinuity Problems," J. Lightwave Technol. 20, 463- (2002) 
*/
class  TPZPlanarWgScatt :
  public TPZScalarField
{
public:
  //! Parent class's constructor
  using TPZScalarField::TPZScalarField;
  
  TPZPlanarWgScatt * NewMaterial() const override;
    
  std::string Name() const override { return "TPZPlanarWgScatt"; }
  
  /**
     @name ContributeMethods
     @{
  */
  void Contribute(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                  TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef) override;

  void Contribute(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                  TPZFMatrix<CSTATE> &ef) override {}//nothing to be done
  
  void ContributeBC(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                    TPZBndCondT<CSTATE> &bc) override;
  /**@}*/

  [[nodiscard]] int ClassId() const override;
protected:
  TPZPlanarWgScatt() = default;
};

#endif
