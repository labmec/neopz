/**
 * @file TPZMatHelmholtz.h
 * @brief Header file for class TPZMatHelmholtz.\n
 */

#ifndef TPZMATHELMHOLTZ_H
#define TPZMATHELMHOLTZ_H

#include "TPZMatHCurlProjection.h"
/**
 * @ingroup material
 * @brief This class implements the inhomogeneous Helmholtz wave equation in 2D.
 */
class TPZMatHelmholtz : public TPZMatHCurlProjection {
  protected:
    const STATE fC;
  public:
    TPZMatHelmholtz(int dim, int id, const STATE &, const REAL &scale = 1.);

    /** @brief Default constructor */
    TPZMatHelmholtz();

    /** @brief Returns the name of the material */
    std::string Name() override { return "TPZMatHelmholtz"; }

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector
     * at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    void Contribute(TPZMaterialData &data, REAL weight,
                            TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
};

#endif
