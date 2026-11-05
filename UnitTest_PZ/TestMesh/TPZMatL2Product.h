#ifndef TPZMATL2PRODUCT_H
#define TPZMATL2PRODUCT_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"

class TPZMatL2Product :
  public TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>>
{
  using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>>;
    
public:
    TPZMatL2Product(int matId, int dim) : TBase(matId), fDim(dim) {}

    [[nodiscard]] int Dimension() const override { return fDim;}

    [[nodiscard]] int NStateVariables() const override { return 1;}

    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight,
                    TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;
    
    void ContributeVecShape(const TPZMaterialDataT<STATE> &data,
                            REAL weight, TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef);

    void ContributeScalarShape(const TPZMaterialDataT<STATE> &data,
                              REAL weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef);
    
    inline void ContributeBC(const TPZMaterialDataT<STATE> &data,
                             REAL weight,
                             TPZFMatrix<STATE> &ek,
                             TPZFMatrix<STATE> &ef,
                             TPZBndCondT<STATE> &bc) override {}
private:
    const int fDim{-1};
};
#endif