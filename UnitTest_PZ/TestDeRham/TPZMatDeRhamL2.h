#ifndef TPZMATDERHAML2_H
#define TPZMATDERHAML2_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"


class TPZMatDeRhamL2 :
  public TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>>
{
  using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>>;
    
public:
    TPZMatDeRhamL2(int matId, int dim) : TBase(matId), fDim(dim) {}

    [[nodiscard]] int Dimension() const override { return fDim;}

    [[nodiscard]] int NStateVariables() const override { return 1;}

    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight,
                    TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;
    
    inline void ContributeBC(const TPZMaterialDataT<STATE> &data,
                             REAL weight,
                             TPZFMatrix<STATE> &ek,
                             TPZFMatrix<STATE> &ef,
                             TPZBndCondT<STATE> &bc) override {}
private:
    const int fDim{-1};
};
#endif