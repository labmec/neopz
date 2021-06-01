#ifndef TPZMATDERHAMH1_H
#define TPZMATDERHAMH1_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"


class TPZMatDeRhamH1 :
  public TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>>
{
  using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>>;
    
public:
    TPZMatDeRhamH1(int matId, int dim) : TBase(matId), fDim(dim) {}

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