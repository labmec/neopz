//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZHybridDarcyFlow.h"
#include "pzaxestools.h"

TPZHybridDarcyFlow::TPZHybridDarcyFlow() : TPZRegisterClassId(&TPZHybridDarcyFlow::ClassId),
                               TPZMatCombinedSpacesT<STATE>(), TPZDarcyFlow() {}

TPZHybridDarcyFlow::TPZHybridDarcyFlow(int id, int dim) : TPZRegisterClassId(&TPZHybridDarcyFlow::ClassId),
                                        TPZDarcyFlow(id,dim) {}







int TPZHybridDarcyFlow::ClassId() const {
    return Hash("TPZHybridDarcyFlow") ^ TPZDarcyFlow::ClassId() << 1;
}

TPZMaterial *TPZHybridDarcyFlow::NewMaterial() const {
    return new TPZHybridDarcyFlow(*this);
}

void TPZHybridDarcyFlow::Print(std::ostream &out) const {
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << TPZDarcyFlow::Id() << "\n";
    out << "Dimension: " << TPZDarcyFlow::Dimension() << "\n\n";
}

/** @name Contribute */
/** @{ */
/**
 * @brief It computes a contribution to the stiffness matrix
 * and load vector at one integration point.
 * @param[in] datavec stores all input data
 * @param[in] weight is the weight of the integration rule
 * @param[out] ek is the element matrix
 * @param[out] ef is the rhs vector
 */
void TPZHybridDarcyFlow::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                        REAL weight,TPZFMatrix<STATE> &ek,
                        TPZFMatrix<STATE> &ef)
{
    const TPZFMatrix<REAL> &phi = datavec[0].fPhi;
    const TPZFMatrix<REAL> &dphi = datavec[0].dphix;
    const TPZVec<REAL> &x = datavec[0].x;
    const TPZFMatrix<REAL> &axes = datavec[0].axes;
    const TPZFMatrix<REAL> &jacinv = datavec[0].jacinv;
    auto phr = dphi.Cols();

    const STATE perm = GetPermeability(datavec[0].x);

    STATE source_term = 0;
    if (this->HasForcingFunction()) {
        TPZManVector<STATE, 1> res(1);
        fForcingFunction(x, res);
        source_term = -res[0];
    }
    
    // Darcy's equation
    for (int in = 0; in < phr; in++) {
        ef(in, 0) -= weight * source_term * (phi(in, 0));
        for (int jn = 0; jn < phr; jn++) {
            for (int kd = 0; kd < fDim; kd++) {
                ek(in, jn) += weight * (dphi(kd, in) * perm * dphi(kd, jn));
            }
        }
    }
    auto nspaces = datavec.size();
    for(int ispace = 2; ispace < nspaces; ispace+= 2)
    {
        for (int in = 0; in <phr; in++) {
            ek(in,phr+(ispace-2)*2) += weight*phi(in,0);
            ek(phr+(ispace-2)*2,in) += weight*phi(in,0);
        }
        ek(phr+(ispace-2)*2,phr+(ispace-2)*2+1) -= weight;
        ek(phr+(ispace-2)*2+1,phr+(ispace-2)*2) -= weight;
    }
}
/**@}*/


/**
 * @brief It computes a contribution to the stiffness matrix
 * and load vector at one BC integration point.
 * @param[in] datavec stores all input data
 * @param[in] weight is the weight of the integration rule
 * @param[out] ek is the element matrix
 * @param[out] ef is the rhs vector
 * @param[in] bc is the boundary condition material
 */
void TPZHybridDarcyFlow::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                          REAL weight, TPZFMatrix<STATE> &ek,
                          TPZFMatrix<STATE> &ef,
                          TPZBndCondT<STATE> &bc)
{
    TPZDarcyFlow::ContributeBC(datavec[1],weight,ek,ef,bc);
}
