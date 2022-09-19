//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZHybridDarcyFlow.h"
#include "pzaxestools.h"

TPZHybridDarcyFlow::TPZHybridDarcyFlow() : TPZRegisterClassId(&TPZHybridDarcyFlow::ClassId),
                               TPZMatCombinedSpacesT<STATE>(), TPZDarcyFlow() {}

TPZHybridDarcyFlow::TPZHybridDarcyFlow(int id, int dim) : TPZRegisterClassId(&TPZHybridDarcyFlow::ClassId),
                                        TPZDarcyFlow(id,dim) {}




int TPZHybridDarcyFlow::VariableIndex(const std::string &name) const {

    if (!strcmp("Solution", name.c_str())) return 1;
    if (!strcmp("Pressure", name.c_str())) return 1;
    if (!strcmp("Derivative", name.c_str())) return 2;
    if (!strcmp("GradU", name.c_str())) return 2;
    if (!strcmp("KDuDx", name.c_str())) return 3;
    if (!strcmp("KDuDy", name.c_str())) return 4;
    if (!strcmp("KDuDz", name.c_str())) return 5;
    if (!strcmp("NormKDu", name.c_str())) return 6;
    if (!strcmp("MinusKGradU", name.c_str())) return 7;
    if (!strcmp("Flux", name.c_str())) return 7;
    if (!strcmp("POrder", name.c_str())) return 8;
    if (!strcmp("ExactPressure", name.c_str())) return 9;
    if (!strcmp("ExactSolution", name.c_str())) return 9;
    if (!strcmp("ExactFlux", name.c_str())) return 10;
    if (!strcmp("Div", name.c_str())) return 11;
    if (!strcmp("Divergence", name.c_str())) return 11;
    if (!strcmp("ExactDiv", name.c_str())) return 12;
    if (!strcmp("ExactDivergence", name.c_str())) return 12;
    if (!strcmp("FluxL2", name.c_str())) return 13;

    return TPZDarcyFlow::VariableIndex(name);
}

int TPZHybridDarcyFlow::NSolutionVariables(int var) const {

    if (var == 1) return 1;      // Solution/Pressure
    if (var == 2) return fDim;   // Derivative/GradU
    if (var == 3) return 1;      // KDuDx;
    if (var == 4) return 1;      // KDuDy;
    if (var == 5) return 1;      // KDuDz;
    if (var == 6) return 1;      // NormKDu;
    if (var == 7) return fDim;   // MinusKGradU/Flux;
    if (var == 8) return 1;      // POrder
    if (var == 9) return 1;      // ExactPressure/ExactSolution
    if (var == 10) return fDim;  // ExactFlux
    if (var == 11) return 1;     // Div/Divergence
    if (var == 12) return 1;     // ExactDiv/ExactDivergence
    if (var == 13) return fDim;  // FluxL2

    return TPZDarcyFlow::NSolutionVariables(var);
}



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
    for(int ispace = 1; ispace < nspaces; ispace+= 2)
    {
        for (int in = 0; in <phr; in++) {
            ek(in,phr+(ispace-1)*2) += weight*phi(in,0);
            ek(phr+(ispace-1)*2,in) += weight*phi(in,0);
        }
        ek(phr+(ispace-1)*2,phr+(ispace-1)*2+1) -= weight;
        ek(phr+(ispace-1)*2+1,phr+(ispace-1)*2) -= weight;
    }
}
/**@}*/
