#include "TPZMatHCurlProjection.h"
#include "pzaxestools.h"
#include "TPZCompElHCurl.h"
#include "pzbndcond.h"
#include "pzlog.h"

TPZMatHCurlProjection::TPZMatHCurlProjection(int dim, int id, const REAL &scale)
        : TPZVecL2(id), fScale(scale) {
    TPZVecL2::fDim = dim;
    fCurlDim = 2 * fDim - 3;
}


/** @brief Default constructor */
TPZMatHCurlProjection::TPZMatHCurlProjection() : TPZVecL2(), fScale(1.) {
    TPZVecL2::fDim = 3;
    fCurlDim = 2 * fDim - 3;
}

void TPZMatHCurlProjection::Contribute(TPZMaterialData &data, REAL weight,
                                 TPZFMatrix<STATE> &ek,
                                 TPZFMatrix<STATE> &ef) {
    TPZFNMatrix<30,REAL> phiHCurl;
    TPZHCurlAuxClass::ComputeShape(data.fVecShapeIndex,data.phi,data.fDeformedDirections,phiHCurl);
    TPZFMatrix<REAL> &curlPhi = data.curlphi;

    //*****************GET FORCING FUNCTION****************//
    TPZManVector<STATE, 3> force(3);
    force.Fill(0.);
    if (fForcingFunction) {
        fForcingFunction->Execute(data.x, force);
    } else {
        DebugStop(); // RHS not set!
    }

    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//

    const auto nHCurlFunctions = phiHCurl.Rows();
    for (auto iVec = 0; iVec < nHCurlFunctions; iVec++) {
        STATE load = 0.;
        for(auto x = 0; x < fDim; x++)   load += phiHCurl(iVec, x) * force[x];
        ef(iVec,0) += load * weight;
        for (auto jVec = 0; jVec < nHCurlFunctions; jVec++) {

            STATE phiIdotPhiJ = 0.;
            for(auto x = 0; x < fDim; x++)   phiIdotPhiJ += phiHCurl(iVec, x) * phiHCurl(jVec, x);

            ek(iVec, jVec) += phiIdotPhiJ * weight;
        }
    }
}

void TPZMatHCurlProjection::ContributeBC(TPZMaterialData &data, REAL weight,
                                   TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
    TPZFNMatrix<30,REAL> phiHCurl;
    TPZHCurlAuxClass::ComputeShape(data.fVecShapeIndex,data.phi,data.fDeformedDirections,phiHCurl);
    const auto nHCurlFunctions = phiHCurl.Rows();
    const auto phiDim = phiHCurl.Cols();
    const auto &BIG = TPZMaterial::gBigNumber;

    const STATE v2 = bc.Val2()(0, 0);
    switch (bc.Type()) {
        case 0:
            for (int i = 0; i < nHCurlFunctions; i++) {
                STATE rhs = 0, stiff = 0;
                for(auto x = 0; x <phiDim; x++) rhs+= phiHCurl(i, x) * BIG * v2;
                ef(i, 0) += rhs * weight;
                for (int j = 0; j < nHCurlFunctions; j++) {
                    for(auto x = 0; x <phiDim; x++) stiff+= phiHCurl(i, x) * phiHCurl(j, x) * BIG;
                    ek(i, j) += stiff * weight;
                }
            }
            break;
        case 1:
            DebugStop();
            break;
        case 2:
            DebugStop();
            break;
    }
}

int TPZMatHCurlProjection::IntegrationRuleOrder(int elPMaxOrder) const {
    return 2 + elPMaxOrder * 2;
}

int TPZMatHCurlProjection::VariableIndex(const std::string &name) {
    if (strcmp(name.c_str(), "E") == 0)
        return 0;
    if (strcmp(name.c_str(), "curlE") == 0)
        return 1;
    return TPZVecL2::VariableIndex(name);
}

/**
 * @brief Returns the number of variables associated with the variable indexed
 * by var.
 * @param var Index variable into the solution, is obtained by calling
 * VariableIndex
 */
int TPZMatHCurlProjection::NSolutionVariables(int var) {
    switch (var) {
        case 0: // E
            return fDim;
        case 1: // curlE
            return 2*fDim-3;
        default:
            return TPZVecL2::NSolutionVariables(var);
    }
}

/** @brief Returns the solution associated with the var index based on the
 * finite element approximation */
void TPZMatHCurlProjection::Solution(TPZMaterialData &data, int var,
                               TPZVec<STATE> &Solout) {
    switch (var) {
        case 0: // E
        {
            for(auto x = 0; x < fDim; x++) Solout[x] = data.sol[0][x];
        } break;
        case 1: // curlE
        {
            for(auto x = 0; x < fCurlDim; x++) Solout[x] = data.curlsol[0][x];
        } break;
        default:
            TPZVecL2::Solution(data,var,Solout);
            break;
    }
}

void TPZMatHCurlProjection::Errors(TPZVec<REAL> &x, TPZVec<STATE> &u,
                TPZFMatrix<STATE> &curlU, TPZFMatrix<REAL> &axes,
                TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &curlU_exact,
                TPZVec<REAL> &values) {
    values.Fill(0.0);

    // values[0] : E error using HCurl norm (values[1]+values[2])
    // values[1] : E error using L2 norm
    // values[2] : curlE error using L2 norm

    // values[1] : E error using L2 norm
    for (int id = 0; id < fDim; id++) {
        const STATE diffE = u[id]-u_exact[id];
        values[1] += std::norm(std::conj(diffE)*diffE);
    }

    // values[2] : curlE error using L2 norm
    for (int x = 0; x < fCurlDim; x++) {
        const STATE diffE = curlU[x] - curlU_exact(x,0);
        values[2] += std::norm(std::conj(diffE)*diffE);
    }

    // values[0] : E error using HCurl norm (values[1]+values[2])
    values[0] = values[1] + values[2];
}
