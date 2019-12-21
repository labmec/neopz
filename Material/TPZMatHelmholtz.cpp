#include "TPZMatHelmholtz.h"
#include "pzaxestools.h"
#include "TPZCompElHCurl.h"
#include "pzbndcond.h"
#include "pzlog.h"

TPZMatHelmholtz::TPZMatHelmholtz(int dim, int id, const STATE &c,const REAL &scale)
    : TPZVecL2(id), fC(c) , fScale(scale) {
    TPZVecL2::fDim = dim;
    fCurlDim = 2 * fDim - 3;
}


/** @brief Default constructor */
TPZMatHelmholtz::TPZMatHelmholtz() : TPZVecL2(), fC(1.) , fScale(1.) {
    TPZVecL2::fDim = 3;
    fCurlDim = 2 * fDim - 3;
}

void TPZMatHelmholtz::Contribute(TPZMaterialData &data, REAL weight,
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

    const int nHCurlFunctions = phiHCurl.Rows();
    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        STATE load = 0.;
        for(auto x = 0; x < fDim; x++)   load += phiHCurl(iVec, x) * force[x];
        ef(iVec) += load * weight;
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {

            STATE phiIdotPhiJ = 0.;
            for(auto x = 0; x < fDim; x++)   phiIdotPhiJ += phiHCurl(iVec, x) * phiHCurl(jVec, x);

            STATE curlPhiIvecCurlPhiJ = 0.;
            for(auto x = 0; x < fCurlDim; x++)   curlPhiIvecCurlPhiJ += curlPhi(x, iVec) * curlPhi(x, jVec);

            const STATE stiff = curlPhiIvecCurlPhiJ + fC * fScale * fScale * phiIdotPhiJ;
            ek(iVec, jVec) += stiff * weight;
        }
    }
}

void TPZMatHelmholtz::ContributeBC(TPZMaterialData &data, REAL weight,
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
    TPZFNMatrix<30,REAL> phiHCurl;
    TPZHCurlAuxClass::ComputeShape(data.fVecShapeIndex,data.phi,data.fDeformedDirections,phiHCurl);
    TPZFMatrix<REAL> &curlPhi = data.curlphi;
    const int nHCurlFunctions = phiHCurl.Rows();
    const int phiDim = phiHCurl.Cols();
    const REAL BIG = TPZMaterial::gBigNumber;

    // const STATE v1 = bc.Val1()(0,0);
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

int TPZMatHelmholtz::IntegrationRuleOrder(int elPMaxOrder) const {
    return 2 + elPMaxOrder * 2;
}

int TPZMatHelmholtz::VariableIndex(const std::string &name) {
    if (strcmp(name.c_str(), "E") == 0)
        return 0;
	if (strcmp(name.c_str(), "curlE") == 0)
		return 1;
    return TPZMaterial::VariableIndex(name);
}

/**
 * @brief Returns the number of variables associated with the variable indexed
 * by var.
 * @param var Index variable into the solution, is obtained by calling
 * VariableIndex
 */
int TPZMatHelmholtz::NSolutionVariables(int var) {
    switch (var) {
    case 0: // E
        return fDim;
        break;
    case 1: // curlE
        return 2*fDim-3;
        break;
    }
    return TPZMaterial::NSolutionVariables(var);
}

/** @brief Returns the solution associated with the var index based on the
 * finite element approximation */
void TPZMatHelmholtz::Solution(TPZMaterialData &data, int var,
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
            TPZMaterial::Solution(data,var,Solout);
            break;
    }
}

void TPZMatHelmholtz::ErrorsHdiv(TPZMaterialData &data,
                                   TPZVec<STATE> &u_exact,
                                   TPZFMatrix<STATE> &curlU_exact,
                                   TPZVec<REAL> &values) {
    values.Fill(0.0);
    TPZVec<STATE> u(3, 0.), curlU(2*fDim-3, 0.);

    Solution(data, 0, u);     // E
    Solution(data, 1, curlU); // curlE

    // values[0] : E error using HCurl norm (values[1]+values[2])
    // values[1] : E error using L2 norm
    // values[2] : curlE error using L2 norm

	// values[1] : E error using L2 norm
    for (int id = 0; id < fDim; id++) {
        STATE diffE = u[id]-u_exact[id];
		values[1] += std::norm(std::conj(diffE)*diffE);
    }
	
	// values[2] : curlE error using L2 norm
	STATE diffCurl{0};
    for(auto x = 0; x < fCurlDim; x++) diffCurl+= curlU[x] - curlU_exact(x,0);
	values[2] = std::norm(std::conj(diffCurl)*diffCurl);
	
	// values[0] : E error using HCurl norm (values[1]+values[2])
    values[0] = values[1] + values[2];
}
