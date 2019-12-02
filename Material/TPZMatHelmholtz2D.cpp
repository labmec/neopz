#include "TPZMatHelmholtz2D.h"
#include "pzaxestools.h"

#include "pzbndcond.h"
#include "pzlog.h"

TPZMatHelmholtz2D::TPZMatHelmholtz2D(int id,
                                     const STATE &c,
                                    const REAL &scale)
    : TPZVecL2(id), fC(c) , fScale(scale) {
    TPZVecL2::fDim = 2;
}

TPZMatHelmholtz2D::TPZMatHelmholtz2D(int id) : TPZVecL2(id), fC(1.) , fScale(1.) {
    TPZVecL2::fDim = 2;
}

/** @brief Default constructor */
TPZMatHelmholtz2D::TPZMatHelmholtz2D() : TPZVecL2(), fC(1.) , fScale(1.) {
    TPZVecL2::fDim = 2;
}

TPZMatHelmholtz2D::TPZMatHelmholtz2D(const TPZMatHelmholtz2D &mat)
    : TPZVecL2(mat), fC(mat.fC) , fScale(mat.fScale) {
    fDim = mat.fDim;
}

TPZMatHelmholtz2D::~TPZMatHelmholtz2D() {}

void TPZMatHelmholtz2D::Contribute(TPZMaterialData &data, REAL weight,
                                   TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef) {
    DebugStop();
//
//    /*********************CREATE HCURL FUNCTIONS*****************************/
//    TPZFNMatrix<36, REAL> phiHCurlAxes = data.phi;
//    TPZFNMatrix<40, REAL> curlPhiDAxes = data.dphix;
//
//    TPZFNMatrix<40, REAL> curlPhi, phiHCurl;
//
//    TPZAxesTools<REAL>::Axes2XYZ(phiHCurlAxes, phiHCurl, data.axes, false);
//
//    TPZManVector<REAL, 3> ax1(3), ax2(3), elNormal(3);
//    for (int i = 0; i < 3; i++) {
//        ax1[i] = data.axes(0, i);
//        ax2[i] = data.axes(1, i);
//    }
//    Cross(ax1, ax2, elNormal);
//    TPZFNMatrix<3, REAL> normalVec(1, 3);
//    normalVec(0, 0) = elNormal[0];
//    normalVec(0, 1) = elNormal[1];
//    normalVec(0, 2) = elNormal[2];
//    TPZAxesTools<REAL>::Axes2XYZ(curlPhiDAxes, curlPhi, normalVec);
//
//    TPZManVector<REAL, 3> x = data.x;
//
//    //*****************GET FORCING FUNCTION****************//
//    TPZManVector<STATE, 3> force(3);
//    force.Fill(0.);
//    if (fForcingFunction) {
//        fForcingFunction->Execute(data.x, force);
//    } else {
//        DebugStop(); // RHS not set!
//    }
//
//    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//
//
//    const int nHCurlFunctions = phiHCurl.Rows();
//
//    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
//        STATE load = 0.;
//        load += phiHCurl(iVec, 0) * force[0];// * (fScale*fScale);
//        load += phiHCurl(iVec, 1) * force[1];// * (fScale*fScale);
//        load += phiHCurl(iVec, 2) * force[2];// * (fScale*fScale);
//        ef(iVec) += load * weight;
//        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
//            STATE stiff = 0.;
//
//            STATE phiIdotPhiJ = 0.;
//            phiIdotPhiJ += phiHCurl(iVec, 0) * phiHCurl(jVec, 0);
//            phiIdotPhiJ += phiHCurl(iVec, 1) * phiHCurl(jVec, 1);
//            phiIdotPhiJ += phiHCurl(iVec, 2) * phiHCurl(jVec, 2);
//
//            STATE curlPhiIvecCurlPhiJ = 0.;
//            curlPhiIvecCurlPhiJ += curlPhi(0, iVec) * curlPhi(0, jVec);
//            curlPhiIvecCurlPhiJ += curlPhi(1, iVec) * curlPhi(1, jVec);
//            curlPhiIvecCurlPhiJ += curlPhi(2, iVec) * curlPhi(2, jVec);
//
//            stiff = curlPhiIvecCurlPhiJ + fC * fScale * fScale * phiIdotPhiJ;
//            ek(iVec, jVec) += stiff * weight;
//        }
//    }
}

void TPZMatHelmholtz2D::ContributeBC(TPZMaterialData &data, REAL weight,
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
	TPZFMatrix<REAL> &phiHCurl = data.phi;

    int nHCurlFunctions = phiHCurl.Rows();
    REAL BIG = TPZMaterial::gBigNumber;

    // const STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de
    // condicao mista
    const STATE v2 = bc.Val2()(0, 0); // sera posto no vetor F

    switch (bc.Type()) {
    case 0:
        for (int i = 0; i < nHCurlFunctions; i++) {
            const STATE rhs = phiHCurl(i, 0) * BIG * v2;
            ef(i, 0) += rhs * weight;
            for (int j = 0; j < nHCurlFunctions; j++) {
                const STATE stiff = phiHCurl(i, 0) * phiHCurl(j, 0) * BIG;
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

int TPZMatHelmholtz2D::IntegrationRuleOrder(int elPMaxOrder) const {
    return 2 + elPMaxOrder * 2;
}

int TPZMatHelmholtz2D::VariableIndex(const std::string &name) {
    if (strcmp(name.c_str(), "E") == 0)
        return 0;
	if (strcmp(name.c_str(), "curlE") == 0)
		return 1;
    DebugStop();
    return 1;
}

/**
 * @brief Returns the number of variables associated with the variable indexed
 * by var.
 * @param var Index variable into the solution, is obtained by calling
 * VariableIndex
 */
int TPZMatHelmholtz2D::NSolutionVariables(int var) {
    switch (var) {
    case 0: // E
        return 2;
        break;
    case 1: // curlE
        return 2;
        break;
    default:
        DebugStop();
        break;
    }
    return 1;
}

/** @brief Returns the solution associated with the var index based on the
 * finite element approximation */
void TPZMatHelmholtz2D::Solution(TPZMaterialData &data, int var,
                                 TPZVec<STATE> &Solout) {
    switch (var) {
    case 0: // E
    {
        Solout = data.sol[0];
    } break;
	case 1: // curlE
	{
		Solout[0] = data.dsol[0](2,0);
	} break;
    default:
        DebugStop();
        break;
    }
}

void TPZMatHelmholtz2D::ErrorsHdiv(TPZMaterialData &data,
                                   TPZVec<STATE> &u_exact,
                                   TPZFMatrix<STATE> &curlU_exact,
                                   TPZVec<REAL> &values) {
    values.Fill(0.0);
    TPZVec<STATE> u(3, 0.), curlU(1, 0.);

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
	STATE diffCurl = curlU[0] - curlU_exact(0,0);
	values[2] = std::norm(std::conj(diffCurl)*diffCurl);
	
	// values[0] : E error using HCurl norm (values[1]+values[2])
    values[0] = values[1] + values[2];
}
