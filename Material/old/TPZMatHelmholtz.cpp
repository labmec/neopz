#include "TPZMatHelmholtz.h"
#include "pzaxestools.h"
#include "TPZCompElHCurl.h"
#include "pzbndcond.h"
#include "pzlog.h"

TPZMatHelmholtz::TPZMatHelmholtz(int dim, int id, const STATE &c,const REAL &scale)
    : TPZMatHCurlProjection(dim,id,scale), fC(c)  {
}


/** @brief Default constructor */
TPZMatHelmholtz::TPZMatHelmholtz() : TPZMatHCurlProjection(), fC(1.) {
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

    const auto nHCurlFunctions = phiHCurl.Rows();
    for (auto iVec = 0; iVec < nHCurlFunctions; iVec++) {
        STATE load = 0.;
        for(auto x = 0; x < fDim; x++)   load += phiHCurl(iVec, x) * force[x];
        ef(iVec,0) += load * weight;
        for (auto jVec = 0; jVec < nHCurlFunctions; jVec++) {

            STATE phiIdotPhiJ = 0.;
            for(auto x = 0; x < fDim; x++)   phiIdotPhiJ += phiHCurl(iVec, x) * phiHCurl(jVec, x);

            STATE curlPhiIvecCurlPhiJ = 0.;
            for(auto x = 0; x < fCurlDim; x++)   curlPhiIvecCurlPhiJ += curlPhi(x, iVec) * curlPhi(x, jVec);

            const STATE stiff = curlPhiIvecCurlPhiJ + fC * fScale * fScale * phiIdotPhiJ;
            ek(iVec, jVec) += stiff * weight;
        }
    }
}