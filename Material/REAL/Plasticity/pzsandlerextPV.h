//
//  pzsandlerextPV.h
//  PZ
//
//  Created by Diogo Cecilio on 9/3/13.
//
//

#ifndef __PZ__pzsandlerextPV__
#define __PZ__pzsandlerextPV__

#include <iostream>
#include "TPZTensor.h"
#include "TPZElasticResponse.h"
#include "TPZPlasticState.h"


class TPZSandlerExtended
{
public:
    
    
    
    TPZSandlerExtended(REAL A, REAL B,REAL C, REAL D,REAL K,REAL G,REAL W,REAL R,REAL Phi,REAL N,REAL Psi);
    TPZSandlerExtended(const TPZSandlerExtended & copy);
    TPZSandlerExtended();
    ~TPZSandlerExtended();
    
    REAL F(REAL x,REAL phi);
    REAL X(REAL k);
    REAL EpsEqX(REAL X);
    REAL EpsEqk(REAL k);
    REAL ResL(TPZTensor<REAL>::TPZDecomposed pt, REAL theta,REAL beta,REAL k,REAL kprev );
    REAL DistF1(TPZTensor<REAL>::TPZDecomposed pt,REAL xi,REAL beta);
    REAL DistF2(TPZTensor<REAL>::TPZDecomposed pt,REAL theta,REAL beta,REAL k);
    void Firstk(REAL &epsp,REAL &k);
    TPZFMatrix<REAL> DDistFunc1(TPZTensor<REAL>::TPZDecomposed pt,REAL xi,REAL beta);
    TPZFMatrix<REAL> DDistFunc2(TPZTensor<REAL>::TPZDecomposed pt,REAL theta,REAL beta,REAL k,REAL kprev);
    TPZFMatrix<REAL> D2DistFunc1(TPZTensor<REAL>::TPZDecomposed pt,REAL xi,REAL beta);
    TPZFMatrix<REAL> D2DistFunc2(TPZTensor<REAL>::TPZDecomposed pt,REAL theta,REAL beta,REAL k);
    void GetRotMatrix(TPZFMatrix<REAL> &Rot);
    
    TPZManVector<REAL> FromHWCylToCart(TPZManVector<REAL> &HWCylCoords);
    //TPZManVector<REAL> FromHWCartToHWCyl(TPZManVector<REAL>&HWCartCoords);
    TPZManVector<REAL> FromPrincipalToHWCyl(TPZTensor<REAL>::TPZDecomposed &PrincipalCoords);
    TPZManVector<REAL> FromPrincipalToHWCart(TPZTensor<REAL>::TPZDecomposed &PrincipalCoords);
    
    REAL D2DisFuncF2ThetaBeta();
    
    
    TPZManVector<REAL> F1Cyl(REAL xi,REAL beta);
    TPZManVector<REAL> F2Cyl(REAL theta,REAL beta,REAL k);
    TPZManVector<REAL> FromHWCylToPrincipal(TPZManVector<REAL> &HWCylCoords);
    
    void F1(REAL xi,REAL beta,TPZManVector<REAL> &sol);
    
    void F2(REAL theta,REAL beta,REAL k,TPZManVector<REAL> &sol);
    
    void YieldFunction(TPZTensor<REAL>::TPZDecomposed &sigma, TPZVec<REAL> &yield,REAL &kprev);
    
    void ProjectF1(TPZTensor<REAL>::TPZDecomposed &sigmatrial, TPZTensor<REAL>::TPZDecomposed &sigproj);
    
    void ProjectF2(TPZTensor<REAL>::TPZDecomposed &sigmatrial, TPZTensor<REAL>::TPZDecomposed &sigproj,REAL &kprev);
    
    void ProjectRing(TPZTensor<REAL>::TPZDecomposed &sigmatrial, TPZTensor<REAL>::TPZDecomposed &sigproj,REAL &kprev);
    
    /**
	 * Imposes the specified strain tensor and returns the correspondent stress state.
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 * @param[out] sigma Resultant stress
	 */
	virtual void ApplyStrainComputeSigma(TPZPlasticState<REAL> &plasticstate, TPZTensor<REAL> &sigma);
    
    void ApplyStrainComputeSigma(TPZTensor<REAL> &eps, TPZTensor<REAL> &sigma,REAL &epspv);
    
    
public:
    
    REAL fA,fB,fC,fD,fW,fK,fR,fG,fPhi,fN,fPsi,fk0;
    bool fIsonCap;
    
    
};


#endif /* defined(__PZ__pzsandlerextPV__) */
