/* 
 * File:   TPZYCCamClayPV.h
 * Author: quinelato
 *
 * Created on August 14, 2017, 5:18 PM
 */

#ifndef TPZYCCAMCLAYPV_H
#define TPZYCCAMCLAYPV_H

#include "pzreal.h"
#include "pzvec.h"
#include "pzfilebuffer.h"
#include "TPZElasticResponse.h"
#include "pzfmatrix.h"

class TPZYCCamClayPV {
public:

    enum {
        NYield = 2
    };

    TPZYCCamClayPV();
    TPZYCCamClayPV(const TPZYCCamClayPV& other);
    void SetUp(const TPZElasticResponse &ER, REAL a, REAL gamma, REAL m, REAL pt);
    void SetElasticResponse(const TPZElasticResponse &ER);
    void Read(TPZStream &buf);
    void Write(TPZStream &buf) const;
    REAL bFromP(REAL p) const;
    REAL bFromTheta(REAL theta) const;
    void Phi(TPZVec<REAL> sigvec, REAL alpha, TPZVec<REAL> &phi) const;
    
    /**
     * Computes the cylindrical coordinates of the point on the yield surface with the given xi and beta.
     * @param xi Hydrostatic component
     * @param beta Lode angle
     * @param returnValue Cylindrical coordinates of the point on the yield surface
     */
    void SurfaceInCyl(const REAL theta, const REAL beta, TPZVec<REAL> &returnValue) const;
    REAL DistanceToSurface(const TPZVec<REAL> &pt, const REAL theta, const REAL beta) const;
    
    /**
     * Computes the derivative of the distance function to the yield surface as a function of theta and beta
     * @param pt
     * @param xi
     * @param beta
     * @param fxn
     */
    void DDistanceToSurface(const TPZVec<STATE> &pt, const STATE theta, const STATE beta, TPZFMatrix<STATE> &fxn) const;
    
    /**
     * Computes the second derivative of the distance as a function of theta and beta
     * @param pt
     * @param xi
     * @param beta
     * @param jac
     */
    void D2DistanceToSurface(const TPZVec<STATE> &pt, const STATE theta, const STATE beta, TPZFMatrix<STATE> &jac) const;
    
    void ProjectToSurface(const TPZVec<REAL> &sigma_trial, const REAL kprev, TPZVec<REAL> &sigma, REAL &kproj, const REAL tol) const;
    void ProjectSigma(const TPZVec<REAL> &sigma_trial, const REAL kprev, TPZVec<REAL> &sigma, REAL &kproj) const;
    void ProjectSigmaDep(const TPZVec<REAL> &sigma_trial, const REAL kprev, TPZVec<REAL> &sigma, REAL &kproj, TPZFMatrix<REAL> &GradSigma) const;
    virtual ~TPZYCCamClayPV();
private:
    TPZElasticResponse fER;

    REAL fA;
    REAL fGamma;
    REAL fM;
    REAL fPt;
    
};

#endif /* TPZYCCAMCLAYPV_H */

