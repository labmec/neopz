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
    
    enum{NYield = 2};
    
    /// Constructor, with all parameters which define the Sandler DiMaggio model
    TPZSandlerExtended(STATE A, STATE B,STATE C, STATE D,STATE K,STATE G,STATE W,STATE R,STATE Phi,STATE N,STATE Psi);
	/// Copy constructor
    TPZSandlerExtended(const TPZSandlerExtended & copy);
	/// Empty constructor
    TPZSandlerExtended();
    
    TPZSandlerExtended & operator=(const TPZSandlerExtended & source)
	{
        fA=source.fA;
        fB=source.fB;
        fC=source.fC;
        fD=source.fD;
        fK=source.fK;
        fG=source.fG;
        fW=source.fW;
        fR=source.fR;
        fPhi=source.fPhi;
        fN=source.fN;
        fPsi=source.fPsi;
        fE=source.fE;
        fnu=source.fnu;
        fElasticResponse =source.fElasticResponse;
    
		return *this;
	}
    
	/// Desctructor
    virtual ~TPZSandlerExtended();

    STATE GetX(STATE k);
    
    /// Compute k as a function of epsp using Newton iterations
    void Firstk(STATE &epsp,STATE &k) const;
    
    TPZElasticResponse GetElasticResponse();
    
    STATE GetR();
    
    void YieldFunction(const TPZVec<STATE>  &sigma, STATE kprev, TPZVec<STATE> &yield) const;

    
    template<class T>
    T F(T x) const;
    
    STATE GetF(STATE x) const;
    
    void Read(TPZStream &buf);
    
    void Write(TPZStream &buf) const;
    

    
private:
    /// The function which defines the plastic surface
  
    
    /// Auxiliary function for Associating the position of the cap with the damage variable
    template<class T>
    T X(T k) const;

    /// compute the damage variable as a function of the X function
    template<class T>
    T EpsEqX(T X) const;
    /// Compute the damage variable as a function of the position of the cap k
    template<class T>
    T EpsEqk(T k) const;
    /// Compute the residual of the equation which defines the update of the damage variable
    template<class T>
    T ResLF2(const TPZVec<T> &pt, T theta,T beta,T k,STATE kprev ) const;    
    /// Compute the residual of the equation which defines the update of the damage variable
    template<class T>
    T ResLF2IJ(const TPZVec<T> &sigtrIJ, T theta,T k,STATE kprev ) const;

    /// Compute the residual of the equation which defines the update of the damage variable
    STATE ResLF1(const TPZVec<STATE> &sigtrial, TPZVec<STATE> &sigproj,STATE k,STATE kprev ) const;
    /// Compute the distance of sigtrial to the point on the yield surface
    STATE DistF1(const TPZVec<STATE> &pt,STATE xi,STATE beta) const;
    /// Compute the distance of sigtrial to the point on the cap
    STATE DistF2(const TPZVec<STATE> &pt,STATE theta,STATE beta,STATE k) const;
    
    /// Compute the distance considering the sigtrial is given as a funcion of I1, sqJ2
    STATE DistF2IJ(const TPZVec<STATE> &sigtrialIJ,STATE theta,STATE k) const;

    
    /// Compute the derivative of the distance function to the yield surface as a function of xi and beta
    void DDistFunc1(const TPZVec<STATE> &pt,STATE xi,STATE beta, TPZFMatrix<STATE> &ddistf1) const;
    
    /// Compute the derivative of the distance function to the yield surface as a function of xi and beta
    
    
    
    /// Compute the derivative of the distance function to the cap function and the result of ResL
    template<class T>
//    void DDistFunc2(const TPZVec<STATE> &pt,STATE theta,STATE beta,STATE k,STATE kprev, TPZVec<STATE> &ddistf2) const;
    void DDistFunc2(const TPZVec<T> &pt,T theta,T beta,T k,T kprev, TPZVec<T> &ddistf2) const;
    /// Compute the value of the equation which determines the orthogonality of the projection
    template<class T>
    void DDistF2IJ(TPZVec<T> &sigtrialIJ, T theta, T L, STATE Lprev, TPZVec<T> &ddistf2) const;

    
    /// Compute the second derivative of the distance as a function of xi and beta
    void D2DistFunc1(const TPZVec<STATE> &pt,STATE xi,STATE beta, TPZFMatrix<STATE> &d2distf1) const;
    
    /// Compute the second derivative of the distance as a function of theta, beta and k
    
    void D2DistFunc2(const TPZVec<STATE> &pt,STATE theta,STATE beta,STATE k, TPZFMatrix<STATE> &d2distf2)const;
    
    
    /// Compute the derivative of the equation which determines the evolution of k
    // the derivative are given in terms of theta, beta and k
    void DResLF2(const TPZVec<STATE> &pt, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &dresl) const;
    /// Compute the derivative of the equation which determines the evolution of k
    // the derivative are given in terms of k
    STATE DResLF1(const TPZVec<STATE> &sigtrial, const TPZVec<STATE> &sigproj, STATE k, STATE kprev) const;
    /// Compute the rotation matrix
    static void GetRotMatrix(TPZFMatrix<STATE> &Rot);
    /// Transform from Haigh Westergaard cylindrical coordinates to Haigh Westergaard cartesian coordinates
    static void FromHWCylToHWCart(const TPZVec<STATE> &HWCylCoords, TPZVec<STATE> &Cart);
    //TPZManVector<STATE> FromHWCartToHWCyl(TPZManVector<STATE>&HWCartCoords);
    /// Transform from eigenvalues to HW Cylindrical coordinates
    static void FromPrincipalToHWCyl(const TPZVec<STATE> &PrincipalCoords, TPZVec<STATE> &HWCyl);
    /// Transform from eigenvalues to HW cartesian coordinates
    static void FromPrincipalToHWCart(const TPZVec<STATE> &PrincipalCoords, TPZVec<STATE> &HWCart);

    /// Transform from HW Cylindrical coordinates to eigenvalues
    static void FromHWCylToPrincipal(const TPZVec<STATE> &HWCylCoords, TPZVec<STATE> &PrincipalCoords);
    
    template<class T>
    void FromThetaKToSigIJ(const T &theta, const T &K, TPZVec<T> &sigIJ) const;
    
    /// Compute the derivative of the stress (principal s;tresses) as a function of xi and beta
    void DF1Cart(STATE xi, STATE beta, TPZFMatrix<STATE> &DF1) const;
    
    /// Compute the derivative of the stress (principal s;tresses) as a function of xi and beta
    void DF2Cart(STATE theta, STATE beta, STATE k, TPZFMatrix<STATE> &DF1) const;

    
    /// Compute the derivative of the residual with respect to sigtrial
    void GradF1SigmaTrial(const TPZVec<STATE> &sigtrial, STATE xi, STATE beta, TPZFMatrix<STATE> &deriv) const;
    
    /// Compute the derivative of the F2 residual with respecto do sigtrial
    void GradF2SigmaTrial(const TPZVec<STATE> &sigtrial, STATE theta, STATE beta, STATE k, STATE kprev, TPZFMatrix<STATE> &deriv) const;
    
    /// Compute the point on F1 in HW Cylindrical coordinates
    void F1Cyl(STATE xi,STATE beta, TPZVec<STATE> &f1cyl) const;
    
    void F2Cyl(STATE theta,STATE beta,STATE k, TPZVec<STATE> &f2cyl) const;
    
public:
    
    void Phi(TPZVec<REAL> sigma,STATE alpha,TPZVec<STATE> &phi)const;
    
    void SurfaceParamF1(TPZVec<STATE> &sigproj, STATE &xi, STATE &beta) const;
    
    void SurfaceParamF2(TPZVec<STATE> &sigproj, STATE k, STATE &theta, STATE &beta) const;
    

    
    void ProjectF1(const TPZVec<STATE>  &sigmatrial, STATE kprev, TPZVec<STATE>  &sigproj, STATE &kproj) const;
    
    void ProjectF2(const TPZVec<STATE>  &sigmatrial, STATE kprev, TPZVec<STATE>  &sigproj, STATE &kproj) const;
    
    void ProjectRing(const TPZVec<STATE>  &sigmatrial, STATE kprev, TPZVec<STATE>  &sigproj,STATE &kproj) const;
    
    void ProjectBetaConstF2(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigproj,STATE &kproj) const;
    
    /**
	 * Imposes the specified strain tensor and returns the correspondent stress state.
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 * @param[out] sigma Resultant stress
	 */
//	virtual void ApplyStrainComputeSigma(TPZPlasticState<STATE> &plasticstate, TPZVec<STATE> &sigma);
    
//    void ApplyStrainComputeSigma(const TPZVec<STATE> &eps, STATE kprev, TPZVec<STATE> &sigma,STATE &kproj) const;
    void ApplyStrainComputeSigma(TPZVec<STATE> &epst,TPZVec<STATE> &epsp,STATE & kprev,TPZVec<STATE> &epspnext,TPZVec<STATE> &stressnext,STATE & knext) const;
    

    
    void ComputeI1(TPZVec<STATE> stress, STATE &I1)const;
    void ComputeJ2(TPZVec<STATE> stress, STATE &J2)const;
    void ApplyStrainComputeElasticStress(TPZVec<STATE> &strain,TPZVec<STATE> &stress)const;
    void ApplyStressComputeElasticStrain(TPZVec<STATE> &stress,TPZVec<STATE> &strain)const;
    
    void ProjectSigma(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigmaproj,STATE &kproj) const;
    
    void ProjectSigmaDep(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigmaproj,STATE &kproj, TPZFMatrix<STATE> &GradSigma) const;
    
    /**
     * Set of methods to verify if the tangent matrices are computed correctly
     */
    
    void TaylorCheckDistF1(const TPZVec<STATE> &sigmatrial, STATE xi, STATE beta, TPZVec<STATE> &xnorm,
                           TPZVec<STATE> &errnorm) const;
    void TaylorCheckDDistF1(const TPZVec<STATE> &sigmatrial, STATE xi, STATE beta, TPZVec<STATE> &xnorm,
                           TPZVec<STATE> &errnorm) const;
    void TaylorCheckDDistF1DSigtrial(const TPZVec<STATE> &sigmatrial, STATE xi, STATE beta, TPZVec<STATE> &xnorm,
                            TPZVec<STATE> &errnorm) const;
    void TaylorCheckDistF2(const TPZVec<STATE> &sigmatrial, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &xnorm,
                           TPZVec<STATE> &errnorm) const;
    void TaylorCheckDDistF2(const TPZVec<STATE> &sigmatrial, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &xnorm,
                            TPZVec<STATE> &errnorm) const;
    void TaylorCheckDDistF2DSigtrial(const TPZVec<STATE> &sigmatrial, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &xnorm,
                                     TPZVec<STATE> &errnorm) const;
    void TaylorCheckDF1Cart(STATE xi, STATE beta,TPZVec<STATE> &xnorm,
                            TPZVec<STATE> &errnorm) const;
    
    void TaylorCheckDF2Cart(STATE theta, STATE beta, STATE k, TPZVec<STATE> &xnorm,
                            TPZVec<STATE> &errnorm) const;
    
    /// verifies the validity of dxi/dsigtrial and dbeta/dsigtrial
    void TaylorCheckParamF1Sigtrial(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const;
    
    void TaylorCheckProjectF1(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const;
    
    /// verify D(theta,beta,k)/D(sigtrial)
    void TaylorCheckDtbkDsigtrial(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const;
    
    void TaylorCheckProjectF2(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const;
    
    void TaylorCheckProjectSigma(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const;

    
    static void ConvergenceRate(TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm, TPZVec<STATE> &convergence);
    
    static void CheckCoordinateTransformation(TPZVec<STATE> &cart);
    
public:
    void SetUp(STATE A, STATE B,STATE C, STATE D,STATE K,STATE G,STATE W,STATE R,STATE Phi,STATE N,STATE Psi);
    static void MCormicRanchSand(TPZSandlerExtended &mat);
    static void ReservoirSandstone(TPZSandlerExtended &mat);
    static void SalemLimestone(TPZSandlerExtended &mat);
    static void PreSMat(TPZSandlerExtended &mat);// em MPa
    

    STATE A() { return fA;}
    STATE B() { return fB;}
    STATE C() { return fC;}
    STATE D() { return fD;}
    STATE W() { return fW;}
    STATE K() { return fK;}
    STATE R() { return fR;}
    STATE G() { return fG;}
    STATE E() { return fE;}
    STATE N() { return fN;}
    STATE Poisson() { return fnu; }
    STATE ftol;

private:
    STATE fA,fB,fC,fD,fW,fK,fR,fG,fPhi,fN,fPsi,fE,fnu;//,fk0;

//    bool fIsonCap;
    TPZElasticResponse    fElasticResponse;
    
    
};


#endif /* defined(__PZ__pzsandlerextPV__) */
