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
#include "TPZPlasticCriterion.h"

class TPZSandlerExtended : public TPZPlasticCriterion {
public:

    enum {
        NYield = 2
    };

    /// Constructor, with all parameters which define the Sandler DiMaggio model
    TPZSandlerExtended(STATE A, STATE B, STATE C, STATE D, STATE K, STATE G, STATE W, STATE R, STATE Phi, STATE N, STATE Psi, STATE kappa_0);
    /// Copy constructor
    TPZSandlerExtended(const TPZSandlerExtended & copy);
    /// Empty constructor
    TPZSandlerExtended();

    TPZSandlerExtended & operator=(const TPZSandlerExtended & source) {
        ftol = source.ftol;
        fA = source.fA;
        fB = source.fB;
        fC = source.fC;
        fD = source.fD;
        fK = source.fK;
        fG = source.fG;
        fW = source.fW;
        fR = source.fR;
        fPhi = source.fPhi;
        fN = source.fN;
        fPsi = source.fPsi;
        fE = source.fE;
        fnu = source.fnu;
        fkappa_0 = source.fkappa_0;
        fElasticResponse = source.fElasticResponse;

        return *this;
    }

    /// Desctructor
    virtual ~TPZSandlerExtended();

    STATE GetX(STATE k);

    /// Compute k as a function of epsp using Newton iterations
    void Firstk(STATE &epsp, STATE &k) const;

    TPZElasticResponse GetElasticResponse();

    void SetElasticResponse(const TPZElasticResponse &ER);

    virtual TPZElasticResponse GetElasticResponse() const;

    STATE GetR();

    virtual void YieldFunction(const TPZVec<STATE> &sigma, STATE kprev, TPZVec<STATE> &yield) const override;

    virtual int GetNYield() const override {
        return as_integer(NYield);
    }

    virtual void Print(std::ostream &out) const override;
    
    template<class T>
    T F(const T x) const;
    
    template<class T>
    T DF(const T x) const;

    STATE GetF(STATE x) const;
    
    virtual int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    
private:
    /// The function which defines the plastic surface


    /// Auxiliary function for Associating the position of the cap with the damage variable
    template<class T>
    T X(const T k) const;
    
    /// compute the damage variable as a function of the X function
    template<class T>
    T EpsEqX(T X) const;
    /// Compute the damage variable as a function of the position of the cap k
    template<class T>
    T EpsEqk(const T k) const;
    /// Compute the residual of the equation which defines the update of the damage variable
    template<class T>
    T ResLF2(const TPZVec<T> &pt, T theta, T beta, T k, STATE kprev) const;
    
    /// Compute the residual of the equation which defines the update of the damage variable
    template<class T>
    T ResLF2IJ(const TPZVec<T> &sigtrIJ, T theta, T k, STATE kprev) const;
    
    /// Compute the normal function to the failure surface based on a reference point (I1_ref,f1(I1_ref))
    STATE NormalToF1(STATE I1, STATE I1_ref) const;
    
    /// Compute the residual of the equation which defines the update of the damage variable
    STATE ResLF1(const TPZVec<STATE> &sigtrial, const TPZVec<STATE> &sigproj, const STATE k, const STATE kprev) const;
    /// Compute the distance of sigtrial to the point on the yield surface
    STATE DistF1(const TPZVec<STATE> &pt, const STATE xi, const STATE beta) const;
    /// Compute the distance of sigtrial to the point on the cap
    STATE DistF2(const TPZVec<STATE> &pt, const STATE theta, const STATE beta, const STATE k) const;

    /// Compute the distance considering the sigtrial is given as a funcion of I1, sqJ2
    STATE DistF2IJ(const TPZVec<STATE> &sigtrialIJ, STATE theta, STATE k) const;


    /// Compute the derivative of the distance function to the yield surface as a function of xi and beta
    void DDistFunc1(const TPZVec<STATE> &pt, STATE xi, STATE beta, TPZFMatrix<STATE> &ddistf1) const;
    
    /// Compute the derivative of the distance function to the failure function and the result of Residue 1 (failure)
    template<class T>
    void Res1(const TPZVec<T> &trial_stress, T i1, T beta, T k, T kprev, TPZVec<T> & residue_1) const;
    
    /// Compute the derivative of the distance function to the cap function and the result of Residue 2 (Cap)
    template<class T>
    void Res2(const TPZVec<T> &trial_stress, T theta, T beta, T k, T kprev, TPZVec<T> & residue_2) const;
    
    /// Compute the derivative of the distance function to the covertex cap function and the result of covertex Residue
    template<class T>
    void Res2CoVertex(const TPZVec<T> &trial_stress, T beta, T k, T kprev, TPZVec<T> & residue_covertex) const;
    
    /// Compute the derivative of the distance function to the vertex cap function and the result of vertex Residue
    template<class T>
    void Res2Vertex(const TPZVec<T> &trial_stress, T k, T kprev, T & residue_vertex) const;
    
    /// Compute the value of the equation which determines the orthogonality of the projection
    template<class T>
    void DDistF2IJ(TPZVec<T> &sigtrialIJ, T theta, T L, STATE Lprev, TPZVec<T> &ddistf2) const;


    /// Compute the second derivative of the distance as a function of xi and beta
    void D2DistFunc1(const TPZVec<STATE> &pt, STATE xi, STATE beta, TPZFMatrix<STATE> &d2distf1) const;

    /// Compute the jacobian function of the f1 (failure) distance as a function of i1, beta and k
    void Jacobianf1(const TPZVec<STATE> &trial_stress, STATE i1, STATE beta, STATE k, TPZFMatrix<STATE> &jacobianf1)const;
    
    /// Compute the jacobian function of the f2 (cap) distance as a function of theta, beta and k
    void Jacobianf2(const TPZVec<STATE> &trial_stress, STATE theta, STATE beta, STATE k, TPZFMatrix<STATE> &jacobianf2)const;
    
    /// Compute the jacobian function of the f2 (cap) distance as a function of beta and k
    void Jacobianf2CoVertex(const TPZVec<STATE> &trial_stress, STATE beta, STATE k, TPZFMatrix<STATE> &jacobianf2_covertex)const;
    
    /// Compute the jacobian function of the vertex on f2 (cap) distance as a function of k
    void Jacobianf2Vertex(const TPZVec<STATE> &trial_stress, STATE k, STATE &jacobianf2_vertex)const;

    /// Compute the jacobian of the distance function to the cap vertex function and the result of Vertex residue (Cap)
    void JacobianVertex(const TPZVec<STATE> &trial_stress, STATE k, STATE &jacobian_vertex)const;
    
    /// Compute the jacobian of the distance function to the cap covertex function and the result of Covertex Residue (Cap intersection with Failure)
    void JacobianCoVertex(const TPZVec<STATE> &trial_stress, STATE beta, STATE k, TPZFMatrix<STATE> &jacobian_covertex)const;

    /// Compute the derivative of the equation which determines the evolution of k
    // the derivative are given in terms of theta, beta and k
    void DResLF2(const TPZVec<STATE> &pt, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &dresl) const;
    /// Compute the derivative of the equation which determines the evolution of k
    // the derivative are given in terms of k
    STATE DResLF1(const TPZVec<STATE> &sigtrial, const TPZVec<STATE> &sigproj, const STATE k, const STATE kprev) const;
    
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
    void F1Cyl(STATE xi, STATE beta, TPZVec<STATE> &f1cyl) const;

    void F2Cyl(const STATE theta, const STATE beta, const STATE k, TPZVec<STATE> &f2cyl) const;

public:
    
    /// Compute initial damage variable from the given principal stress state
    REAL InitialDamage(const TPZVec<REAL> &stress_p) const;

    void Phi(TPZVec<REAL> sigma, STATE alpha, TPZVec<STATE> &phi)const;

    void SurfaceParamF1(TPZVec<STATE> &sigproj, STATE &xi, STATE &beta) const;

    void SurfaceParamF2(const TPZVec<STATE> &sigproj, const STATE k, STATE &theta, STATE &beta) const;

    STATE NormalFunctionToF1(STATE & I1, STATE & k) const;

    void ProjectApex(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const;
    
    void ProjectF1(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const;

    void ProjectF2(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const;
    
    void ProjectCapVertex(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const;
    
    void ProjectCapCoVertex(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const;

    void ProjectCoVertex(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const;
    
    void ProjectVertex(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const;
    
    void ProjectRing(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const;

    void ProjectBetaConstF2(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj) const;

    /**
     * Imposes the specified strain tensor and returns the correspondent stress state.
     *
     * @param[in] epsTotal Imposed total strain tensor
     * @param[out] sigma Resultant stress
     */
    //	virtual void ApplyStrainComputeSigma(TPZPlasticState<STATE> &plasticstate, TPZVec<STATE> &sigma);

    //    void ApplyStrainComputeSigma(const TPZVec<STATE> &eps, STATE kprev, TPZVec<STATE> &sigma,STATE &kproj) const;
    void ApplyStrainComputeSigma(TPZVec<STATE> &epst, TPZVec<STATE> &epsp, STATE & kprev, TPZVec<STATE> &epspnext, TPZVec<STATE> &stressnext, STATE & knext, TPZFMatrix<REAL> * tangent = NULL) const;



    void ComputeI1(TPZVec<STATE> stress, STATE &I1)const;
    void ComputeJ2(TPZVec<STATE> stress, STATE &J2)const;
    void ApplyStrainComputeElasticStress(TPZVec<STATE> &strain, TPZVec<STATE> &stress)const;
    void ApplyStressComputeElasticStrain(TPZVec<STATE> &stress, TPZVec<STATE> &strain)const;

    void ProjectSigma(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigmaproj, STATE &kproj, int &m_type, TPZFMatrix<REAL> * gradient = NULL) const;

    void ProjectSigmaDep(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigmaproj, STATE &kproj, TPZFMatrix<STATE> &GradSigma) const;

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
    void TaylorCheckDF1Cart(STATE xi, STATE beta, TPZVec<STATE> &xnorm,
            TPZVec<STATE> &errnorm) const;

    void TaylorCheckDF2Cart(STATE theta, STATE beta, STATE k, TPZVec<STATE> &xnorm,
            TPZVec<STATE> &errnorm) const;

    /// verifies the validity of dxi/dsigtrial and dbeta/dsigtrial
    void TaylorCheckParamF1Sigtrial(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const;

    void TaylorCheckProjectF1(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const;

    /// verify D(theta,beta,k)/D(sigtrial)
    void TaylorCheckDtbkDsigtrial(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const;

    void TaylorCheckProjectF2(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const;

    static void ConvergenceRate(TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm, TPZVec<STATE> &convergence);

    static void CheckCoordinateTransformation(TPZVec<STATE> &cart);
    
    /// Compute the derivative of the projected stresses respect to trial stresses (tangent) over the cap
    void ComputeCapTangent(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj, TPZFMatrix<REAL> * gradient) const;
    
    /// Compute the derivative of the projected stresses respect to trial stresses (tangent) over the cap
    void ComputeCapVertexTangent(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj, TPZFMatrix<REAL> * gradient) const;

    /// Compute the derivative of the projected stresses respect to trial stresses (tangent) over the cap
    void ComputeCapCoVertexTangent(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj, TPZFMatrix<REAL> * gradient) const;
    
    /// Compute the derivative of the projected stresses respect to trial stresses (tangent) over the failure
    void ComputeFailureTangent(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &projected_stress, STATE &kproj, TPZFMatrix<REAL> * gradient) const;
    
    /// Compute the approximation rate for the derivative of the projected stresses respect to trial stresses
    void TaylorCheckProjectSigma(const TPZVec<STATE> &trial_stress, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const;
    

public:
    
    // Set up charecteristic parameters
    void SetUp(STATE A, STATE B, STATE C, STATE D, STATE K, STATE G, STATE W, STATE R, STATE Phi, STATE N, STATE Psi);
    
    // Set up the initial damage
    void SetInitialDamage(STATE kappa_0);
    
    static void MCormicRanchSand(TPZSandlerExtended &mat);
    static void ReservoirSandstone(TPZSandlerExtended &mat);
    static void SalemLimestone(TPZSandlerExtended &mat);
    static void PreSMat(TPZSandlerExtended &mat); // em MPa

    STATE A() {
        return fA;
    }
    
    void SetA(STATE A)
    {
        fA = A;
    }

    STATE B() {
        return fB;
    }
    
    void SetB(STATE B)
    {
        fB = B;
    }

    STATE C() {
        return fC;
    }
    
    void SetC(STATE C)
    {
        fC = C;
    }

    STATE D() {
        return fD;
    }

    STATE W() {
        return fW;
    }

    STATE K() {
        return fK;
    }

    STATE R() {
        return fR;
    }

    STATE G() {
        return fG;
    }

    STATE E() {
        return fE;
    }

    STATE N() {
        return fN;
    }

    STATE Poisson() {
        return fnu;
    }
    
    STATE InitialDamage() {
        return fkappa_0;
    }
    
    STATE Apex() const {
        STATE apex = log(fA/fC)/fB;
        return apex;
    }
    
    STATE X_0() const {
        STATE X_0 = this->X(fkappa_0);
        return X_0;
    }
    
    STATE CPerturbation() const {
        STATE CK = fE/(3.0*(1.0 - 2.0 *fnu));
        STATE C_per = (fD*fC)/(3.0*CK);
        return C_per;
    }
    
    STATE ftol;

private:

    STATE fA, fB, fC, fD, fW, fK, fR, fG, fPhi, fN, fPsi, fE, fnu, fkappa_0;

    //    bool fIsonCap;
    TPZElasticResponse fElasticResponse;


};


#endif /* defined(__PZ__pzsandlerextPV__) */
