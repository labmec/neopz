//
//  TPZMatElasticity2D.h
//  PZ
//
//  Created by Omar on 10/27/14.
//  Changed by Nathalia on 05/12/14
//

#ifndef __PZ__TPZMatElasticity2D__
#define __PZ__TPZMatElasticity2D__

#include <stdio.h>
#include "pzmaterial.h"
#include "pzvec.h"
#include "math.h"
#include <iostream>
#include <cmath>


/**
 * @ingroup material
 * @author Omar Duran
 * @since 10/27/2014.
 * @brief Material to solve a 2D linear elasticity
 * @brief Here is used H1.
 */




/**
 * @ingroup material
 * @brief Description of Biot's (1941) Linear Poroelastic system
 */
/**
 **@ingroup Linear Elastic Equation
 * \f$  div(T(u)) + b = 0  ==> Int{Grad(v).T(u)}dx - Int{v.gN}ds  = Int{b.v}dx  \f$ (Eq. 1)
 *
 *\f$ T(u) =  lambda*Trace(E(u)I + 2*mu*(E(u)) - \f$
 *
 *\f$ E(u) =  (1/2)(Grad(u) + Transpose(Grad(u)) \f$
 *
 *@ingroup	Diffusion equation for monophasic slightly compressible flow (e.g. oil)
 *
 *\f$ -(k/mu)*Div(Grad(p))  = d/dt{Se*p + alpha*Div(u)} (Eq. 2)  \f$
 *
 */

class TPZMatElasticity2D : public TPZMaterial {
    
protected:
    
    /** @brief Forcing vector */
    TPZManVector<REAL,2>  ff;
    
    /** @brief Elasticity modulus */
    REAL fE;
    
    /** @brief Poison coeficient */
    REAL fnu;
    
    /** @brief first Lame Parameter */
    REAL flambda;
    
    /** @brief Second Lame Parameter */
    REAL fmu;
    
    /** @brief Initial Stress */
    REAL fPreStressXX;
    REAL fPreStressXY;
    REAL fPreStressYY;
    REAL fPreStressZZ;
    REAL fPreStressXZ;
    REAL fPreStressYZ;
    
    /** @brief Uses plain stress
     * @note \f$fPlaneStress = 1\f$ => Plain stress state
     * @note \f$fPlaneStress != 1\f$ => Plain Strain state
     */
    int fPlaneStress;
    
    /** @brief Uses inclined borehole
     * @note \f$fInclinedWell = 1\f$ => Inclined Borehole Elasticity
     * @note \f$fInclinedWell != 1\f$ => Generalized Elasticity
     */
    int fInclinedWell;
    
    
    /** @brief Wellbore Direction/Azimuth angle (rad) */
    REAL falpha;
    
    /** @brief Wellbore Inclination angle (rad) */
    REAL fbeta;
    
    /** @brief InSitu Vertical and Horizontal Stresses */
    REAL fPreStressHH;
    REAL fPreStresshh;
    REAL fPreStressVV;
    
    
    /** @brief Wellbore pressure */
    REAL fPw;
    
    /** @brief Wellbore radius (m)  */
    REAL frw;
    
    /** @brief Uses Analytical Solution
     * @note \f$fAnalytics = 1\f$ => Uses Analytical Solution
     * @note \f$fAnalytics != 1\f$ => Does not use Analytical Solution for fPreStresses
     */
    int fAnalytics;
    
    
    /** @brief Uses Projection
     * @note \f$fProjection = 1\f$ =>  Uses Projection
     * @note \f$fProjection != 1\f$ => Does not use Projection
     */
    int fProjection;
    
    
    /** @brief Material Constants Sandler-DiMaggio
     * @note \fA  =>  A
     * @note \fB  =>  B
     * @note \fC  =>  C
     */
    REAL fA;
    REAL fB;
    REAL fC;
    
    
    /** @brief Material Constants Mohr-COulomb and Mogi-Coulomb
     * @note \fc         =>  cohesion angle
     * @note \ffriction  =>  friction angle
     */
    REAL fc;
    REAL ffriction;

    
public:
    TPZMatElasticity2D();
    
    /**
     * @brief Creates an elastic material with:
     * @param id material id
     * @param E elasticity modulus
     * @param nu poisson coefficient
     * @param fx forcing function \f$ -x = fx \f$
     * @param fy forcing function \f$ -y = fy \f$
     * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
     */
    TPZMatElasticity2D(int matid, REAL E, REAL nu, REAL fx, REAL fy, int plainstress = 1);
    
    TPZMatElasticity2D(int matid);
    
    TPZMatElasticity2D &operator=(const TPZMatElasticity2D &copy);
    
    virtual ~TPZMatElasticity2D();
    
    /** @brief Copy constructor */
    TPZMatElasticity2D(const TPZMatElasticity2D &cp);    
    
    virtual void Print(std::ostream & out);
    
    virtual std::string Name() { return "TPZMatElasticity2D"; }
    
    int Dimension() const {return 2;}
    
    virtual int NStateVariables();

    /**
     * @brief Set parameters of elastic material:
     * @param First  Lame Parameter Lambda
     * @param Second Lame Parameter Mu -> G
     * @param fx forcing function \f$ -x = fx \f$
     * @param fy forcing function \f$ -y = fy \f$
     * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
     */
    void SetParameters(REAL Lambda, REAL mu, REAL fx, REAL fy)
    {
        fE = (mu*(3.0*Lambda+2.0*mu))/(Lambda+mu);
        fnu = (Lambda)/(2*(Lambda+mu));
        
        flambda = Lambda;
        fmu = mu;
        ff[0] = fx;
        ff[1] = fy;
    }
    
    
    /**
     * @brief Set parameters of elastic material:
     * @param First  Lame Parameter Lambda
     * @param Second Lame Parameter Mu -> G
     * @param fx forcing function \f$ -x = 0 \f$
     * @param fy forcing function \f$ -y = 0 \f$
     */
    void SetElasticParameters(REAL Eyoung, REAL nu, REAL fbx, REAL fby)
    {
        this->SetElasticity(Eyoung,nu, fbx, fby);
    }
    
    /**
     * @brief Set parameters of elastic material:
     * @param First  Lame Parameter Lambda
     * @param Second Lame Parameter Mu -> G
     * @param fx forcing function \f$ -x = fx \f$
     * @param fy forcing function \f$ -y = fy \f$
     * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
     */
    void SetElasticity(REAL Ey, REAL nu, REAL fx, REAL fy)
    {
        fE = Ey;
        fnu = nu;
        flambda = (Ey*nu)/((1+nu)*(1-2*nu));
        fmu = Ey/(2*(1+nu));
        ff[0] = fx;
        ff[1] = fy;

    }
    
    /** @brief Set plane problem
     * planestress = 1 => Plain stress state
     * planestress != 1 => Plain Strain state
     */
    void SetfPlaneProblem(int planestress)
    {
        fPlaneStress = planestress;
    }
    
    /** @brief Set plane problem
     * planestress = 1 => Plain stress state
     * planestress != 1 => Plain Strain state
     */
    void SetPlaneStrain()
    {
        fPlaneStress = 0;
    }
    
    
    /** @brief Set Inclined Well Problem
     * @note \f$fInclinedWell = 1\f$ => Inclined Borehole Elasticity
     * @note \f$fInclinedWell != 1\f$ => Generalized Elasticity
     */
    void SetInclinedWellProblem(int inclinedwell)
    {
        fInclinedWell = inclinedwell;
        
        if (inclinedwell == 1) {
            SetPlaneStrain();
        }

    }
    
    /** @brief Set Initial Inclined Wellbore Stress */
    void SetInclinedWellborePreStress(REAL &SigmaXX, REAL &SigmaXY, REAL &SigmaYY, REAL &SigmaZZ)
    {
            // Convert the insitu stresses in local inicial stresses
            this->CalculateLocalInSituStresses(SigmaXX, SigmaXY, SigmaYY, SigmaZZ);
    }
    
    
    /** @brief Set Initial Stress */
    void SetPreStress(REAL &SigmaXX, REAL &SigmaXY, REAL &SigmaYY, REAL &SigmaZZ)
    {
        if ((fInclinedWell = 1)) {
            SetInclinedWellborePreStress(SigmaXX, SigmaXY, SigmaYY, SigmaZZ);
            
        }
    
        fPreStressXX = SigmaXX;
        fPreStressXY = SigmaXY;
        fPreStressYY = SigmaYY;
        fPreStressZZ = SigmaZZ;
        
    }
    
    
    /**
     * @brief Set parameters of a inclined wellbore:
     * @param Maximum Horizontal Stress - SigmaH
     * @param Minimum Horizontal Stress - Sigmah
     * @param Vertical Stress - SigmaV
     * @param Wellbore Direction/Azimuth angle (rad) - alpha
     * @param Wellbore Inclination angle (rad) - beta
     * @param Wellbore State - Inclined (1) or Not Inclined (!=1)
     * @param Wellbore Mud Pressure - Pw
     * @param Wellbore Radius - rw
     * @param Uses Analytical Solution as PreStress - Uses (1) or not (!=1)
     * @param Wellbore Projection - Uses projection in Post Process (1) or Not (!=1)
     */
    void SetInclinedWellboreParameters(REAL SigmaH, REAL Sigmah, REAL SigmaV, REAL alpha, REAL beta, int wellborestate, REAL Pw, REAL rw, int analytics, int projection)
    {
        
        fPreStressHH  = SigmaH;
        fPreStresshh  = Sigmah;
        fPreStressVV  = SigmaV;
        falpha        = alpha;
        fbeta         = beta;
        fInclinedWell = wellborestate;
        fPw = Pw;
        frw = rw;
        fAnalytics = analytics;
        fProjection = projection;
        
        // Sets inclined wellbore PreStress
        REAL SigXX, SigXY, SigYY, SigZZ;
        SetInclinedWellborePreStress(SigXX, SigXY, SigYY, SigZZ);
        fPreStressXX = SigXX;
        fPreStressXY = SigXY;
        fPreStressYY = SigYY;
        fPreStressZZ = SigZZ;
        
        SetInclinedWellProblem(wellborestate);
        
    }
    
    
    /** @brief Calculates the Initial Stress in Local Coordinates */
    void CalculateLocalInSituStresses(REAL &SigmaXX, REAL &SigmaXY, REAL &SigmaYY, REAL &SigmaZZ)
    {
      
        /**** Rotation Matrix ******/
        // alpha = direction/azimuth
        // beta  = wellbore inclination
       
        REAL lxx = 0., lxy = 0., lxz = 0., lyx =0., lyy = 0., lyz = 0., lzx = 0., lzy = 0., lzz = 0.;
        
        // x-diretion
        lxx = cos(falpha)*cos(fbeta);
        lxy = sin(falpha)*cos(fbeta);
        lxz = -sin(fbeta);
        // y-direction
        lyx = -sin(falpha);
        lyy = cos(falpha);
        lyz = 0;
        // z-direction
        lzx = cos(falpha)*sin(fbeta);
        lzy = sin(falpha)*sin(fbeta);
        lzz = cos(fbeta);
        
        
        // Local Inicial Stresses after Inclination
        SigmaXX      = ((lxx*lxx) * fPreStressHH) + ((lxy*lxy) * fPreStresshh) + ((lxz*lxz) * fPreStressVV);
        SigmaYY      = ((lyx*lyx) * fPreStressHH) + ((lyy*lyy) * fPreStresshh) + ((lyz*lyz) * fPreStressVV);
        SigmaZZ      = ((lzx*lzx) * fPreStressHH) + ((lzy*lzy) * fPreStresshh) + ((lzz*lzz) * fPreStressVV);
        SigmaXY      = ((lxx*lyx) * fPreStressHH) + ((lxy*lyy) * fPreStresshh) + ((lxz*lyz) * fPreStressVV);
        fPreStressYZ = ((lyx*lzx) * fPreStressHH) + ((lyy*lzy) * fPreStresshh) + ((lyz*lzz) * fPreStressVV);
        fPreStressXZ = ((lzx*lxx) * fPreStressHH) + ((lzy*lxy) * fPreStresshh) + ((lzz*lxz) * fPreStressVV);
    
        
    }

    
    // Calculates analytical solution to be used as fPreStresses
    void AnalyticalWellboreSolution(REAL &SigmaX, REAL &SigmaY, REAL &SigmaXY, REAL &SigmaZ, REAL &SigmaXZ, REAL &SigmaYZ, REAL theta, REAL R)
    {
        
        
        SigmaX   = (2.0*fPreStressXX*pow(R,4.0) - (3.0*fPreStressXX - fPreStressYY - 2.0*fPw)*pow(frw,2.0)*pow(R,2.0)*cos(2.0*theta) +
                (fPreStressXX - fPreStressYY)*pow(frw,2.0)*(3*pow(frw,2.0) - 2.0*pow(R,2.0))*cos(4.0*theta) - 4*fPreStressXY*pow(frw,2.0)*pow(R,2.0)*sin(2.0*theta) +
                6.0*fPreStressXY*pow(frw,4.0)*sin(4.0*theta) - 4.0*fPreStressXY*pow(frw,2.0)*pow(R,2.0)*sin(4.0*theta))/(2.*pow(R,4.0));
        
        SigmaY   = -(-2.0*fPreStressYY*pow(R,4.0) + (fPreStressXX - 3.0*fPreStressYY + 2.0*fPw)*pow(frw,2.0)*pow(R,2.0)*cos(2.0*theta) +
                 (fPreStressXX - fPreStressYY)*pow(frw,2.0)*(3.0*pow(frw,2.0) - 2.0*pow(R,2.0))*cos(4.0*theta) + 4.0*fPreStressXY*pow(frw,2.0)*pow(R,2.0)*sin(2.0*theta) +
                 6.0*fPreStressXY*pow(frw,4.0)*sin(4.0*theta) - 4.0*fPreStressXY*pow(frw,2.0)*pow(R,2.0)*sin(4.0*theta))/(2.*pow(R,4.0));
        
        SigmaXY  = (2.0*fPreStressXY*pow(R,4.0) + fPreStressXY*(-6.0*pow(frw,4.0) + 4.0*pow(frw,2.0)*pow(R,2.0))*cos(4.0*theta) -
                 (fPreStressXX + fPreStressYY - 2.0*fPw)*pow(frw,2.0)*pow(R,2.0)*sin(2.0*theta) + 3.0*fPreStressXX*pow(frw,4.0)*sin(4.0*theta) -
                 3.0*fPreStressYY*pow(frw,4.0)*sin(4.0*theta) - 2.0*fPreStressXX*pow(frw,2.0)*pow(R,2.0)*sin(4.0*theta) +
                 2.0*fPreStressYY*pow(frw,2.0)*pow(R,2.0)*sin(4.0*theta))/(2.*pow(R,4.0));
        
        SigmaZ   = fPreStressZZ - (2.0*fnu*pow(frw,2.0)*((fPreStressXX - fPreStressYY)*cos(2.0*theta) + 2.0*fPreStressXY*sin(2.0*theta)))/pow(R,2.0);
        
        SigmaXZ  = fPreStressXZ - (fPreStressXZ*pow(frw,2)*cos(2*theta))/pow(R,2) - (fPreStressYZ*pow(frw,2)*sin(2*theta))/pow(R,2);
        
        SigmaYZ  = fPreStressYZ + (fPreStressYZ*pow(frw,2)*cos(2*theta))/pow(R,2) - (fPreStressXZ*pow(frw,2)*sin(2*theta))/pow(R,2);

        
        
    }
    
    
    /** @brief Uses Analytical Solution
     * @note \f$fAnalytics = 1\f$ => Uses Analytical Solution
     * @note \f$fAnalytics != 1\f$ => Does not use Analytical Solution for fPreStresses
     */
    void SetAnalyticalSolution(int analytics)
    {
        fAnalytics = analytics;
        
    }
    
    
    /** @brief Uses Projection
     * @note \f$fProjection = 1\f$ =>  Uses Projection
     * @note \f$fProjection != 1\f$ => Does not use Projection
     */
    void SetProjection(int projection)
    {
        fProjection = projection;
        
    }
    
    
    // Get Initial Stress */
    void GetPreStress(REAL &SigmaXX, REAL &SigmaXY, REAL &SigmaYY, REAL &SigmaZZ)
    {
        
        SigmaXX = fPreStressXX;
        SigmaXY = fPreStressXY;
        SigmaYY = fPreStressYY;
        SigmaZZ = fPreStressZZ;
        
    }

    
    // Get Elastic Materials Parameters
    void GetElasticParameters(REAL &Ey, REAL &nu, REAL &Lambda, REAL &G)
    {
        Ey = fE;
        nu =  fnu;
        Lambda =  flambda;
        G = fmu;
    }
    
    
    // Get Wellbore angles
    void GetWellboreAngles(REAL alpha, REAL beta, int projection)
    {
        alpha = falpha;
        beta = fbeta;
        projection = fProjection;
        
    }
    
    
    /** @brief Set Material Constants Sandler-DiMaggio
     * @note \fA  =>  A
     * @note \fB  =>  B
     * @note \fC  =>  C
     */
    void SetSandlerDiMaggioParameters(REAL A, REAL B, REAL C)
    {
        fA = A;
        fB = B;
        fC = C;
        
    }
    
    
    
    /** @brief Material Constants Mohr-COulomb and Mogi-Coulomb
     * @note \fc         =>  cohesion angle
     * @note \ffriction  =>  friction angle
     */
    void SetMogiAndMohrCoulombParameters(REAL c, REAL frict)
    {
        fc        = c;
        ffriction = frict;
        
    }
    
    
    
    /** @brief Get Eyoung and Poisson
     * fE young modulus
     * fnu Poisson ratio
     */
    STATE GetEyoung() {return fE;}
    STATE GetNu() {return fnu;}
    
    /** @brief Get lame parameters
     * Lambda first lame
     * Mu Second lame
     */
    STATE GetLambda() {return flambda;}
    STATE GetMu() {return fmu;}

    
    virtual void FillDataRequirements(TPZMaterialData &data);
    
    virtual void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
//    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    virtual int VariableIndex(const std::string &name);
    
    virtual int NSolutionVariables(int var);
    
    //public:
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    virtual void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right) {
        DebugStop();
    }
    
    /**
     * Save the element data to a stream
     */
    void Write(TPZStream &buf, int withclassid);
    
    /**
     * Read the element data from a stream
     */
    void Read(TPZStream &buf, void *context);

};





#endif /* defined(__PZ__TPZMatElasticity2D__) */
