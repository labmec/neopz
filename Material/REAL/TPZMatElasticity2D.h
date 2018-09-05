//
//  TPZMatElasticity2D.h
//  PZ
//
//  Created by Omar on 10/27/14.
//
//

#ifndef __PZ__TPZMatElasticity2D__
#define __PZ__TPZMatElasticity2D__

#include <stdio.h>
#include "TPZMaterial.h"
#include "pzvec.h"
#include <iostream>


/**
 * @ingroup material
 * @brief Description Linear elastic equations
 */
/**
 **@ingroup Linear Elastic Equation
 * \f$  div(T(u)) + b = 0  ==> Int{Grad(v).T(u)}dx - Int{v.gN}ds  = Int{b.v}dx  \f$ (Eq. 1)
 *
 *\f$ T(u) =  lambda*Trace(E(u)I + 2*mu*(E(u)) - \f$
 *
 *\f$ E(u) =  (1/2)(Grad(u) + Transpose(Grad(u)) \f$
 *
 */

class TPZMatElasticity2D : public TPZMaterial {
    
protected:
    
    /** @brief Forcing vector */
    TPZManVector<STATE,2>  m_f;
    
    /** @brief Elasticity modulus */
    REAL m_E;
    
    /** @brief Poison coeficient */
    REAL m_nu;
    
    /** @brief first Lame Parameter */
    REAL m_lambda;
    
    /** @brief Second Lame Parameter */
    REAL m_mu;
    
    /** @brief Initial Stress */
    REAL m_s0_xx;
    REAL m_s0_xy;
    REAL m_s0_yy;
    REAL m_s0_zz;
    
    /** @brief plain stress directive
     */
    int m_plane_stress;
    
    
public:
    
virtual int ClassId() const;

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
    
    virtual TPZMaterial *NewMaterial()
    {
        return new TPZMatElasticity2D(*this);
    }
    
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
        m_E = (mu*(3.0*Lambda+2.0*mu))/(Lambda+mu);
        m_nu = (Lambda)/(2*(Lambda+mu));
        
        m_lambda = Lambda;
        m_mu = mu;
        m_f[0] = fx;
        m_f[1] = fy;
    }
    
    /**
     * @brief Set parameters of elastic material:
     * @param Eyoung  Young modulus
     * @param nu poisson ratio
     */
    void SetElasticParameters(REAL Eyoung, REAL nu)
    {
        this->SetElasticity(Eyoung,nu);
    }
    
    /**
     * @brief Set parameters of elastic material:
     * @param Eyoung  Young modulus
     * @param nu poisson ratio
     */
    void SetElasticity(REAL Ey, REAL nu)
    {
        m_E = Ey;
        m_nu = nu;
        m_lambda = (Ey*nu)/((1+nu)*(1-2*nu));
        m_mu = Ey/(2*(1+nu));

    }    
    
    /** @brief Set plane stress problem
     */
    void SetPlaneStress()
    {
        m_plane_stress = 1;
    }
    
    /** @brief Set plane strain problem
     */
    void SetPlaneStrain()
    {
        m_plane_stress = 0;
    }    
    
    /** @brief Set Initial Stress */
    void SetPreStress(REAL SigmaXX, REAL SigmaXY, REAL SigmaYY, REAL SigmaZZ)
    {
        m_s0_xx = SigmaXX;
        m_s0_xy = SigmaXY;
        m_s0_yy = SigmaYY;
        m_s0_zz = SigmaZZ;
    }

    /// compute the stress tensor as a function of the solution gradient
    void ComputeSigma(const TPZFMatrix<STATE> &dudx, TPZFMatrix<STATE> &sigma);
    
    // Get Elastic Materials Parameters
    void GetElasticParameters(REAL &Ey, REAL &nu, REAL &Lambda, REAL &G)
    {
        Ey = m_E;
        nu =  m_nu;
        Lambda =  m_lambda;
        G = m_mu;
    }
    
    /** @brief Get Eyoung and Poisson
     * m_E young modulus
     * m_nu Poisson ratio
     */
    STATE GetEyoung() {return m_E;}
    
    STATE GetNu() {return m_nu;}
    
    /** @brief Get lame parameters
     * Lambda first lame
     * Mu Second lame
     */
    STATE GetLambda() {return m_lambda;}
    
    STATE GetMu() {return m_mu;}

    
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
    void ContributeVec(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void ContributeVec(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);

    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    virtual int VariableIndex(const std::string &name);
    
    virtual int NSolutionVariables(int var);
    
    //public:
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    virtual void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right) {
        DebugStop();
    }
    
    /**
     * @brief Computes the error due to the difference between the interpolated flux \n
     * and the flux computed based on the derivative of the solution
     */
    virtual void Errors(TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,
                        TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
                        TPZVec<STATE> &uexact, TPZFMatrix<STATE> &duexact,
                        TPZVec<REAL> &val);

    
    /**
     * Save the element data to a stream
     */
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /**
     * Read the element data from a stream
     */
    void Read(TPZStream &buf, void *context);

};

#endif /* defined(__PZ__TPZMatElasticity2D__) */
