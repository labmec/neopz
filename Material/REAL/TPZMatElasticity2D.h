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
#include "pzmaterial.h"
#include "pzvec.h"
#include <iostream>


/**
 * @ingroup material
 * @author Omar Duran
 * @since 10/27/2014.
 * @brief Material to solve a 2D linear elasticcity
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
    TPZVec<REAL>  ff;
    
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
    
    /** @brief Uses plain stress
     * @note \f$fPlaneStress = 1\f$ => Plain stress state
     * @note \f$fPlaneStress != 1\f$ => Plain Strain state
     */
    int fPlaneStress;
    
    
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
    void SetElasticParameters(REAL Eyoung, REAL nu)
    {
        this->SetElasticity(Eyoung,nu);
    }
    
    /**
     * @brief Set parameters of elastic material:
     * @param First  Lame Parameter Lambda
     * @param Second Lame Parameter Mu -> G
     * @param fx forcing function \f$ -x = fx \f$
     * @param fy forcing function \f$ -y = fy \f$
     * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
     */
    void SetElasticity(REAL Ey, REAL nu)
    {
        fE = Ey;
        fnu = nu;
        flambda = (Ey*nu)/((1+nu)*(1-2*nu));
        fmu = Ey/(2*(1+nu));

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
    
    /** @brief Set Initial Stress */
    void SetPreStress(REAL SigmaXX, REAL SigmaXY, REAL SigmaYY, REAL SigmaZZ)
    {
        fPreStressXX = SigmaXX;
	fPreStressXY = SigmaXY;
	fPreStressYY = SigmaYY;
	fPreStressZZ = SigmaZZ;
    }    
    
    
    // Get Elastic Materials Parameters
    void GetElasticParameters(REAL &Ey, REAL &nu, REAL &Lambda, REAL &G)
    {
        Ey = fE;
        nu =  fnu;
        Lambda =  flambda;
        G = fmu;
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
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
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
