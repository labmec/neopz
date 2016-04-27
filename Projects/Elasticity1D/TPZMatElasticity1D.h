//
//  TPZMatElasticity1D.h
//  PZ
//
//  Created by Nathalia on 4/15/16.
//
//

#ifndef TPZMatElasticity1D_hpp
#define TPZMatElasticity1D_hpp

#include <stdio.h>
#include "pzmaterial.h"
#include "pzvec.h"
#include <iostream>


/**
 * @ingroup material
 * @author Nathalia
 * @since 4/15/2016.
 * @brief Material to solve a 1D linear elasticity
 * @brief Here is used H1.
 */




/**
 * @ingroup material
 * @brief Description of Linear Poroelastic
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

class TPZMatElasticity1D : public TPZMaterial {
    
protected:
    
    /** @brief Body forces vector */
    TPZVec<REAL>  fb;
    
    /** @brief Mass density */
    REAL frho;
    
    /** @brief first Lame Parameter */
    REAL flambda;
    
    /** @brief Second Lame Parameter */
    REAL fmu;

    
    /** @brief Uses line stress
     * @note \f$fLineStrain = 1\f$ => Line strain state
     * @note \f$fLineStrain != 1\f$ => Line stress state
     */
    int fLineStrain;
    
    
public:
    
    TPZMatElasticity1D();
    
    /**
     * @brief Creates an elastic material with:
     * @param id material id
     * @param Lambda first Lamé
     * @param Mu second Lamé
     * @param b body force function \f$ b = rho * g * s \f$
     * @param LineStrain \f$ LineStrain = 1 \f$ indicates use of linestrain
     */
    TPZMatElasticity1D(int matid, REAL Lambda, REAL Mu, REAL b, int LineStrain = 1);
    
    /**
     * @brief Creates an elastic material with:
     * @param id material id
     */
    TPZMatElasticity1D(int matid);
    
    TPZMatElasticity1D &operator=(const TPZMatElasticity1D &copy);
    
    ~TPZMatElasticity1D();
    
    /** @brief Copy constructor */
    TPZMatElasticity1D(const TPZMatElasticity1D &cp);
    
    
    /**
     * @brief Set parameters of elastic material:
     * @param First  Lame Parameter Lambda
     * @param Second Lame Parameter Mu -> G
     * @param b body force function \f$ b = rho * g * s \f$
     * @param LineStrain \f$ LineStrain = 1 \f$ indicates use of "linestrain"
     */
    void SetParameters(REAL Lambda, REAL mu, REAL b, int LineStrain)
    {
        fb.resize(1);
        flambda = Lambda;
        fmu = mu;
        fb[0] = b;
        fLineStrain = LineStrain;
    
    }
    
    
    /**
     * @brief Get parameters of elastic material:
     * @param First  Lame Parameter Lambda
     * @param Second Lame Parameter Mu -> G
     * @param b body force function \f$ b = rho * g * s \f$
     * @param LineStrain \f$ LineStrain = 1 \f$ indicates use of "linestrain"
     */
    void GetParameters(REAL &Lambda, REAL &mu, REAL &b, int &LineStrain)
    {
        Lambda = flambda;
        mu =  fmu;
        b =  fb[0];
        LineStrain = fLineStrain;
    }
    

    
    std::string Name() { return "TPZMatElasticity1D"; }
    
    int Dimension() const {return 1;}
    
    virtual int NStateVariables();

    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    
    void FillDataRequirements(TPZMaterialData &data);
    
    
    void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data);
    
    
    void Print(std::ostream & out);
    
    
    virtual int VariableIndex(const std::string &name);
    
    
    /**
     * Save the element data to a stream
     */
    void Write(TPZStream &buf, int withclassid);
    
    
    /**
      * Read the element data from a stream
      */
    void Read(TPZStream &buf, void *context);
    
    
    
    virtual int NSolutionVariables(int var);
    
    
    //public:
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
  
    
    
    virtual void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right) {
        DebugStop();
    }
    
    
    
};

#endif /* TPZMatElasticity1D_hpp */
