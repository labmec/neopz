//
//  TPZMonoPhaseWell.h
//  PZ
//
//  Created by omar duran on 25/05/2015.
//
//

#ifndef __PZ__TPZMonoPhaseWell__
#define __PZ__TPZMonoPhaseWell__

#include <stdio.h>
#include "TPZMaterial.h"


/**
 * @ingroup material
 * @brief Implements a Multiphycis mixed formulation for monophasic wel
 */
class TPZMonoPhaseWell : public TPZMaterial {

public:
    
    /**
     * @name Constructors and Destructors methods
     * @{
     */
    
    /** @brief Constructor for given convection */
    TPZMonoPhaseWell(int id);
    
    /** @brief Destructor */
    ~TPZMonoPhaseWell();
    
    /** @brief Copy constructor */
    TPZMonoPhaseWell(TPZMonoPhaseWell & copy);
    
    /** @} */
    
private:
    
    /**
     * @name Well atributes and methods
     * @{
     */
    
    /** @brief Well cross sectional area */
    REAL fAp;
    
    /** @brief Well diameter */
    REAL fd;
    
    /** @brief Pipe roughness */
    REAL fepsilon;
    
    /** @brief Pipe inclination */
    REAL ftheta;
    
    /** @brief atribute for derivative computations */
    REAL fdelta;
    
    /** @brief time step */
    REAL fdt;
    
    /** @brief Gravity constant */
    REAL fg;
    
    /** @brief time step */
    bool  fNextStep;
    
    /** @brief friction factor using Haaland friction factor because it is a explicit expression of f */
    void friction(REAL &f, REAL P, REAL w, REAL &dfdw, REAL &dfdP);
    
    
    /** @brief Fluid density  */
    void Rho(REAL &rho, REAL P, REAL &drhodw, REAL &drhodP);
    
    
    /** @brief Fluid density  */
    void Mu(REAL &mu, REAL P, REAL &dmudw, REAL &dmudP);
    
    /** @} */
    
    /** @brief Returns the name of the material */
    std::string Name() { return "TPZMonoPhaseWell";}
    
    /** @brief Returns the integrable dimension of the material */
    int Dimension() const {return 1;}
    
    /** @brief Returns the number of state variables associated with the material */
    int NStateVariables() {return 1;}
    
    /** @brief returns the minimum number of load cases for this material */
    int MinimumNumberofLoadCases()
    {
        return 1;
    }
    
    /**
     * @name Contribute methods (weak formulation)
     * @{
     */
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * to multiphysics simulation.
     * @param datavec [in]  stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 18, 2011
     */
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /** @} */
    
    /** @name Contribute methods not used here
     * @{
     */
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);

    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }

    
    /**
     * @brief It computes a contribution to the residual vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the residual vector
     * @since April 16, 2007
     */
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
        DebugStop();
    }

    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }

    /** @} */
    
    
    /**
     * @brief Fill material data parameter with necessary requirements for the
     * @since April 10, 2007
     */
    /**
     * Contribute method. Here, in base class, all requirements are considered as necessary.
     * Each derived class may optimize performance by selecting only the necessary data.
     */
    void FillDataRequirements(TPZMaterialData &data);
    
    /**
     * @brief Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered as necessary.
     * Each derived class may optimize performance by selecting only the necessary data.
     */
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    /**
     * @brief Fill material data parameter with necessary requirements for the
     * ContributeBC method. Here, in base class, all requirements are considered as necessary.
     * Each derived class may optimize performance by selecting only the necessary data.
     */
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    /** @brief Print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** @brief Returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);

	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);

    
    /** @{
     * @name Save and Load methods
     */
    
    /** @brief Unique identifier for serialization purposes */
    public:
virtual int ClassId() const;

    
    /** @brief Saves the element data to a stream */
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /** @brief Reads the element data from a stream */
    virtual void Read(TPZStream &buf, void *context);
    
    /** @} */
    
public:
    
    /** @{
     * @name Set and Get methods
     */
    
    /** @brief Well cross sectional area */
    void SetAp(REAL A){ fAp = A;}
    REAL GetAp(){ return fAp;}
    
    /** @brief Well diameter */
    void Setd(REAL d){ fd = d;}
    REAL Getd(){ return fd;}
    
    /** @brief Pipe roughness */
    void Setepsilon(REAL epsilon){ fepsilon = epsilon;}
    REAL Getepsilon(){ return fepsilon;}
    
    /** @brief Pipe inclination */
    void Settheta(REAL theta){ ftheta = theta;}
    REAL Gettheta(){ return ftheta;}
    
    /** @brief atribute for derivative computations */
    void Setdelta(REAL delta){ fdelta = delta;}
    REAL Getdelta(){ return fdelta;}
    
    /** @brief time step */
    void Setdt(REAL dt){ fdt = dt;}
    REAL Getdt(){ return fdt;}
    
    /** @brief time step */
    void SetNextStep(bool nextstep) {fNextStep = nextstep;}
    bool GetNextStep() {return fNextStep;}
    
    /** @} */
    
    
};


#endif /* defined(__PZ__TPZMonoPhaseWell__) */
