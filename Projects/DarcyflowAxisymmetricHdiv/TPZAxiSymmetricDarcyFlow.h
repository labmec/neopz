/*
 *  TPZAxiSymmetricDarcyFlow.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/11/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "pzdiscgal.h"
#include "pzfmatrix.h"
#include "pzbndcond.h"
#include "pzlog.h"

#include "tpzautopointer.h"
#include "ReservoirData.h"


#ifndef TPZDARCYFLOW
#define TPZDARCYFLOW



class TPZAxiSymmetricDarcyFlow : public TPZDiscontinuousGalerkin {
    
private:
    
    TPZAutoPointer<ReservoirData> fReservoirdata;
    
public:
    
    
    /**
     * Empty Constructor
     */
    TPZAxiSymmetricDarcyFlow();
    
    /** Creates a material object and inserts it in the vector of
     *  material pointers of the mesh.
     */
    TPZAxiSymmetricDarcyFlow(int matid);
    
    
    /** Creates a material object based on the referred object and
     *  inserts it in the vector of material pointers of the mesh.
     */
    TPZAxiSymmetricDarcyFlow(const TPZAxiSymmetricDarcyFlow &mat);
    
    /**
     * Destructor
     */
    ~TPZAxiSymmetricDarcyFlow();
    
    /** Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     * @since April 10, 2007
     */
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec);
    
    /** returns the name of the material */
    std::string Name() {
        return "TPZAxiSymmetricDarcyFlow";
    }
    
    /** returns the integrable dimension of the material */
    int Dimension() const {return 2;}
    
    /** returns the number of state variables associated with the material */
    int NStateVariables() {return 1;} // for hdiv are 3
    
    /** print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    /** returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex */
    int NSolutionVariables(int var);
    
    /** Computes the divergence over the parametric space */
    void ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi);
    
    /** returns the solution associated with the var index based on
     * the finite element approximation */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
    // Contribute Methods
    
    /** @brief Not used contribute methods */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){DebugStop();}
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){DebugStop();}
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){DebugStop();}
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){DebugStop();}
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){DebugStop();}
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ef){DebugStop();}
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){DebugStop();}
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){DebugStop();}
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * Unique identifier for serialization purposes
     */
    int ClassId() const;
    
    /**
     * Save the element data to a stream
     */
    void Write(TPZStream &buf, int withclassid);
    
    /**
     * Read the element data from a stream
     */
    void Read(TPZStream &buf, void *context);
    
    
    /**
     * Set the simulation data,
     */
    void SetReservoirData(TPZAutoPointer<ReservoirData> ResData){fReservoirdata = ResData;}
    
    /**
     * Get the simulation data,
     */
    TPZAutoPointer<ReservoirData> GetReservoirData() {return fReservoirdata;}
    
};

#endif
