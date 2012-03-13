//
//  pzvoidflux.h
//  PZ
//
//  Created by Philippe Devloo on 5/3/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//
#ifndef PZVOIDFLUX
#define PZVOIDFLUX

#include "pzmaterial.h"
#include "pzdiscgal.h"

class TPZVoidFlux : public TPZDiscontinuousGalerkin
{
    
public:
    /// Constructor
    /**
     * @param conductivity conductivity of the porous material
     */
    TPZVoidFlux(int materialid, REAL conductivity, REAL bridgesize) : TPZDiscontinuousGalerkin(materialid), fConductivity(conductivity),
        fBridgeSize(bridgesize)
    {
        
    }
    
    /// copy constructor
    TPZVoidFlux(const TPZVoidFlux &cp) : TPZDiscontinuousGalerkin(cp), fConductivity(cp.fConductivity), fBridgeSize(cp.fBridgeSize)
    {
        
    }
    
    TPZVoidFlux() : TPZDiscontinuousGalerkin(), fConductivity(0.), fBridgeSize(0.)
    {
        
    }
    
    /// Name of the class as a string
    virtual std::string Name()
    {
        return "TPZVoidFlux";
    }
    
    /// returns the integrable dimension of the material
    virtual int Dimension()
    {
        return 2;
    }
    

    
    /** Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     * @since April 10, 2007
     */
    virtual void FillDataRequirements(TPZMaterialData &data);
    
    /** Fill material data parameter with necessary requirements for the
     * ContributeInterface method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     * @since April 10, 2007
     */
    virtual void FillDataRequirementsInterface(TPZMaterialData &data);
    
    /// returns the number of state variables associated with the material
    virtual int NStateVariables();
    
    ///Metodos Contribute
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the residual vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the residual vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ef, TPZBndCond &bc);
    
    
    /// computes a contribution to stiffness matrix and load vector at one integration point
    /**
     * @param data [in] all data needed to compute the stiffness matrix
     * @param weight [in] weight of the integration point
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
                                     REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);
    
    /**
     * It computes a contribution to residual vector at one integration point
     * @param data [in]
     * @param weight [in]
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
                                     REAL weight, TPZFMatrix &ef);
    
    /**
     * It computes a contribution to stiffness matrix and load vector at one BC integration point
     * @param data [in]
     * @param weight [in]
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition object
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, 
                                       REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
    
    /**
     * It computes a contribution to residual vector at one BC integration point
     * @param data [in]
     * @param weight [in]
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition object
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix &ef,TPZBndCond &bc);
    
    /**
     * Dicontinuous galerkin materials implement contribution of discontinuous elements and interfaces.
     * Interfaces may be conservative or not conservative. It is important to agglomeration techniques
     * when using multigrid pre-conditioner. Conservative interfaces into agglomerate elements do not
     * need to be computed. However non-conservative interfaces must be computed in all multigrid levels.
     * Default is non-conservative, because of the computation of a conservative interface into an agglomerate
     * does not ruin the solution.
     * @since Feb 05, 2004
     */
    virtual int IsInterfaceConservative();
        
    /** print out the data associated with the material*/
    virtual void Print(std::ostream &out = std::cout);
    
    /**returns the variable index associated with the name*/
    virtual int VariableIndex(const std::string &name);
    
    /** returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex*/
    virtual int NSolutionVariables(int var);
    
    /**returns the solution associated with the var index based on
     * the finite element approximation*/
    virtual void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<REAL> &Solout);

    /**
     * Unique identifier for serialization purposes
     */
    virtual int ClassId() const;
    
    /**
     * Save the element data to a stream
     */
    virtual void Write(TPZStream &buf, int withclassid);
    
    /**
     * Read the element data from a stream
     */
    virtual void Read(TPZStream &buf, void *context);

    /// create another material of the same type
    virtual TPZAutoPointer<TPZMaterial> NewMaterial();
    
    /// Read data of the material from a istream (file data)
    virtual void SetData(std::istream &data);
    
protected:
    
    /// Conductivity of the porous material
    REAL fConductivity;
    
    /// Thickness of the bridges between the voids
    REAL fBridgeSize;
};


#endif