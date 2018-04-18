/**
 * @file
 * @author Philippe Devloo
 * @since 5/3/2011.
 */

#ifndef PZVOIDFLUX
#define PZVOIDFLUX

#include "TPZMaterial.h"
#include "pzdiscgal.h"

class TPZVoidFlux : public TPZDiscontinuousGalerkin
{
    
public:
    /// Constructor
    /**
	 * @param materialid id of the material
     * @param conductivity conductivity of the porous material
	 * @param bridgesize size
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
    virtual int Dimension() const
    {
        return 2;
    }
    

    
    /** 
	 * @brief Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     * @since April 10, 2007
     */
    virtual void FillDataRequirements(TPZMaterialData &data);
    
    /**
	 * @brief Fill material data parameter with necessary requirements for the
     * ContributeInterface method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     * @since April 10, 2007
     */
    virtual void FillDataRequirementsInterface(TPZMaterialData &data);
    
    /** @brief returns the number of state variables associated with the material */
    virtual int NStateVariables();
    
    /**
	 * @name Contribute Methods
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
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * @brief It computes a contribution to the residual vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the residual vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    
    /**
	 * @brief computes a contribution to stiffness matrix and load vector at one integration point
     * @param data [in] all data needed to compute the stiffness matrix
     * @param dataleft [in] data needed to compute the stiffness matrix from left element
     * @param dataright [in] data needed to compute the stiffness matrix from right element
     * @param weight [in] weight of the integration point
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
                                     REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief It computes a contribution to residual vector at one integration point
     * @param data [in]
     * @param dataleft [in]
     * @param dataright [in]
     * @param weight [in]
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
                                     REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
     * @param data [in]
     * @param dataleft [in]
     * @param weight [in]
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition object
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, 
                                       REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /**
     * @brief It computes a contribution to residual vector at one BC integration point
     * @param data [in]
     * @param dataleft [in]
     * @param weight [in]
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition object
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);

    /** @} */
	
    /**
     * @brief Dicontinuous galerkin materials implement contribution of discontinuous elements and interfaces.
     * Interfaces may be conservative or not conservative. It is important to agglomeration techniques
     * when using multigrid pre-conditioner. Conservative interfaces into agglomerate elements do not
     * need to be computed. However non-conservative interfaces must be computed in all multigrid levels.
     * Default is non-conservative, because of the computation of a conservative interface into an agglomerate
     * does not ruin the solution.
     * @since Feb 05, 2004
     */
    virtual int IsInterfaceConservative();
        
    /** @brief print out the data associated with the material*/
    virtual void Print(std::ostream &out = std::cout);
    
    /**returns the variable index associated with the name*/
    virtual int VariableIndex(const std::string &name);
    
    /** @brief returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex*/
    virtual int NSolutionVariables(int var);
    
    /**
	 * @brief returns the solution associated with the var index based on
     * the finite element approximation*/
    virtual void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<STATE> &Solout);

    /** @brief Unique identifier for serialization purposes */
    public:
virtual int ClassId() const;

    
    /** @brief Save the element data to a stream */
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /** @brief Read the element data from a stream */
    virtual void Read(TPZStream &buf, void *context);

    /** @brief create another material of the same type */
    virtual TPZMaterial * NewMaterial();
    
    /** @brief Read data of the material from a istream (file data) */
    virtual void SetData(std::istream &data);
    
protected:
    
    /** @brief Conductivity of the porous material */
    REAL fConductivity;
    
    /** @brief Thickness of the bridges between the voids */
    REAL fBridgeSize;
};


#endif