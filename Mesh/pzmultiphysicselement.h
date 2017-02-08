/**
 * @file
 * @brief Contains the declaration of the TPZMultiphysicsElement class. This class is abstract.
 */

#ifndef PZMULTIPHYSICSELEMENTH
#define PZMULTIPHYSICSELEMENTH 

#include <iostream>

#include "pzcompel.h"
#include "pzgeoelbc.h"
#include "pzfunction.h"

class TPZMultiphysicsInterfaceElement;

class TPZMultiphysicsElement : public TPZCompEl {
	
    /// list of restraints applied to one shape function
    std::list<TPZOneShapeRestraint> fRestraints;
    
public:
	/** @brief Default constructor */
	TPZMultiphysicsElement() : TPZCompEl()
	{
	}
	/**
	 * @brief Constructor
	 * @param mesh Multiphysics mesh where will be created the element
	 * @param ref geometric element reference
	 * @param index Index of the element created
	 */
	TPZMultiphysicsElement(TPZCompMesh &mesh, TPZGeoEl *ref, long &index) : TPZCompEl(mesh, ref, index)
	{
	}
  
  /** @brief Put a copy of the element in the referred mesh */
  TPZMultiphysicsElement(TPZCompMesh &mesh, const TPZMultiphysicsElement &copy);
  
	/** @brief Default destructor */
	virtual ~TPZMultiphysicsElement()
	{
        TPZGeoEl *reference = Reference();
        if (reference) {
            reference->ResetReference();
        }
	}
	
	virtual void AddElement(TPZCompEl *cel, long mesh) = 0;
    
    virtual void AddElement(const TPZCompElSide &cel, long mesh) = 0;
    
    virtual TPZCompEl *Element(long elindex) = 0;
	
	virtual TPZCompEl *ReferredElement(long mesh) = 0;
    
    virtual long NMeshes() = 0;
	
	virtual void SetConnectIndexes(TPZVec<long> &indexes) = 0;
	
	virtual void AffineTransform(TPZManVector<TPZTransform<> > &trVec) const = 0;
	
	virtual void InitMaterialData(TPZVec<TPZMaterialData > &dataVec) = 0;	
    
    virtual void ComputeRequiredData(TPZVec<REAL> &point, TPZVec<TPZTransform<> > &trvec, TPZVec<TPZMaterialData> &datavec);
    
	/**
	 * @brief Performs an error estimate on the elemen
	 * @param fp function pointer which computes the exact solution
	 * @param errors (output) each norm or true error of the error of the solution at each physics
	 * @param flux (input) value of the interpolated flux values
	 */
	virtual void EvaluateError(  void (*fp)(const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv),
                               TPZVec<REAL> &errors,TPZBlock<REAL> * flux );

    virtual void EvaluateError(TPZFunction<STATE> &func,
                               TPZVec<REAL> &errors);

	virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef) = 0;
	
	virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension)=0;
	
	//virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, std::set<int> dimension, std::set<int> dimension)=0;	
    
    void CreateInterfaces();
    
    bool ExistsInterface(int side);
    
    TPZMultiphysicsInterfaceElement * CreateInterface(int side);
    
    void RemoveInterfaces();
    
    void RemoveInterface(int side);
    
   	/** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
    virtual const TPZIntPoints &GetIntegrationRule() const = 0;
    
    /** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
    virtual TPZIntPoints &GetIntegrationRule() = 0;
    
    /// After adding the elements initialize the integration rule
    virtual void InitializeIntegrationRule() = 0;

    virtual int IntegrationOrder() = 0;
    
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<long> &connectindexes) const
    {
        std::cout << "To Be Implemented\n";
        DebugStop();
    }

    virtual void TransferMultiphysicsElementSolution();
    
    virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) = 0;
    
    /// Add a shape restraint (meant to fit the pyramid to restraint
    virtual void AddShapeRestraint(TPZOneShapeRestraint restraint)
    {
        fRestraints.push_back(restraint);
    }
    
    /// Return a list with the shape restraints
    virtual std::list<TPZOneShapeRestraint> GetShapeRestraints() const
    {
        return fRestraints;
    }

};

#endif
