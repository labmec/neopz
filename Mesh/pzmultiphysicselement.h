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
	TPZMultiphysicsElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index) : TPZCompEl(mesh, ref, index)
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
	
	virtual void AddElement(TPZCompEl *cel, int64_t mesh) = 0;
    
    virtual void AddElement(const TPZCompElSide &cel, int64_t mesh) = 0;
    
    virtual TPZCompEl *Element(int64_t elindex) = 0;
	
	virtual TPZCompEl *ReferredElement(int64_t mesh) = 0;
    
    virtual int64_t NMeshes() = 0;
	
	virtual void SetConnectIndexes(TPZVec<int64_t> &indexes) = 0;
	
	virtual void AffineTransform(TPZVec<TPZTransform<> > &trVec) const = 0;
	
	virtual void InitMaterialData(TPZVec<TPZMaterialData > &dataVec, TPZVec<int64_t> *indices = 0) = 0;
    
    virtual void ComputeRequiredData(TPZVec<REAL> &point, TPZVec<TPZTransform<> > &trvec, TPZVec<TPZMaterialData> &datavec, TPZVec<int64_t> *indices = 0);
    
    /** @brief Compute and fill data with requested attributes */
    virtual void ComputeRequiredData(TPZMaterialData &data,
                                     TPZVec<REAL> &qsi)
    {
        
    }
	/**
	 * @brief Performs an error estimate on the elemen
	 * @param fp function pointer which computes the exact solution
	 * @param errors (output) each norm or true error of the error of the solution at each physics
	 * @param flux (input) value of the interpolated flux values
	 */
    virtual void EvaluateError( std::function<void(const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv)> func,
                               TPZVec<REAL> &errors, bool store_error ) override;

    virtual void EvaluateError(TPZFunction<STATE> &func,
                               TPZVec<STATE> &errors, bool store_error);

	virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef) override = 0 ;
	
	virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) override =0 ;
	
	//virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, std::set<int> dimension, std::set<int> dimension)=0;	
    
    void CreateInterfaces();
    
    bool ExistsInterface(int side);
    
    TPZMultiphysicsInterfaceElement * CreateInterface(int side);
    
    void RemoveInterfaces();
    
    void RemoveInterface(int side);
    
   	/** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
    virtual const TPZIntPoints &GetIntegrationRule() const override = 0 ;
    
    /** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
    virtual TPZIntPoints &GetIntegrationRule() = 0;
    
    /// After adding the elements initialize the integration rule
    virtual void InitializeIntegrationRule() override = 0 ;

    virtual int IntegrationOrder() = 0;
    
    virtual int ComputeIntegrationOrder() const override;
    
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override
    {
        std::cout << "To Be Implemented\n";
        DebugStop();
    }

    virtual void TransferMultiphysicsElementSolution() override;
    
    virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override  = 0 ;
    
    /// Add a shape restraint (meant to fit the pyramid to restraint
    virtual void AddShapeRestraint(TPZOneShapeRestraint restraint) override
    {
        fRestraints.push_back(restraint);
    }
    
    /// Return a list with the shape restraints
    virtual std::list<TPZOneShapeRestraint> GetShapeRestraints() const override
    {
        return fRestraints;
    }

};

#endif
