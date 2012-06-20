/**
 * @file
 * @brief Contains the declaration of the TPZMultiphysicsElement class. This class is abstract.
 */

#ifndef PZMULTIPHYSICSELEMENTH
#define PZMULTIPHYSICSELEMENTH 

#include <iostream>

#include "pzcompel.h"
#include "pzgeoelbc.h"

class TPZMultiphysicsInterfaceElement;

class TPZMultiphysicsElement : public TPZCompEl {
	
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
	TPZMultiphysicsElement(TPZCompMesh &mesh, TPZGeoEl *ref, int &index) : TPZCompEl(mesh, ref, index)
	{
	}
	/** @brief Default destructor */
	virtual ~TPZMultiphysicsElement()
	{
	}
	
	virtual void AddElement(TPZCompEl *cel, int mesh) = 0;
    
    virtual TPZCompEl *Element(int elindex) = 0;
	
	virtual TPZCompEl *ReferredElement(int mesh) = 0;
    
    virtual int NMeshes() = 0;
	
	virtual void SetConnectIndexes(TPZVec<int> &indexes) = 0;
	
	virtual void AffineTransform(TPZManVector<TPZTransform> &trVec) const = 0;
    
    virtual void ComputeRequiredData(TPZVec<REAL> &point, TPZVec<TPZTransform> &trvec, TPZVec<TPZMaterialData> &datavec);
	
	virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef) = 0;
	
	virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension)=0;
    
    void CreateInterfaces();
    
    bool ExistsInterface(int side);
    
    TPZMultiphysicsInterfaceElement * CreateInterface(int side);
    
    void RemoveInterfaces();
    
    void RemoveInterface(int side);

    virtual int IntegrationOrder() = 0;
    
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int> &connectindexes) const
    {
        std::cout << "To Be Implemented\n";
        DebugStop();
    }
    

	
};

#endif