/*
 *  pzmultiphysiccompel.h
 *  PZ
 *
 *  Created by Agnaldo on 9/12/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef PZMULTIPHYSICCOMPELH
#define PZMULTIPHYSICCOMPELH 

#include <iostream>

#include "pzmultiphysicselement.h"
#include "pzmaterialdata.h"

#include "pzelctemp.h"

class TPZTransform;

/** @brief class to create a compute element multiphysics */
template <class TGeometry>
class TPZMultiphysicsCompEl : public TPZMultiphysicsElement {
	
protected:
	
	/** @brief List of pointers to computational elements */
	TPZManVector<TPZCompEl *,5>		fElementVec;
	
	/** @brief Indexes of the connects of the element */
	TPZVec<int> fConnectIndexes;

public:
	/**
	 * @brief Creates a multiphysic computational element within mesh. 
	 * @param mesh mesh multiphysic where will be created the element
	 * @param gel geometric element for which the computational element will be created
	 * @param index new elemen index
	 */
	TPZMultiphysicsCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);
	
	TPZMultiphysicsCompEl();
	
	virtual ~TPZMultiphysicsCompEl();
	
	/**
	 * @brief Returns a reference to the element pointers vector
	 */
	TPZManVector<TPZCompEl *,5>   &ElementVec() {return fElementVec;}
	
	//void SetElementVec(TPZManVector<TPZCompEl *> elemVec);
	
	/**
	 * @brief Compute the map of a paramenter point in the multiphysic element to a parameter point in the super element
	 * @param trVec Transform 
	**/
	virtual void AffineTransform(TPZManVector<TPZTransform> &trVec) const;
	
	
	/**
	 * @brief Method to obtain an reference index set of multiphysics computational elements.
	 * @param cmeshVec Vector of computational meshes
	 * @param refIndexVec
	 **/
	void GetReferenceIndexVec(TPZManVector<TPZCompMesh *> cmeshVec, std::set<int> &refIndexVec);
	
	
	/** @brief Method for creating a copy of the element */
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const;
	
	/**
	 * @brief Method for creating a copy of the element in a patch mesh
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the computational elements
     */
	/**
	 * Otherwise of the previous clone function, this method don't
	 * copy entire mesh. Therefore it needs to map the connect index
	 * from the both meshes - original and patch
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
									std::map<int,int> & gl2lcConMap,
									std::map<int,int> & gl2lcElMap) const;
	
	/** @brief Returns the number of nodes of the element */
	virtual int NConnects() const;
	
	/**
	 * @brief Returns the index of the ith connectivity of the element
	 * @param i connectivity index who want knows
	 */
	virtual int ConnectIndex(int i) const ;
	
	/** @brief Dimension of the element */
	virtual int Dimension() const;
	
	/**
	 * @brief Post processing method which computes the solution for the var post processed variable.
	 * @param qsi coordinate of the point in master element space where the solution will be evaluated
	 * @param var variable which will be computed
	 * @param sol (output) solution computed at the given point
	 * @see TPZMaterial::VariableIndex
	 * @see TPZMaterial::NSolutionVariables
	 * @see TPZMaterial::Solution
	 */
	/** The var index is obtained by calling the TPZMaterial::VariableIndex method with a post processing name */
	virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<REAL> &sol);

	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi master element coordinate
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 * @param axes axes associated with the derivative of the solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes);
	
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi. \n
	 * This method will function for both volumetric and interface elements
	 * @param qsi master element coordinate of the interface element
	 * @param normal vector
	 * @param leftsol finite element solution
	 * @param dleftsol solution derivatives
	 * @param leftaxes axes associated with the left solution
	 * @param rightsol finite element solution
	 * @param drightsol solution derivatives
	 * @param rightaxes axes associated with the right solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZVec<REAL> &normal,
                    TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
					TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes);
	
	/**
	 * @brief Computes solution and its derivatives in local coordinate qsi
	 * @param qsi master element coordinate
	 * @param phi matrix containing shape functions compute in qsi point
	 * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
	 * @param axes [in] axes indicating the direction of the derivatives
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
								 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);
	
	/**
	 * @brief Set the index i to node inode
	 * @param inode node to set index
	 * @param index index to be seted
	 */
	virtual void SetConnectIndex(int inode, int index);
	
	
	/** @brief Sets create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh){
		mesh->SetAllCreateFunctionsMultiphysicElem();
	}
    
    /** @brief add an element to the datastructure */
    virtual void AddElement(TPZCompEl *cel, int meshindex)
    {
            if (fElementVec.size() <= meshindex) 
		{
			fElementVec.resize(meshindex+1);
		}
		fElementVec[meshindex] = cel;
    }
	
	/**@brief Returns referred element of this*/
	virtual TPZCompEl *ReferredElement(int mesh)
	{
	
#ifdef DEBUG
		if (fElementVec.size() <= mesh) {
			PZError << "Error at " << __PRETTY_FUNCTION__ << " index does not exist!\n";
			DebugStop();
		};
#endif
		return fElementVec[mesh];
	}


	/**
	 * @brief Sets indexes of the connects of the element
	 * @param indexes List of the connects of the element
	 */
	virtual void SetConnectIndexes(TPZVec<int> &indexes)
	{
		fConnectIndexes = indexes;
	}
		
	/**
	 * @brief Prints element data
	 * @param out Indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream &out = std::cout) const;
	
	/**
	* @brief Computes the element stiffness matrix and right hand side
	* @param ek element matrix
	* @param ef element right hand side
	*/
	virtual void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);

	/** @brief Initialize element matrix in which is computed CalcStiff */
	void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);
	
	/** 
	 * @brief Initialize a material data vector and its attributes based on element dimension, number
	 * of state variables and material definitions
	 */
	void InitMaterialData(TPZVec<TPZMaterialData > &dataVec);
	
	virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension);

};

//template <class TGeometry>
//inline void TPZMultiphysicsCompEl<TGeometry>::CreateGraphicalElement(TPZGraphMesh &, int) 
//{
//	std::cout << "TPZMultiphysicsCompEl::CreateGrafEl called\n";
//	this->Print(std::cout);
//}

/** @brief Creates computational point element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational linear element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational quadrilateral element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational triangular element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational cube element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational prismal element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational pyramidal element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational tetrahedral element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

#endif
