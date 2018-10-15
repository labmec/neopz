/**
 * @file
 * @brief It has the declaration of the TPZMultiphysicsCompEl class.
 */

#ifndef PZMULTIPHYSICCOMPELH
#define PZMULTIPHYSICCOMPELH 

#include <iostream>

#include "pzmultiphysicselement.h"
#include "pzmaterialdata.h"

#include "pzelctemp.h"
#include "pzreducedspace.h"

template<class T>
class TPZTransform;

/** @brief class to create a compute element multiphysics */
template <class TGeometry>
class TPZMultiphysicsCompEl : public TPZMultiphysicsElement {
	
protected:
	
	/** @brief List of pointers to computational elements */
	TPZManVector<TPZCompElSide ,5>		fElementVec;
	
	/** @brief Indexes of the connects of the element */
	TPZVec<int64_t> fConnectIndexes;
    
    /// Integration rule associated with the element
    typename TGeometry::IntruleType fIntRule;
	
public:
	/**
	 * @brief Creates a multiphysic computational element within mesh. 
	 * @param mesh mesh multiphysic where will be created the element
	 * @param gel geometric element for which the computational element will be created
	 * @param index new elemen index
	 */
	TPZMultiphysicsCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	/** @brief Default constructor */
	TPZMultiphysicsCompEl();
  
	/** @brief Put a copy of the element in the referred mesh */
    TPZMultiphysicsCompEl(TPZCompMesh &mesh, const TPZMultiphysicsCompEl<TGeometry> &copy);
    /** @brief Constructor used to generate patch mesh... generates a map of connect index from global mesh to clone mesh */
	TPZMultiphysicsCompEl(TPZCompMesh &mesh,
              const TPZMultiphysicsCompEl<TGeometry> &copy,
              std::map<int64_t,int64_t> & gl2lcConMap,
              std::map<int64_t,int64_t> & gl2lcElMap);
  
	/** @brief Default destructor */
	virtual ~TPZMultiphysicsCompEl();
	
	/** @brief Returns a reference to the element pointers vector */
	TPZManVector<TPZCompElSide,5>   &ElementVec() { return fElementVec; }
	
	/**
	 * @brief Compute the map of a paramenter point in the multiphysic element to a parameter point in the super element
	 * @param trVec Transform 
	 **/
	virtual void AffineTransform(TPZVec<TPZTransform<> > &trVec) const override;

    /**
	 * @brief Performs an error estimate on the elemen
	 * @param fp function pointer which computes the exact solution
	 * @param errors (output) each norm or true error of the error of the solution at each physics
	 * @param flux (input) value of the interpolated flux values
	 */
    virtual void EvaluateError(std::function<void(const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv)> func,
                               TPZVec<REAL> &errors, bool store_error ) override;
    
    virtual void EvaluateError(TPZFunction<STATE> &func,
                               TPZVec<STATE> &errors, bool store_error) override;
    

	/**
	 * @brief Method to obtain an reference index set of multiphysics computational elements.
	 * @param cmeshVec Vector of computational meshes
	 * @param refIndexVec
	 **/
	void GetReferenceIndexVec(TPZManVector<TPZCompMesh *> cmeshVec, std::set<int64_t> &refIndexVec);

	/** @brief Method for creating a copy of the element */
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override;
	
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
									std::map<int64_t,int64_t> & gl2lcConMap,
									std::map<int64_t,int64_t> & gl2lcElMap) const override;
	
	/** @brief Returns the number of nodes of the element */
	virtual int NConnects() const override;
	
    /** @brief Return the number of meshes associated with the element */
    virtual int64_t NMeshes() override
    {
        return fElementVec.size();
    }
    
	/**
	 * @brief Returns the index of the ith connectivity of the element
	 * @param i connectivity index who want knows
	 */
	virtual int64_t ConnectIndex(int i) const override;
    
    virtual int64_t ConnectIndex(int elem, int connect) const ;
	
	/** @brief Dimension of the element */
	virtual int Dimension() const override;

    /** Returns the maximum interpolation order of all connected elements */
    virtual int IntegrationOrder() override;

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
	virtual void Integrate(int variable, TPZVec<STATE> & value) override;
	
	/** @brief Compute and fill data with requested attributes for each of the compels in fElementVec*/
    virtual void ComputeRequiredData(TPZVec<REAL> &point, TPZVec<TPZTransform<> > &trvec, TPZVec<TPZMaterialData> &datavec);

    virtual void ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &point) override
    {
        TPZMultiphysicsElement::ComputeRequiredData(data, point);
    }

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
	virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override;
	
    /**
     * @brief Compute the integral of a variable
     */
    virtual TPZVec<STATE> IntegrateSolution(int var) const override;

    
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi master element coordinate
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 * @param axes axes associated with the derivative of the solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes) override;
	
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
								 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes) override;
	
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
								 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol) override;
	
	/**
	 * @brief Set the index i to node inode
	 * @param inode node to set index
	 * @param index index to be seted
	 */
	virtual void SetConnectIndex(int inode, int64_t index) override;
	
	
	/** @brief Sets create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh) override;
    
    /** @brief add an element to the datastructure */
    virtual void AddElement(TPZCompEl *cel, int64_t meshindex) override
    {
		if (fElementVec.size() <= meshindex) 
		{
			fElementVec.resize(meshindex+1);
		}
        if (cel)
        {
            TPZGeoEl *gel = cel->Reference();
            TPZCompElSide celside(cel,gel->NSides()-1);
            fElementVec[meshindex] = celside;
        }
        else
        {
            fElementVec[meshindex] = TPZCompElSide();
        }
    }
    
    /** @brief add an element to the datastructure */
    virtual void AddElement(const TPZCompElSide &celside, int64_t meshindex) override
    {
        if (fElementVec.size() <= meshindex)
        {
            fElementVec.resize(meshindex+1);
        }
        fElementVec[meshindex] = celside;
    }
    
    virtual TPZCompEl *Element(int64_t elindex) override
    {
        return fElementVec[elindex].Element();
    }
	
	/**@brief Returns referred element of this*/
	virtual TPZCompEl *ReferredElement(int64_t mesh) override
	{
		
#ifdef PZDEBUG
		if (fElementVec.size() <= mesh) {
			PZError << "Error at " << __PRETTY_FUNCTION__ << " index does not exist!\n";
			DebugStop();
		};
#endif
		
		return fElementVec[mesh].Element();
	}
	
	
	/**
	 * @brief Sets indexes of the connects of the element
	 * @param indexes List of the connects of the element
	 */
	virtual void SetConnectIndexes(TPZVec<int64_t> &indexes) override
	{
		fConnectIndexes = indexes;
	}
	
	/**
	 * @brief Prints element data
	 * @param out Indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream &out = std::cout) const override;
	
	/**
	 * @brief Computes the element stiffness matrix and right hand side
	 * @param ek element matrix
	 * @param ef element right hand side
	 */
	virtual void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef) override;
	
    /**
     * @brief Computes the element stiffness matrix and right hand side
     * @param ek element matrix
     * @param ef element right hand side
     */
    virtual void CalcResidual(TPZElementMatrix &ef) override;
    
	/** @brief Initialize element matrix in which is computed CalcStiff */
	void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);
	
    /** @brief Initialize element matrix in which is computed CalcStiff */
    void InitializeElementMatrix(TPZElementMatrix &ef);
    
	/**
	 * @brief Initialize a material data vector and its attributes based on element dimension, number
	 * of state variables and material definitions
	 */
	void InitMaterialData(TPZVec<TPZMaterialData > &dataVec, TPZVec<int64_t> *indices = 0) override;
	
	virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) override;
	
	//virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, std::set<int> dimension, std::set<int> MaterialID);
  
    virtual void SetIntegrationRule(int int_order) override;
    
    /// After adding the elements initialize the integration rule
    virtual void InitializeIntegrationRule() override;
    
    /** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
    virtual const TPZIntPoints &GetIntegrationRule() const override;
    
    /** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
    virtual TPZIntPoints &GetIntegrationRule() override;
	
	/** @brief Return the size of the elementvec in multiphysics, if it is not multiphysics, just return 1 */
	virtual int NumberOfCompElementsInsideThisCompEl() override {
		return fElementVec.NElements();
	}	
    public:
virtual int ClassId() const override;

};


/** @brief Creates computational point element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational linear element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational quadrilateral element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational triangular element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational cube element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational prismal element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational pyramidal element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational tetrahedral element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

//--------------------- WITH MEMORY ----------------------

/** @brief Creates computational point element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsPointElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational linear element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsLinearElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational quadrilateral element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsQuadElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational triangular element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsTriangleElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational cube element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsCubeElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational prismal element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsPrismElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational pyramidal element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsPyramElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Creates computational tetrahedral element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsTetraElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

#include "pzcmesh.h"

template<class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::SetCreateFunctions(TPZCompMesh* mesh) {
    mesh->SetAllCreateFunctionsMultiphysicElem();
}

#endif
