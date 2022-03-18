/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */

#ifndef PZELCHDIVHTT
#define PZELCHDIVHTT

#include "pzelctemp.h"
#include "TPZOneShapeRestraint.h"
#include "TPZEnumApproxFamily.h"


/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDiv : public TPZIntelGen<TSHAPE> {
	
    /// vector which defines whether the normal is outward or not
    TPZManVector<int, TSHAPE::NFacets> fSideOrient;
    
    /// Data structure which defines the restraints
    std::list<TPZOneShapeRestraint> fRestraints;

protected:
  ///! Indexes of the connects associated with the elements
  TPZManVector<int64_t,TSHAPE::NFacets+1> fConnectIndexes =
    TPZManVector<int64_t,TSHAPE::NFacets+1>(TSHAPE::NFacets+1,-1);
    /** @brief To append vectors */
	void Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12);
    
    /// Family of the HDiv space being used. Changing this will change the shape generating class
    HDivFamily fhdivfam = DefaultFamily::fHDivDefaultValue;

public:
	
    //Constructors and destructor
	TPZCompElHDiv(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam = DefaultFamily::fHDivDefaultValue);
	
	TPZCompElHDiv(TPZCompMesh &mesh, const TPZCompElHDiv<TSHAPE> &copy);
	
	/**
	 * @brief Constructor used to generate patch mesh... generates a map of connect index from
	 * global mesh to clone mesh
	 */
	TPZCompElHDiv(TPZCompMesh &mesh,
				  const TPZCompElHDiv<TSHAPE> &copy,
				  std::map<int64_t,int64_t> & gl2lcConMap,
				  std::map<int64_t,int64_t> & gl2lcElMap);
	
	TPZCompElHDiv();
	
	virtual ~TPZCompElHDiv();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const  override {
		return new TPZCompElHDiv<TSHAPE> (mesh, *this);
	}
    
    const HDivFamily GetHDivFamily() const { return fhdivfam;}
	
	/**
	 * @brief Create a copy of the given element. The clone copy have the connect indexes
	 * mapped to the local clone connects by the given map
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t>&gl2lcElMap) const override
	{
		return new TPZCompElHDiv<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
	
    /** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh) override;
	
    /** @brief Prints the relevant data of the element to the output stream */
	virtual void Print(std::ostream &out = std::cout) const override;
	

	virtual MElementType Type() override;
	
	virtual int NConnects() const override;
	
	virtual void SetConnectIndex(int i, int64_t connectindex) override;

  //! Reference to the connect vector
  inline const TPZVec<int64_t> & ConnectVec() const override{
    return fConnectIndexes;
  }
	
    /// return the first one dof restraint
    int RestrainedFace();
    
    /**
     * @brief Number of shapefunctions of the connect associated
     * @param connect connect number
     * @return number of shape functions
     */
	virtual int NConnectShapeF(int connect, int order) const override;
	
	virtual int Dimension() const  override {
		return TSHAPE::Dimension;
	}
	
	virtual int NCornerConnects() const override {
		return 0;
	}
	
	virtual int NSideConnects(int side) const override;
    
	/** 
     * @brief return the local index for connect
	 **/
	virtual int SideConnectLocId(int node, int side) const override;
    
    /** 
     * @brief return the local index for side
     **/
	virtual int ConnectSideLocId(int connect) const;
	
	virtual int64_t ConnectIndex(int con) const override;
    
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
    
    /// Return a list with the shape restraints
    virtual void ResetShapeRestraints() override
    {
        fRestraints.clear();
    }

    /**
     * @brief It returns the normal orientation of the reference element by the side.
     * Only side that has dimension larger than zero and smaller than me.
     * @param side: side of the reference elemen
     */
    virtual int GetSideOrient(int side) override;
    
    /**
     * @brief It set the normal orientation of the element by the side.
     * Only side that has dimension equal to my dimension minus one.
     * @param side: side of the reference elemen
     */
    virtual void SetSideOrient(int side, int sideorient) override;
    
	
	virtual void SetIntegrationRule(int ord) override;
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord) override;
	
	/** @brief Returns the preferred order of the polynomial along side iside*/
	virtual int PreferredSideOrder(int iside) override;
	
	/*
     * @brief Sets the preferred interpolation order along a side \n
	 * This method only updates the datastructure of the element
	 * In order to change the interpolation order of an element, use the method PRefine
	 */
	virtual void SetPreferredOrder(int order) override;
	
	/** @brief Sets the interpolation order of side to order*/
	virtual void SetSideOrder(int side, int order) override;
	
	/** @brief Returns the actual interpolation order of the polynomial along the side*/
	virtual int EffectiveSideOrder(int side) const override;
	
    /**
     * @brief return the interpolation order of the polynomial for connect
     **/
	virtual int ConnectOrder(int connect) const override;
	/**
     * @brief return the number of continuous functions 
     **/
	int NShapeContinuous(TPZVec<int> &order);
    
    /// Return the maximum order of the HDiv functions
    virtual int MaxOrder() override;
    
    /// the orientation of the face
    int SideOrient(int face)
    {
#ifdef PZDEBUG
        if (face < 0 || face >= TSHAPE::NFacets) {
            DebugStop();
        }
#endif
        return fSideOrient[face];
    }
	
	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions */
	virtual void InitMaterialData(TPZMaterialData &data) override;

	/** @brief Computes the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) override ;
	
    /**
     * @brief Computes the shape functions at the point qsi corresponding to x.
     * @param qsi point in master element coordinates
     * @param phi vector of values of shapefunctions, dimension (numshape,1)
     * @param dphi matrix of derivatives of shapefunctions in master element coordinates, dimension (dim,numshape)
     */
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) override;

    /**
     * @brief Computes the shape functions at the point qsi corresponding to x.
     * The values correspond to the function in the deformed element.    *
     * @param qsi point in master element coordinates
     * @param data fields regarding the geometric transformation and data.phi and data.divphi will be filled
     */
    void ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data) override;
    
    /** @brief Compute the solution for a given variable */
	virtual void Solution( TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override;
	
	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) override;
	
	
	/** Jorge 09/06/2001
	 * @brief Returns the transformation which transform a point from the side to the interior of the element
	 */
	TPZTransform<> TransformSideToElement(int side) override;
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
    int ClassId() const override;

	/** @brief Save the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Read the element data from a stream */
	void Read(TPZStream &buf, void *context) override;
    /** @brief Refinement along the element */
    virtual void PRefine(int order) override;
    
protected:
    //@{
    /** @brief Compute the solution using Hdiv structure */
	void ReallyComputeSolution(TPZMaterialDataT<STATE> &data) override{
        ComputeSolutionHDivT(data);
    }
    void ReallyComputeSolution(TPZMaterialDataT<CSTATE> &data) override{
        ComputeSolutionHDivT(data);
    }
    //@}
    template<class TVar>
    void ComputeSolutionHDivT(TPZMaterialDataT<TVar> &data);
    // template<class TVar>
    // void ComputeRequiredDataT(TPZMaterialDataT<TVar> &data, TPZVec<REAL>&qsi);
};

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::ClassId() const{
    return Hash("TPZCompElHDiv") ^ TPZIntelGen<TSHAPE>::ClassId() << 1;
}

#include "pzcmesh.h"

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetCreateFunctions(TPZCompMesh* mesh) {
    mesh->SetAllCreateFunctionsHDiv();
}


/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateHDivLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateHDivQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateHDivTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational cube element for HDiv approximate space */
TPZCompEl *CreateHDivCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational prismal element for HDiv approximate space */
TPZCompEl *CreateHDivPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational pyramidal element for HDiv approximate space */
TPZCompEl *CreateHDivPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational tetrahedral element for HDiv approximate space */
TPZCompEl *CreateHDivTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational point element for HDiv approximate space */
TPZCompEl *CreateHDivBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateHDivBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateHDivBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateHDivBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);


/** @} */

#endif
