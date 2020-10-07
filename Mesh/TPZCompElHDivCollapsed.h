/**
 * @file
 * @brief Contains declaration of TPZCompElHDivCollapsed class which implements a generic computational element (HDiv scope).
 */

#ifndef PZELCHDIVCOLLAPSEDHTT
#define PZELCHDIVCOLLAPSEDHTT

#include "pzelchdiv.h"
#include "pzelchdivbound2.h"

/**
 * @brief This class implements a "generic" computational element of an HDiv space. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivCollapsed : public TPZCompElHDiv<TSHAPE> {
	
    /// vector which defines whether the normal is outward or not
    TPZCompElHDivBound2<TSHAPE> fBottom,fTop;
protected:
    /** @brief To append vectors */
	void Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12);

public:
	    
	TPZCompElHDivCollapsed(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
	TPZCompElHDivCollapsed(TPZCompMesh &mesh, const TPZCompElHDivCollapsed<TSHAPE> &copy);
	
	/**
	 * @brief Constructor used to generate patch mesh... generates a map of connect index from
	 * global mesh to clone mesh
	 */
	TPZCompElHDivCollapsed(TPZCompMesh &mesh,
				  const TPZCompElHDivCollapsed<TSHAPE> &copy,
				  std::map<int64_t,int64_t> & gl2lcConMap,
				  std::map<int64_t,int64_t> & gl2lcElMap);
	
	TPZCompElHDivCollapsed();
	
	virtual ~TPZCompElHDivCollapsed();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const  override {
		return new TPZCompElHDivCollapsed<TSHAPE> (mesh, *this);
	}
	
	/**
	 * @brief Create a copy of the given element. The clone copy have the connect indexes
	 * mapped to the local clone connects by the given map
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t>&gl2lcElMap) const override
	{
		return new TPZCompElHDivCollapsed<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
	
    /** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh) override;
	
    /** @brief Prints the relevant data of the element to the output stream */
	virtual void Print(std::ostream &out = std::cout) const override;
	

	
	virtual MElementType Type() override;
	
	virtual int NConnects() const override;
	
	virtual void SetConnectIndex(int i, int64_t connectindex) override;
	
    
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
//	virtual int ConnectSideLocId(int connect) const override;
	
	virtual int64_t ConnectIndex(int con) const override;
    
    /// Add a shape restraint (meant to fit the pyramid to restraint
    virtual void AddShapeRestraint(TPZOneShapeRestraint restraint) override
    {
        DebugStop();
    }
    
    /// Return a list with the shape restraints
    virtual std::list<TPZOneShapeRestraint> GetShapeRestraints() const override
    {
        std::list<TPZOneShapeRestraint> loc;
        DebugStop();
        return loc;
    }
    
    /// Return a list with the shape restraints
    virtual void ResetShapeRestraints() override
    {
        DebugStop();
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
    
	
//	virtual void SetIntegrationRule(int ord) override;
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord) override;
	
	/** @brief Returns the preferred order of the polynomial along side iside*/
	virtual int PreferredSideOrder(int iside) override;
	
	/*
     * @brief Sets the preferred interpolation order along a side \n
	 * This method only updates the datastructure of the element
	 * In order to change the interpolation order of an element, use the method PRefine
	 */
//	virtual void SetPreferredOrder(int order) override;
	
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
    
    /// Fill the polynomial order needed from the continuous shape functions
    void FillOrder(TPZVec<int> &order) const ;
	
    /// Return the maximum order??
//    virtual int MaxOrder() override;
    
    /// the orientation of the face
    int SideOrient(int face)
    {
#ifdef PZDEBUG
        if (face < 0 || face >= TSHAPE::NFacets+2) {
            DebugStop();
        }
#endif
        if(face < TSHAPE::NFacets)
            return SideOrient(face);
        else if(face == TSHAPE::NFacets)
        {
            return fBottom.GetSideOrient(face);
        }
        else
        {
            return fTop.GetSideOrient(face);
        }
    }
	
	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions */
	virtual void InitMaterialData(TPZMaterialData &data) override;
    
	/** @brief Compute and fill data with requested attributes */
	virtual void ComputeRequiredData(TPZMaterialData &data,
									 TPZVec<REAL> &qsi) override;

	/** @brief Compute the correspondence between the normal vectors and the shape functions */
	void ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int64_t> &shapeindex);
	
	/** 
	 * @brief Returns the vector index  of the first index shape associate to to each side 
	 * Special implementation to Hdiv
	 */
//	void FirstShapeIndex(TPZVec<int64_t> &Index) const;
    
	
	/** @brief Computes the values of the shape function of the side without considering the side orientation*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) override;
	
    /** @brief Computes the values of the shape function of the side NSides-1 considering the side orientation*/
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) override;
    
    /** @brief Compute the solution for a given variable */
//	virtual void Solution( TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override;
	
public:
//    virtual	void ComputeSolution(TPZVec<REAL> &qsi, TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes) override;
    
public:
    
    /** 
	 * @brief Compute shape functions based on master element in the classical FEM manne. 
	 * @param[in] qsi point in master element coordinates 
	 * @param[in] data stores all input data
	 */
//    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data) override;
	
//    void ComputeSolutionHDiv(TPZVec<REAL> &qsi, TPZMaterialData &data);
//    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
//                                 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol) override;
    
    /**
     * @brief Computes solution and its derivatives in the local coordinate qsi.
     * @param qsi master element coordinate of the interface element
     * @param normal unit normal vector
     * @param leftsol finite element solution
     * @param dleftsol solution derivatives
     * @param leftaxes axes associated with the left solution
     * @param rightsol finite element solution
     * @param drightsol solution derivatives
     * @param rightaxes axes associated with the right solution
     */
    /**
     * This method will function for both volumetric and interface elements
     */
    virtual void ComputeSolution(TPZVec<REAL> &qsi,
                                 TPZVec<REAL> &normal,
                                 TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
                                 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes) override
    {
        DebugStop();
    }

    
    /** @brief Compute the solution using Hdiv structure */
	void ComputeSolutionHDiv(TPZMaterialData &data);
	
//	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) override;
	
	
	/** Jorge 09/06/2001
	 * @brief Returns the transformation which transform a point from the side to the interior of the element
	 */
//	TPZTransform<> TransformSideToElement(int side) override;
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	public:
int ClassId() const override;

	/** @brief Save the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Read the element data from a stream */
	void Read(TPZStream &buf, void *context) override;
    /** @brief Refinement along the element */
//    virtual void PRefine(int order) override;
	
};

template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::ClassId() const{
    return Hash("TPZCompElHDivCollapsed") ^ TPZIntelGen<TSHAPE>::ClassId() << 1;
}

#include "pzcmesh.h"


/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateHDivColapsedLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateHDivColapsedQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateHDivColapsedTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);


/** @} */

#endif
