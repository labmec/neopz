/**
 * @file
 * @brief Contains declaration of TPZCompElHCurl class which implements a generic computational HCurl-conforming element
 */

#ifndef TPZCOMPELHCURL_H
#define TPZCOMPELHCURL_H

#include <pzelctemp.h>


class TPZHCurlSettings{
public:
    enum EHCurlFamily{EFullOrder = 1};
private:
    static EHCurlFamily hCurlFamily;
public:
    static void SetHCurlFamily(EHCurlFamily fam){
        hCurlFamily = fam;
    }

    static EHCurlFamily GetHCurlFamily(){
        return hCurlFamily;
    }
};

/**
 * @brief This class implements a "generic" computational HCurl-conforming element. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHCurl : public TPZIntelGen<TSHAPE> {
protected:
    /// vector describing the permutation associated with each side
    TPZManVector<int, TSHAPE::NSides - TSHAPE::NCornerNodes> fSidePermutation;

    /// Data structure which defines the restraints
    std::list<TPZOneShapeRestraint> fRestraints;

    //the hcurl vectors already taking into account the sides transformation ids, but on the master element
    TPZFNMatrix<TSHAPE::Dimension * TSHAPE::Dimension * TSHAPE::NSides,REAL> fMasterDirections;

public:
	TPZCompElHCurl(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
	TPZCompElHCurl(TPZCompMesh &mesh, const TPZCompElHCurl<TSHAPE> &copy);
	
	/**
	 * @brief Constructor used to generate patch mesh... generates a map of connect index from
	 * global mesh to clone mesh
	 */
	TPZCompElHCurl(TPZCompMesh &mesh,
				  const TPZCompElHCurl<TSHAPE> &copy,
				  std::map<int64_t,int64_t> & gl2lcConMap,
				  std::map<int64_t,int64_t> & gl2lcElMap);
	
	TPZCompElHCurl();
	
	virtual ~TPZCompElHCurl();

    /** @brief Returns the unique identifier for reading/writing objects to streams */
    int ClassId() const override;
    //since TPZCompElHCurl is an abstract class
    static int StaticClassId();
    /** @brief Set create function in TPZCompMesh to create elements of this type */
	void SetCreateFunctions(TPZCompMesh *mesh) override;

    MElementType Type() override;

	int NConnects() const override;
	/** @brief return the local index for connect */
    int SideConnectLocId(int node, int side) const override;

    int NSideConnects(int side) const override;

    int NCornerConnects() const override {
        return 0;
    }

    int64_t ConnectIndex(int con) const override;

    void SetConnectIndex(int i, int64_t connectindex) override;

    /**
    * @brief Number of shapefunctions of the connect associated
    * @param connect connect number
    * @return number of shape functions
    */
    int NConnectShapeF(int connect, int order) const override = 0;

    /**
    * @brief return the interpolation order of the polynomial for connect
    **/
    int ConnectOrder(int connect) const override;

    /** @brief Returns the actual interpolation order of the polynomial along the side*/
    int EffectiveSideOrder(int side) const override;

    /** @brief Sets the interpolation order of side to order*/
    void SetSideOrder(int side, int order) override;
protected:
    void CreateHCurlConnects(TPZCompMesh &mesh);
//
//	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const  override {
//		return new TPZCompElHCurl<TSHAPE> (mesh, *this);
//	}
//
//	/**
//	 * @brief Create a copy of the given element. The clone copy have the connect indexes
//	 * mapped to the local clone connects by the given map
//	 * @param mesh Patch clone mesh
//	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
//	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
//	 */
//	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t>&gl2lcElMap) const override
//	{
//		return new TPZCompElHCurl<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
//	}
//
//
//    /** @brief Prints the relevant data of the element to the output stream */
//	virtual void Print(std::ostream &out = std::cout) const override;
//
//
//
//	virtual MElementType Type() override;
//

//
//    /// return the first one dof restraint
//    int RestrainedFace();
//
//    /**
//     * @brief Number of shapefunctions of the connect associated
//     * @param connect connect number
//     * @return number of shape functions
//     */
//	virtual int NConnectShapeF(int connect, int order) const override;
//
//	virtual int Dimension() const  override {
//		return TSHAPE::Dimension;
//	}
//
//	virtual int NCornerConnects() const override {
//		return 0;
//	}
//	/**
//     * @brief return the number of shape for flux(just for flux)
//	 **/
//	virtual int NFluxShapeF() const;
//
//	virtual int NSideConnects(int side) const override;
//
//
//    /**
//     * @brief return the local index for side
//     **/
//	virtual int ConnectSideLocId(int connect) const;
//
//	virtual int64_t ConnectIndex(int con) const override;
//
//    /// Add a shape restraint (meant to fit the pyramid to restraint
//    virtual void AddShapeRestraint(TPZOneShapeRestraint restraint) override
//    {
//        fRestraints.push_back(restraint);
//    }
//
//    /// Return a list with the shape restraints
//    virtual std::list<TPZOneShapeRestraint> GetShapeRestraints() const override
//    {
//        return fRestraints;
//    }
//
//    /// Return a list with the shape restraints
//    virtual void ResetShapeRestraints() override
//    {
//        fRestraints.clear();
//    }
//
//    /**
//     * @brief It returns the normal orientation of the reference element by the side.
//     * Only side that has dimension larger than zero and smaller than me.
//     * @param side: side of the reference elemen
//     */
//    virtual int GetSideOrient(int side) override;
//
//    /**
//     * @brief It set the normal orientation of the element by the side.
//     * Only side that has dimension equal to my dimension minus one.
//     * @param side: side of the reference elemen
//     */
//    virtual void SetSideOrient(int side, int sideorient) override;
//
//
//	virtual void SetIntegrationRule(int ord) override;
//
//	/** @brief Identifies the interpolation order on the interior of the element*/
//	virtual void GetInterpolationOrder(TPZVec<int> &ord) override;
//
//	/** @brief Returns the preferred order of the polynomial along side iside*/
//	virtual int PreferredSideOrder(int iside) override;
//
//	/*
//     * @brief Sets the preferred interpolation order along a side \n
//	 * This method only updates the datastructure of the element
//	 * In order to change the interpolation order of an element, use the method PRefine
//	 */
//	virtual void SetPreferredOrder(int order) override;
//

//

//

//	/**
//     * @brief return the number of continuous functions
//     **/
//	int NShapeContinuous(TPZVec<int> &order);
//
//    /// Fill the polynomial order needed from the continuous shape functions
//    void FillOrder(TPZVec<int> &order) const ;
//
//    /// Return the maximum order??
//    virtual int MaxOrder() override;
//
//    /// the orientation of the face
//    int SideOrient(int face)
//    {
//#ifdef PZDEBUG
//        if (face < 0 || face >= TSHAPE::NFaces) {
//            DebugStop();
//        }
//#endif
//        return fSideOrient[face];
//    }
//
//	/** @brief Initialize a material data and its attributes based on element dimension, number
//	 * of state variables and material definitions */
//	virtual void InitMaterialData(TPZMaterialData &data) override;
//
//	/** @brief Compute and fill data with requested attributes */
//	virtual void ComputeRequiredData(TPZMaterialData &data,
//									 TPZVec<REAL> &qsi) override;
//
//	/** @brief Compute the correspondence between the normal vectors and the shape functions */
//	void ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int64_t> &shapeindex);
//
//	/**
//	 * @brief Returns the vector index  of the first index shape associate to to each side
//	 * Special implementation to Hdiv
//	 */
//	void FirstShapeIndex(TPZVec<int64_t> &Index) const;
//
//	/**
//     * @brief Returns a matrix index of the shape and vector  associate to element
//     * @param[in] VectorSide Indicates the side associated with each vector
//     * @param[out] IndexVecShape Indicates for the vector/shape function for the approximation space
//	 * @param[in] pressureorder Order of the pressure (to select shape functions?)
//	 */
//	void IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,int64_t> > & IndexVecShape, int pressureorder);
//
//	/**
//     * @brief Returns a matrix index of the shape and vector  associate to element
//     * @param[in] VectorSide Indicates the side associated with each vector
//     * @param[out] IndexVecShape Indicates for the vector/shape function for the approximation space
//	 * @param[in] pressureorder Order of the pressure (to select shape functions?)
//	 */
//	void IndexShapeToVec(TPZVec<int> &VectorSide, TPZVec<int> &bilinear, TPZVec<int> &direction, TPZVec<std::pair<int,int64_t> > & IndexVecShape, int pressureorder);
//
//	/** @brief Computes the values of the shape function of the side*/
//	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) override;
//
//	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) override;
//
//    /** @brief Compute the solution for a given variable */
//	virtual void Solution( TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override;
//
//public:
//    virtual	void ComputeSolution(TPZVec<REAL> &qsi, TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes) override;
//
//public:
//
//    /**
//	 * @brief Compute shape functions based on master element in the classical FEM manne.
//	 * @param[in] qsi point in master element coordinates
//	 * @param[in] data stores all input data
//	 */
//    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data) override;
//
//    void ComputeSolutionHDiv(TPZVec<REAL> &qsi, TPZMaterialData &data);
//    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
//                                 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol) override;
//
//    /**
//     * @brief Computes solution and its derivatives in the local coordinate qsi.
//     * @param qsi master element coordinate of the interface element
//     * @param normal unit normal vector
//     * @param leftsol finite element solution
//     * @param dleftsol solution derivatives
//     * @param leftaxes axes associated with the left solution
//     * @param rightsol finite element solution
//     * @param drightsol solution derivatives
//     * @param rightaxes axes associated with the right solution
//     */
//    /**
//     * This method will function for both volumetric and interface elements
//     */
//    virtual void ComputeSolution(TPZVec<REAL> &qsi,
//                                 TPZVec<REAL> &normal,
//                                 TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
//                                 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes) override
//    {
//        DebugStop();
//    }
//
//
//    /** @brief Compute the solution using Hdiv structure */
//	void ComputeSolutionHCurl(TPZMaterialData &data);
//
//	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) override;
//
//
//	/** Jorge 09/06/2001
//	 * @brief Returns the transformation which transform a point from the side to the interior of the element
//	 */
//	TPZTransform<> TransformSideToElement(int side) override;
//
//
//	/** @brief Save the element data to a stream */
//	void Write(TPZStream &buf, int withclassid) const override;
//
//	/** @brief Read the element data from a stream */
//	void Read(TPZStream &buf, void *context) override;
//    /** @brief Refinement along the element */
//    virtual void PRefine(int order) override;
	
};


/** @brief Creates computational linear element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational quadrilateral element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational triangular element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational cube element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational prismal element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational pyramidal element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational tetrahedral element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational point element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational linear element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational quadrilateral element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational triangular element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

TPZCompEl * CreateRefHCurlLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefHCurlQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefHCurlTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefHCurlCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefHCurlPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefHCurlPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefHCurlTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);


/** @} */

#endif
