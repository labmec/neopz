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
	
protected:
    /// index of the bottom and top connect
    int64_t fbottom_c_index = -1, ftop_c_index = -1;
    /// sideorient of bottom and top shape functions. By default -1 for bot and +1 for top
    int fbottom_side_orient = -1, ftop_side_orient = 1;
    
    /** @brief To append vectors */
	void Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12);
    template<class TVar>
    void ComputeRequiredDataT(TPZMaterialDataT<TVar> &data,
                              TPZVec<REAL> &qsi);
    template<class TVar>
    void InitMaterialDataT(TPZMaterialDataT<TVar> &data);
public:
	    
	TPZCompElHDivCollapsed(TPZCompMesh &mesh, TPZGeoEl *gel);
	
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
	

	
	virtual int NConnects() const override;
	
	virtual void SetConnectIndex(int i, int64_t connectindex) override;
	
    
    /**
     * @brief Number of shapefunctions of the connect associated
     * @param connect connect number
     * @return number of shape functions
     */
	virtual int NConnectShapeF(int connect, int order) const override;
	
	
   /**
     * @brief return the local index for side
     **/
//	virtual int ConnectSideLocId(int connect) const override;
	
	virtual int64_t ConnectIndex(int con) const override;
    
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
    
	
	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions */
	virtual void InitMaterialData(TPZMaterialData &data) override;
    
    /**
     * @brief Destroy internally allocated data structures
     */
    virtual void CleanupMaterialData(TPZMaterialData &data) override;

    //@{
	/** @brief Compute and fill data with requested attributes */
//	void ComputeRequiredData(TPZMaterialDataT<STATE> &data,
//                             TPZVec<REAL> &qsi) override{
//        ComputeRequiredDataT(data,qsi);
//    }
//    void ComputeRequiredData(TPZMaterialDataT<CSTATE> &data,
//                             TPZVec<REAL> &qsi) override{
//        ComputeRequiredDataT(data,qsi);
//    }
    
    void ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data) override;


    //@}

	
	/** @brief Computes the values of the shape function of the side without considering the side orientation*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) override;
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
    int ClassId() const override;

	/** @brief Save the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Read the element data from a stream */
	void Read(TPZStream &buf, void *context) override;
    /** @brief Refinement along the element */
//    virtual void PRefine(int order) override;

	virtual int ConnectOrder(int connect) const override;
	
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
