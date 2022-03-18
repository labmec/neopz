/**
 * @file
 * @brief Contains declaration of TPZCompElHDivBound2 class which implements a generic computational element (HDiv scope variant).
 */

#ifndef PZELCHDIVBOUND_H_
#define PZELCHDIVBOUND_H_

#include "pzelctemp.h"
#include "TPZOneShapeRestraint.h"
#include "TPZEnumApproxFamily.h"

/** \addtogroup CompElement */
/** @{ */
/**
 * @brief Implements a generic computational element to HDiv scope. \ref CompElement "Computational Element"
 * @author Philippe Devloo
 * @since Sep 29, 2009.
 */
/**
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivBound2 : public TPZIntelGen<TSHAPE> {
  int fSideOrient;
	
	/** @brief Method to append vectors */
	void Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12);

protected:
    ///! Indexes of the connects associated with the elements
    TPZManVector<int64_t,1> fConnectIndexes = TPZManVector<int64_t,1>(1,-1);
    
    /// Family of the HDiv space being used. Changing this will change the shape generating class
    HDivFamily fhdivfam = DefaultFamily::fHDivDefaultValue;
public:
	
	TPZCompElHDivBound2(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam = DefaultFamily::fHDivDefaultValue);
	
	TPZCompElHDivBound2(TPZCompMesh &mesh, const TPZCompElHDivBound2<TSHAPE> &copy);

	/** 
	 * @brief Constructor used to generate patch mesh...\n 
	 * Generates a map of connect index from global mesh to clone mesh
	 */
	TPZCompElHDivBound2(TPZCompMesh &mesh,
						const TPZCompElHDivBound2<TSHAPE> &copy,
						std::map<int64_t,int64_t> & gl2lcConMap,
						std::map<int64_t,int64_t> & gl2lcElMap);
	
	/** @brief Default constructor */
	TPZCompElHDivBound2();
	/** @brief Default destructor */
	virtual ~TPZCompElHDivBound2();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override {
		auto cop = new TPZCompElHDivBound2<TSHAPE> (mesh, *this);
#ifdef PZDEBUG
        if(cop->ConnectIndex(0) != this->ConnectIndex(0)) DebugStop();
#endif
        return cop;
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
		return new TPZCompElHDivBound2<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
	
	/** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh) override;
	
	virtual MElementType Type() override;
	
	virtual int NConnects() const override;
	
	virtual void SetConnectIndex(int i, int64_t connectindex) override;
	
	virtual int NConnectShapeF(int connect, int order) const override;

  inline const TPZVec<int64_t> & ConnectVec() const override{
    return fConnectIndexes;
  }
  
	virtual int Dimension() const  override {
		return TSHAPE::Dimension;
	}
	
	virtual int NCornerConnects() const override {
		return 0;
	}
	
	virtual int NSideConnects(int side) const override;
	
	virtual int SideConnectLocId(int node, int side) const override;
    
    
    virtual void SetSideOrient(int side, int sideorient) override;
    virtual int GetSideOrient(int side) override;
	
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord) override;
	
	/** @brief Returns the preferred order of the polynomial along side iside*/
	virtual int PreferredSideOrder(int iside) override;
	
	/** @brief Sets the interpolation order of side to order*/
	virtual void SetSideOrder(int side, int order) override;
	
	/** @brief Returns the actual interpolation order of the polynomial along the side*/
	virtual int EffectiveSideOrder(int side) const override;
	
	virtual int ConnectOrder(int connect) const override;

	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions */
	virtual void InitMaterialData(TPZMaterialData &data) override;
	
    
    /** @brief Compute Shape for boundary of a hdiv computational element */
    void ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data) override;
    
    /** @brief Compute the correspondence between the normal vectors and the shape functions */
    void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                      REAL &detjac, TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidx) override;
    
	/** @brief Compute the correspondence between the normal vectors and the shape functions */
	void ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int64_t> &shapeindex);
	
	/** @brief Returns the vector index  of the first index shape associate to element */
	/** Special implementation to Hdiv */
	void FirstShapeIndex(TPZVec<int64_t> &Index);
	
	/** @brief Compute the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) override;
	
	/** @brief Compute the shape function at the integration point */
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) override;
	
	/** @brief Returns a matrix index of the shape and vector  associate to element*/
	void IndexShapeToVec(TPZVec<int> &fVectorSide,TPZVec<std::pair<int,int64_t> > & IndexVecShape);
    
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

	/** @brief Returns the unique identifier for reading/writing objects to streams */
	public:
int ClassId() const override;

	/** @brief Saves the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Reads the element data from a stream */
	void Read(TPZStream &buf, void *context) override;
    
    /** @brief Prints the relevant data of the element to the output stream */
    virtual void Print(std::ostream &out) const override;
    
};

template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::ClassId() const{
    return Hash("TPZCompElHDivBound2") ^ TPZIntelGen<TSHAPE>::ClassId() << 1;
}


/** @} */

#endif /* PZELCHDIVBOUND_H_ */
