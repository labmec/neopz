/**
 * @file
 * @brief Contains declaration of TPZCompElHDivBound2 class which implements a generic computational element (HDiv scope variant).
 */

#ifndef PZELCHDIVBOUND_H_
#define PZELCHDIVBOUND_H_

#include "pzelctemp.h"
#include "TPZOneShapeRestraint.h"

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
    
    TPZCompElSide fneighbour;
    
    /// Restraint on a single shape function for pyramid implementation
    TPZOneShapeRestraint fRestraint;
public:
	
	TPZCompElHDivBound2(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
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
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
		return new TPZCompElHDivBound2<TSHAPE> (mesh, *this);
	}
	
	/**
	 * @brief Create a copy of the given element. The clone copy have the connect indexes
	 * mapped to the local clone connects by the given map
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t>&gl2lcElMap) const
	{
		return new TPZCompElHDivBound2<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
	
	/** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh);
	
	virtual MElementType Type();
	
	virtual int NConnects() const;
	
	virtual void SetConnectIndex(int i, int64_t connectindex);
	
	virtual int NConnectShapeF(int connect, int order) const;
	
	virtual int Dimension() const {
		return TSHAPE::Dimension;
	}
	
	virtual int NCornerConnects() const {
		return 0;
	}
	
	virtual int NSideConnects(int side) const;
	
	virtual int SideConnectLocId(int node, int side) const;
    
    
    virtual void SetSideOrient(int side, int sideorient);
    virtual int GetSideOrient(int side);
	
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord);
	
	/** @brief Returns the preferred order of the polynomial along side iside*/
	virtual int PreferredSideOrder(int iside);
	
	/** @brief Sets the interpolation order of side to order*/
	virtual void SetSideOrder(int side, int order);
	
	/** @brief Returns the actual interpolation order of the polynomial along the side*/
	virtual int EffectiveSideOrder(int side) const;
	
	virtual int ConnectOrder(int connect) const;

	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions */
	virtual void InitMaterialData(TPZMaterialData &data);
	
    
    /** @brief Compute Shape for boundary of a hdiv computational element */
    void ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data);
    
    /** @brief Compute the correspondence between the normal vectors and the shape functions */
    void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                      REAL &detjac, TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidx);
    
	/** @brief Compute the correspondence between the normal vectors and the shape functions */
	void ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int64_t> &shapeindex);
	
	/** @brief Returns the vector index  of the first index shape associate to element */
	/** Special implementation to Hdiv */
	void FirstShapeIndex(TPZVec<int64_t> &Index);
	
	/** @brief Compute the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
	
	/** @brief Compute the shape function at the integration point */
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
	
	/** @brief Returns a matrix index of the shape and vector  associate to element*/
	void IndexShapeToVec(TPZVec<int> &fVectorSide,TPZVec<std::pair<int,int64_t> > & IndexVecShape);
    
    /// Add a shape restraint (meant to fit the pyramid to restraint
    virtual void AddShapeRestraint(TPZOneShapeRestraint restraint)
    {
        if (fRestraint.IsInitialized()) {
            std::cout << "****************** Overwriting a constraint\n";
        }
        fRestraint = restraint;
    }
    

    /// Return a list with the shape restraints
    virtual std::list<TPZOneShapeRestraint> GetShapeRestraints() const
    {
        std::list<TPZOneShapeRestraint> loc;
        if (fRestraint.IsInitialized()) {
            loc.push_back(fRestraint);
        }
        return loc;
    }

    /// Return a list with the shape restraints
    virtual void ResetShapeRestraints()
    {
        fRestraint = TPZOneShapeRestraint();
    }

	/** @brief Returns the unique identifier for reading/writing objects to streams */
	public:
virtual int ClassId() const;

	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
    
    /** @brief Prints the relevant data of the element to the output stream */
    virtual void Print(std::ostream &out) const;

};

template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::ClassId() const{
    return Hash("TPZCompElHDivBound2") ^ TPZIntelGen<TSHAPE>::ClassId() << 1;
}

/** @brief Creates computational point element for HDiv approximate space */
TPZCompEl *CreateRefHDivBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateRefHDivBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateRefHDivBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateRefHDivBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @} */

#endif /* PZELCHDIVBOUND_H_ */
