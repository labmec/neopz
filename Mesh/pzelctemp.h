/**
 * @file
 * @brief Contains declaration of TPZIntelGen class which implements a generic computational element.
 */

#ifndef PZELCTEMPH
#define PZELCTEMPH

#include "pzintel.h"
#include "pzquad.h"

/**
 * @brief Implements a generic computational element. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
/**
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZIntelGen : public TPZInterpolatedElement {
	
protected:
	
		TPZManVector<long,TSHAPE::NSides> fConnectIndexes;//fazer resize qdo usar
	
	typename TSHAPE::IntruleType fIntRule;
	
public:
	
	TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel, long &index);
	
	TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel, long &index, int nocreate);
	
	TPZIntelGen(TPZCompMesh &mesh, const TPZIntelGen<TSHAPE> &copy);
	
	/** @brief Constructor used to generate patch mesh... generates a map of connect index from global mesh to clone mesh */
	TPZIntelGen(TPZCompMesh &mesh,
				const TPZIntelGen<TSHAPE> &copy,
				std::map<long,long> & gl2lcConMap,
				std::map<long,long> & gl2lcElMap);
	
	TPZIntelGen();
	
	virtual ~TPZIntelGen();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
		return new TPZIntelGen<TSHAPE> (mesh, *this);
	}
	
	/**
	 * @brief Create a copy of the given element. The clone copy have the connect indexes
	 * mapped to the local clone connects by the given map
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<long,long> & gl2lcConMap,std::map<long,long>&gl2lcElMap) const
	{
		return new TPZIntelGen<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
	
	
	virtual MElementType Type();
	
	virtual int NConnects() const {
		return fConnectIndexes.size();
	}
	
	virtual void SetConnectIndex(int i, long connectindex);
	
	virtual int NConnectShapeF(int connect) const;
	
	virtual int Dimension() const {
		return TSHAPE::Dimension;
	}
	
	virtual int NCornerConnects() const {
		return TSHAPE::NCornerNodes;
	}
	
	virtual int NSideConnects(int side) const;
	
	virtual int SideConnectLocId(int node, int side) const;
	
	virtual long ConnectIndex(int node) const;
	
	virtual void SetIntegrationRule(int ord);
	
	/** @brief Sets the interpolation order for the interior of the element*/
	virtual void SetInterpolationOrder(int order);
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord);
	
	/** @brief Returns the preferred order of the polynomial along side iside*/
	virtual int PreferredSideOrder(int iside);
	
	/** @brief Sets the preferred interpolation order along a side */
	/** This method only updates the datastructure of the element \n
	 * In order to change the interpolation order of an element, use the method PRefine */
	virtual void SetPreferredOrder(int order);
	
	/** @brief Sets the interpolation order of side to order*/
	virtual void SetSideOrder(int side, int order);
	
	/** @brief Returns the actual interpolation order of the polynomial along the side*/
    virtual int SideOrder(int side) const;
	/** @brief Returns the actual interpolation order of the polynomial for a connect*/
	virtual int ConnectOrder(int connect) const;

	/** @brief Compute the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
	
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
	
	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);
	
	/** @brief Returns the transformation which transform a point from the side to the interior of the element */
	TPZTransform TransformSideToElement(int side);
	
	virtual const TPZIntPoints &GetIntegrationRule() const {
		return fIntRule;
	}
	
	virtual TPZIntPoints &GetIntegrationRule() {
		return fIntRule;
	}
	
	/** @brief returns the unique identifier for reading/writing objects to streams */
	virtual int ClassId() const;
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
};

#endif
