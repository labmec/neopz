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
	
    /// Integration rule associated with the topology of the element
	typename TSHAPE::IntruleType fIntRule;
	
public:
	
	TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel);
	
	TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel, int nocreate);
	
	TPZIntelGen(TPZCompMesh &mesh, const TPZIntelGen<TSHAPE> &copy);
	
	/** @brief Constructor used to generate patch mesh... generates a map of connect index from global mesh to clone mesh */
	TPZIntelGen(TPZCompMesh &mesh,
				const TPZIntelGen<TSHAPE> &copy,
				std::map<int64_t,int64_t> & gl2lcConMap,
				std::map<int64_t,int64_t> & gl2lcElMap);
	
	TPZIntelGen();
	
	~TPZIntelGen() = default;
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override = 0;
	
	/**
	 * @brief Create a copy of the given element. The clone copy have the connect indexes
	 * mapped to the local clone connects by the given map
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t>&gl2lcElMap) const override = 0;
	
	
	virtual MElementType Type() override;


  //! Reference to the connect vector
  virtual const TPZVec<int64_t> & ConnectVec() const = 0;
  
	virtual int NConnects() const override = 0;
	
	void SetConnectIndex(int i, int64_t connectindex) override = 0;
	
	int NConnectShapeF(int connect, int order) const override = 0;
	
	virtual int Dimension() const override {
		return TSHAPE::Dimension;
	}
	
	virtual int NCornerConnects() const override{
		return TSHAPE::NCornerNodes;
	}
	
	virtual int NSideConnects(int side) const override = 0;
	
	virtual int SideConnectLocId(int node, int side) const override = 0;
	
	virtual int64_t ConnectIndex(int node) const  override;
    
	virtual void SetIntegrationRule(int ord) override;
	
	/** @brief Sets the interpolation order for the interior of the element*/
	virtual void SetInterpolationOrder(int order);
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord) override = 0;
	
	/** @brief Returns the preferred order of the polynomial along side iside*/
	virtual int PreferredSideOrder(int iside) override;
	
	/** @brief Sets the preferred interpolation order along a side */
	/** This method only updates the datastructure of the element \n
	 * In order to change the interpolation order of an element, use the method PRefine */
	virtual void SetPreferredOrder(int order) override;
	
	/** @brief Sets the interpolation order of side to order*/
	void SetSideOrder(int side, int order) override = 0;
	
	/** @brief Returns the actual interpolation order of the polynomial along the side*/
   int EffectiveSideOrder(int side) const override = 0;
	/** @brief Returns the actual interpolation order of the polynomial for a connect*/
	virtual int ConnectOrder(int connect) const;

	/** @brief Compute the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) override = 0;
	
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) override = 0;
	
	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) override;
	
	/** @brief Returns the transformation which transform a point from the side to the interior of the element */
	TPZTransform<> TransformSideToElement(int side) override;
	
	virtual const TPZIntPoints &GetIntegrationRule() const  override {
    if (this->fIntegrationRule) {
      return *fIntegrationRule;
    }
    else
    {
      return fIntRule;
    }
	}
	
	virtual TPZIntPoints &GetIntegrationRule()  override {
    if (fIntegrationRule) {
      return *fIntegrationRule;
    }
    else
    {
      return fIntRule;
    }
	}
	
	/** @brief returns the unique identifier for reading/writing objects to streams */
	public:
    virtual int ClassId() const override;

	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context) override;
};

template<class TSHAPE>
int TPZIntelGen<TSHAPE>::ClassId() const{
    return Hash("TPZIntelGen") ^ TPZInterpolatedElement::ClassId() << 1 ^ TSHAPE().ClassId() << 2;
}

#endif
