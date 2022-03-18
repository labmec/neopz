#ifndef TPZCOMPELH1_H
#define TPZCOMPELH1_H


#include "pzelctemp.h"
#include "TPZEnumApproxFamily.h"

/**
   @brief TPZCompElH1 implements a H1-conforming approximation space.
*/
template<class TSHAPE>
class TPZCompElH1  : public TPZIntelGen<TSHAPE>{
protected:
  ///! Indexes of the connects associated with the elements
  TPZManVector<int64_t,TSHAPE::NSides> fConnectIndexes =
    TPZManVector<int64_t,TSHAPE::NSides>(TSHAPE::NSides,-1);
    
  /// Family of the HDiv space being used. Changing this will change the shape generating class
  H1Family fh1fam = DefaultFamily::fH1DefaultValue;
public:

  TPZCompElH1() = default;

  TPZCompElH1(TPZCompMesh &mesh, TPZGeoEl *gel, const H1Family h1fam = DefaultFamily::fH1DefaultValue);
	
  TPZCompElH1(TPZCompMesh &mesh, TPZGeoEl *gel, int nocreate, const H1Family h1fam = DefaultFamily::fH1DefaultValue);
	
  TPZCompElH1(TPZCompMesh &mesh, const TPZCompElH1<TSHAPE> &copy);

  TPZCompElH1(const TPZCompElH1 &) = default;
  TPZCompElH1(TPZCompElH1 &&) = default;
  ~TPZCompElH1();
  TPZCompElH1 & operator =(const TPZCompElH1 &) = default;
  TPZCompElH1 & operator =(TPZCompElH1 &&) = default;
  /** @brief Constructor used to generate patch mesh... generates a map of connect index from global mesh to clone mesh */
  TPZCompElH1(TPZCompMesh &mesh,
              const TPZCompElH1<TSHAPE> &copy,
              std::map<int64_t,int64_t> & gl2lcConMap,
              std::map<int64_t,int64_t> & gl2lcElMap);
  
  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const  override {
    return new TPZCompElH1<TSHAPE> (mesh, *this);
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
		return new TPZCompElH1<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
  
    /**
     * @brief Initialize a material data and its attributes based on element dimension, number
     * of state variables and material definitions
     */
virtual void InitMaterialData(TPZMaterialData &data) override;

  /// computes the shape functions in the master element AND its derivatives
  void ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data) override;

  inline int NConnects() const override{
		return TSHAPE::NSides;
	}
  int NSideConnects(int side) const override;
  
  int NConnectShapeF(int connect, int order) const override;

  void SetConnectIndex(int i, int64_t connectindex) override;

  void SetSideOrder(int side, int order) override;

  int EffectiveSideOrder(int side) const override;

  int SideConnectLocId(int node, int side) const override;

  void GetInterpolationOrder(TPZVec<int> &ord) override;

  //! Reference to the connect vector
  inline const TPZVec<int64_t> & ConnectVec() const override{
    return fConnectIndexes;
  }

  	/** @brief Compute the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) override;
	
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) override;
};


#endif /* TPZCOMPELH1_H */
