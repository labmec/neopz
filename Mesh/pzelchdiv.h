//$Id: pzelctemp.h,v 1.14 2008-10-08 02:13:33 phil Exp $

// -*- c++ -*-
#ifndef PZELCHDIVHTT
#define PZELCHDIVHTT

#include "pzelctemp.h"


/// This class implements a "generic" computational element
/**
By varying the classes passed as template arguments, the complete family of computational elements are implemented
@ingroup CompElement
*/
template<class TSHAPE>
class TPZCompElHDiv : public TPZIntelGen<TSHAPE> {
	/**define the interpolation order for pressure variable*/
	int fPressureOrder;

	//metodo para adicionar vetores
	void Append(TPZFMatrix &u1, TPZFMatrix &u2, TPZFMatrix &u12);
public:

  TPZCompElHDiv(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);

  TPZCompElHDiv(TPZCompMesh &mesh, const TPZCompElHDiv<TSHAPE> &copy);

  /**
   * used to generate patch mesh... generates a map of connect index from
   * global mesh to clone mesh
   */
  TPZCompElHDiv(TPZCompMesh &mesh,
              const TPZCompElHDiv<TSHAPE> &copy,
              std::map<int,int> & gl2lcConMap,
              std::map<int,int> & gl2lcElMap);

  TPZCompElHDiv();

  virtual ~TPZCompElHDiv();

  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
    return new TPZCompElHDiv<TSHAPE> (mesh, *this);
  }

  /**
   * Create a copy of the given element. The clone copy have the connect indexes
   * mapped to the local clone connects by the given map
   * @param mesh Patch clone mesh
   * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
   * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
   */
  virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int,int> & gl2lcConMap,std::map<int,int>&gl2lcElMap) const
  {
    return new TPZCompElHDiv<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
  }


  virtual MElementType Type();

  virtual int NConnects() const;

  virtual void SetConnectIndex(int i, int connectindex);

  virtual int NConnectShapeF(int connect) const;

  virtual int Dimension() const {
    return TSHAPE::Dimension;
  }

  virtual int NCornerConnects() const {
   	  return 0;
  }

  virtual int NSideConnects(int side) const;

  virtual int SideConnectLocId(int node, int side) const;

  virtual int ConnectIndex(int node) const;

  virtual void SetIntegrationRule(int ord);
  /**Identifies the interpolation order for pressure variable*/
  virtual void SetPressureOrder(int ord);
	/**return the interpolation order to dual variable */
	int DualOrder();

  /**Identifies the interpolation order on the interior of the element*/
  virtual void GetInterpolationOrder(TPZVec<int> &ord);

  /**return the preferred order of the polynomial along side iside*/
  virtual int PreferredSideOrder(int iside);

  /**Sets the preferred interpolation order along a side
  This method only updates the datastructure of the element
  In order to change the interpolation order of an element, use the method PRefine*/
  virtual void SetPreferredOrder(int order);

  /**sets the interpolation order of side to order*/
  virtual void SetSideOrder(int side, int order);

  /**returns the actual interpolation order of the polynomial along the side*/
	virtual int SideOrder(int side) const;

	virtual int ConnectOrder(int connect) const;

  /**transform a point in the parameter space of the side into a point in the space
     of the master element*/
  //  virtual void SideParameterToElement(int side,TPZVec<REAL> &par,TPZVec<REAL> &point);

  /**transform a point in the parameter space of the master element into a point in the
     space of the side*/
  //  virtual void ElementToSideParameter(int side, TPZVec<REAL> &point, TPZVec<REAL> &par);


	  /** Initialize a material data and its attributes based on element dimension, number
	   * of state variables and material definitions
	   */
	  virtual void InitMaterialData(TPZMaterialData &data);

	  /**
	   * Compute the correspondence between the normal vectors and the shape functions
	   */
	void ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int> &shapeindex);


	//Implementacao especial pra Hdiv
	/** return the vector index  of the first index shape associate to element */
	void FirstShapeIndex(TPZVec<int> &Index);
	/** return a matrix index of the shape and vector  associate to element*/
	void IndexShapeToVec(TPZVec<int> &fVectorSide,TPZVec<std::pair<int,int> > & IndexVecShape);

  /**compute the values of the shape function of the side*/
  virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi);

  /** compute the shape functions corresponding to the dual space
   *
   */
  virtual void ShapeDual(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphi);

  void Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi);
    
    ///Compute the solution for a given variable
	virtual void Solution( TPZVec<REAL> &qsi,int var,TPZVec<REAL> &sol);
    
    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
                                 const TPZFMatrix &axes, TPZVec<REAL> &sol, TPZFMatrix &dsol);

    /**   
	 *Compute the solution using Hdiv structure
	 */
	virtual void ComputeSolution(TPZMaterialData &data);


  void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);

  //virtual void Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol);

  /** Jorge 09/06/2001
   * Returns the transformation which transform a point from the side to the interior of the element
   */
  TPZTransform TransformSideToElement(int side);

  /**
  * returns the unique identifier for reading/writing objects to streams
  */
  virtual int ClassId() const;
  /**
  Save the element data to a stream
  */
  virtual void Write(TPZStream &buf, int withclassid);

  /**
  Read the element data from a stream
  */
  virtual void Read(TPZStream &buf, void *context);

};

TPZCompEl *CreateHDivPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateHDivLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateHDivQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateHDivTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateHDivCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateHDivPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateHDivPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateHDivTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

#endif
