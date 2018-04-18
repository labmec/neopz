/*
 *  pzhdivfull.h
 *  PZ
 *
 *  Created by labmec on 10/23/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __PZ__pzhdvifull__
#define __PZ__pzhdvifull__

#include <iostream>



#include "pzelchdiv.h"


/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
/**
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivFull : public TPZCompElHDiv<TSHAPE> {
    
    
public:
	
	TPZCompElHDivFull(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
	TPZCompElHDivFull();
	
	virtual ~TPZCompElHDivFull();
	/** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh){
		mesh->SetAllCreateFunctionsHDivFull();
	}
	
	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);
	
	virtual MElementType Type();
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	public:
virtual int ClassId() const;

	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	virtual int NConnectShapeF(int connect, int order) const;
	virtual void SetSideOrder(int side, int order);
	virtual void IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,int64_t> > & ShapeAndVec, int pressureorder);
	virtual void FirstShapeIndex(TPZVec<int64_t> &Index);
	virtual void NShapeContinuous(TPZVec<int> &order, int &nshape );
	virtual void InitMaterialData(TPZMaterialData &data);
	virtual int NFluxShapeF() const;
	virtual void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);


    		
};


template<class TSHAPE>
int TPZCompElHDivFull<TSHAPE>::ClassId() const{
    return Hash("TPZCompElHDivFull") ^ TPZCompElHDiv<TSHAPE>::ClassId() << 1;
}

/** @brief Creates computational point element for HDivFull approximate space */
TPZCompEl *CreateHDivFullPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational linear element for HDivFull approximate space */
TPZCompEl *CreateHDivFullLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational quadrilateral element for HDivFull approximate space */
TPZCompEl *CreateHDivFullQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational triangular element for HDivFull approximate space */
TPZCompEl *CreateHDivFullTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational cube element for HDivFull approximate space */
TPZCompEl *CreateHDivFullCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational prismal element for HDivFull approximate space */
TPZCompEl *CreateHDivFullPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational pyramidal element for HDivFull approximate space */
TPZCompEl *CreateHDivFullPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational tetrahedral element for HDivFull approximate space */
TPZCompEl *CreateHDivFullTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational point element for HDivFull approximate space */
TPZCompEl *CreateHDivFullBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational linear element for HDivFull approximate space */
TPZCompEl *CreateHDivFullBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational quadrilateral element for HDivFull approximate space */
TPZCompEl *CreateHDivFullBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational triangular element for HDivFull approximate space */
TPZCompEl *CreateHDivFullBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

#endif

