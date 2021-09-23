/**
 * @file
 * @brief Contains declaration of TPZCompElHDivCollapsed class which implements a generic computational element (HDiv scope).
 */

#pragma once

#include "pzelchdiv.h"
#include "pzelchdivbound2.h"
#include "TPZCompElHDivCollapsed.h"

/**
 * @brief This class implements a "generic" computational element of an HDiv space. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivSBFem : public TPZCompElHDivCollapsed<TSHAPE> {

    /// geometric element representing the collapsed volume
    TPZGeoElSide fGelVolSide;

    /// vector which defines whether the normal is outward or not
    TPZManVector<int, TSHAPE::NFacets> fSideOrient;
    
public:
	    
	TPZCompElHDivSBFem(TPZCompMesh &mesh, TPZGeoEl *gel, TPZGeoElSide &gelside, int64_t &index);

	TPZCompElHDivSBFem(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
	TPZCompElHDivSBFem(TPZCompMesh &mesh, const TPZCompElHDivSBFem<TSHAPE> &copy);
	
	TPZCompElHDivSBFem();
	
	virtual ~TPZCompElHDivSBFem();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const  override {
		return new TPZCompElHDivSBFem<TSHAPE> (mesh, *this);
	}

    /** @brief Returns the unique identifier for reading/writing objects to streams */
    int ClassId() const override;

    void SetGelVolumeSide(TPZGeoElSide &gelside);

    TPZGeoElSide & GetGelVolumeSide();

    virtual void ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi) override;

    void ComputeDeformedDirections(TPZMaterialDataT<STATE> &data);

    void HDivCollapsedDirections(TPZMaterialDataT<STATE> &data, int64_t nshape1d);

    void AdjustAxes3D(const TPZFMatrix<REAL> &axes2D, TPZFMatrix<REAL> &axes3D, TPZFMatrix<REAL> &jac3D, TPZFMatrix<REAL> &jacinv3D, REAL detjac);

    void ExtendShapeFunctions(TPZMaterialDataT<STATE> &data2d, REAL nshape1d);

    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialDataT<STATE> &data);

};

template<class TSHAPE>
int TPZCompElHDivSBFem<TSHAPE>::ClassId() const{
    return Hash("TPZCompElHDivSBFem") ^ TPZIntelGen<TSHAPE>::ClassId() << 1;
}

#include "pzcmesh.h"


/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateHDivColapsedLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateHDivColapsedQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateHDivColapsedTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
