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
    TPZGeoElSide fGeoElVolSide;

    /// element representing the flux element
    TPZCompElHDivBound2<TSHAPE> * fCelFlux;
    
public:
	    
	TPZCompElHDivSBFem(TPZCompMesh &mesh, TPZGeoEl *gel, TPZGeoElSide &gelside);

	TPZCompElHDivSBFem(TPZCompMesh &mesh, TPZGeoEl *gel);
	
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

    void ComputeSBFemVolumeHdivData(TPZMaterialDataT<STATE> &data, int64_t nshape1d);

    void AdjustAxes3D(const TPZFMatrix<REAL> &axes2D, TPZFMatrix<REAL> &axes3D, TPZFMatrix<REAL> &jac3D, TPZFMatrix<REAL> &jacinv3D, REAL detjac);

    void ExtendShapeFunctions(TPZMaterialDataT<STATE> &data2d, REAL nshape1d);

    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialDataT<STATE> &data);

    void SetCompElFlux(TPZCompElHDivBound2<TSHAPE> * cel)
    {
        fCelFlux = cel;
    }

};

template<class TSHAPE>
int TPZCompElHDivSBFem<TSHAPE>::ClassId() const{
    return Hash("TPZCompElHDivSBFem") ^ TPZIntelGen<TSHAPE>::ClassId() << 1;
}

#include "pzcmesh.h"


/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateHDivColapsedLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateHDivColapsedQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateHDivColapsedTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh);
