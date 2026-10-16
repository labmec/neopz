
#ifndef PZELC_CONSTFLUX_HYBRID
#define PZELC_CONSTFLUX_HYBRID

/**
 * @brief This class implements a computational element corresponding to an Unitary Lagrange multiplier 
 * with orientation given by fSideOrient
 * @addtogroup CompElement
 * @{
 */

#include "TPZCompElDisc.h"
#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatLoadCases.h"
#include "TPZMatErrorSingleSpace.h"

class TPZCompElUnitaryLagrange : public TPZCompElDisc{

public:
    /// @brief CompElUnitaryLagrange constructor - Consists in a simplified interface element
    /// with constant flux, currently used for CompElHDivDuplConnects semi-hybridization
    /// @param mesh Computational mesh
    /// @param reference Interface GeoElement
    /// @param wrapSide Wrap neighbour CompElSide 
    /// @param lagrangeSide Lagrange neighbour CompElSide
    TPZCompElUnitaryLagrange(TPZCompMesh &mesh, TPZGeoEl *reference, TPZCompElSide &wrapSide, TPZCompElSide &lagrangeSide, bool allocateNewConnect);

    void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) override{
        CalcStiffInternal(ek,ef);
    }

    virtual int NConnects() const override;

     /**
     * @brief Number of shapefunctions of the connect associated
     * @param connect connect number
     * @return number of shape functions
     */
	virtual int NConnectShapeF(int connect, int order) const override;
	
    virtual int NSideConnects(int side) const override;

    virtual int64_t ConnectIndex(int con) const override;

    /** 
     * @brief return the local index for connect
	 **/
	virtual int SideConnectLocId(int node, int side) const override;
    
    void SetConnectIndex(int i, int64_t connectindex) override;
    
    void SetSideOrient(int sideorient){
        fSideOrient = sideorient;
    };

    int GetSideOrient(){
        return fSideOrient;
    };


protected:
	
	/** @brief It preserves index of connect associated to the element */
	TPZManVector<int64_t,2> fConnectIndexes;

    int fSideOrient = 1;

    template<class TVar>
    void CalcStiffInternal(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef);
};

#endif