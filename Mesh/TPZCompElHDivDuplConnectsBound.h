/**
 * @file
 * @brief Contains declaration of TPZCompElHDivDuplConnectsBound
 */

#ifndef PZELCHDIVBOUND_DUPL_CONNECTS
#define PZELCHDIVBOUND_DUPL_CONNECTS

#include "pzelchdivbound2.h"

/**
 * @brief This class implements a the HDiv-derived boundary element with two connects per facet.
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivDuplConnectsBound : public TPZCompElHDivBound2<TSHAPE> {

    /// vector which defines whether the normal is outward or not
    TPZManVector<int, TSHAPE::NFacets> fSideOrient;
    
    bool fDuplicationActive = false;

public:

    TPZCompElHDivDuplConnectsBound(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam = DefaultFamily::fHDivDefaultValue);

    virtual int NConnects() const override;

     /**
     * @brief Number of shapefunctions of the connect associated
     * @param connect connect number
     * @return number of shape functions
     */
	virtual int NConnectShapeF(int connect, int order) const override;
	

    virtual int NSideConnects(int side) const override;

    virtual int64_t ConnectIndex(int con) const override;

    virtual int SideConnectLocId(int node, int side) const override;

    void SetConnectIndex(int i, int64_t connectindex) override;    
	
	void ActiveDuplConnects(std::map<int64_t,int64_t> &fConnDuplicated);
    void InactiveDuplConnects();
};



#endif
