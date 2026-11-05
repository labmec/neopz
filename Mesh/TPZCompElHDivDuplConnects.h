/**
 * @file
 * @brief Contains declaration of TPZCompElHDivDuplConnects
 */

#ifndef PZELCHDIV_DUPL_CONNECTS
#define PZELCHDIV_DUPL_CONNECTS

#include "pzelchdiv.h"
#include "TPZMaterialData.h"
/**
 * @brief This class implements an HDiv-derived element containing two connects per facet. This approach
 * is usefull for semi-hybridization and applications where we want to separate certain functions in a facet.
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivDuplConnects : public TPZCompElHDiv<TSHAPE> {

    bool fDuplicationActive = false;

public:

    TPZCompElHDivDuplConnects(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam = DefaultFamily::fHDivDefaultValue);

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

    void ActiveDuplConnects(std::map<int64_t,int64_t> &fConnDuplicated);
    
    void InactiveDuplConnects();

    /** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions */
	virtual void InitMaterialData(TPZMaterialData &data) override;

};


/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateHDivDuplConnectsBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
TPZCompEl *CreateHDivDuplConnectsBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
TPZCompEl *CreateHDivDuplConnectsBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
TPZCompEl *CreateHDivDuplConnectsBoundTriangEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);

/** @brief Creates computational Linear element for HDiv approximate space */
TPZCompEl *CreateHDivDuplConnectsLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateHDivDuplConnectsQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateHDivDuplConnectsTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational cube element for HDiv approximate space */
TPZCompEl *CreateHDivDuplConnectsCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational tetrahedral element for HDiv approximate space */
TPZCompEl *CreateHDivDuplConnectsTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);


#endif
