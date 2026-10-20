#ifndef TPZSYMTENSORAPPROXCREATOR_H
#define TPZSYMTENSORAPPROXCREATOR_H

#include "TPZApproxCreator.h"

/**
 * Builds the two-space, strongly symmetric tensor approximation of main-sym.cpp.
 * Stress uses three continuous scalar components, disconnected across triangles
 * except at the macroelement center; displacement has order p-1 in each triangle
 * and order p on its skeleton. This is not an HDiv element family.
 *
 * Insert a TPZMixedSymElasticityND and its boundary materials before calling
 * CreateApproximationSpace(). The returned mesh owns the inserted materials;
 * its two atomic meshes are returned in MeshVector(), in stress/displacement order.
 * As for TPZHDivApproxCreator, the caller manages the atomic meshes' lifetime.
 * Unprepared domain elements must support refinement patterns.
 * Creation modifies the supplied geometry and must be called only once.
 */
class TPZSymTensorApproxCreator : public TPZApproxCreator {
public:
    struct MaterialIds {
        int fDomain = 1;
        int fMacroBoundary = 4;
        int fWrap = 5;
        int fInterfacePositive = 6;
        int fInterfaceNegative = 7;
        int fDisplacement = 8;
        int fInterfacePositiveLow = 9;
        int fInterfaceNegativeLow = 10;
        int fDisplacementLow = 11;
        int fTractionLow = 12;
    };

    explicit TPZSymTensorApproxCreator(TPZGeoMesh *gmesh);
    ~TPZSymTensorApproxCreator() = default;

    void SetMaterialIds(const MaterialIds &ids) {
        fIds = ids;
        fHybridizationData.fWrapMatId = ids.fWrap;
        auto &lagrange = fHybridizationData.fLagrangeData;
        lagrange.Resize(3);
        lagrange[0].fLagrangeMatId = ids.fDisplacement;
        lagrange[0].fLeftInterfaceMatId = ids.fInterfacePositive;
        lagrange[0].fRightInterfaceMatId = ids.fInterfaceNegative;
        lagrange[0].fLeftMultiplier = 1.;
        lagrange[0].fRightMultiplier = -1.;
        lagrange[1].fLagrangeMatId = ids.fTractionLow;
        lagrange[1].fLeftInterfaceMatId = ids.fInterfaceNegativeLow;
        lagrange[1].fRightInterfaceMatId = ids.fInterfacePositiveLow;
        lagrange[1].fLeftMultiplier = -1.;
        lagrange[1].fRightMultiplier = 1.;
        lagrange[2].fLagrangeMatId = ids.fDisplacementLow;
        lagrange[2].fLeftInterfaceMatId = ids.fInterfacePositiveLow;
        lagrange[2].fRightInterfaceMatId = ids.fInterfaceNegativeLow;
        lagrange[2].fLeftMultiplier = 1.;
        lagrange[2].fRightMultiplier = -1.;
    }
    const MaterialIds &GetMaterialIds() const { return fIds; }

    /// -1: single hybridization; 0..p: triple hybridization on macro boundaries.
    void SetLowOrder(int order) { fLowOrder = order; }
    int GetLowOrder() const { return fLowOrder; }

    /// Use an existing triangle fan with exactly one unmarked center per triangle.
    /// Macro boundary markers must already exist; wraps/interfaces must not exist.
    void SetGeometryIsPrepared(bool prepared) { fGeometryIsPrepared = prepared; }

    /// Only elasticity is supported. Hybridization is selected with SetLowOrder.
    void SetProbType(ProblemType prob) override;
    /// @brief Set the hybridization type. Only EStandard and EStandardSquared are supported.
    /// @param hybrid definition of the hybridization type. EStandard is a single hybridization, EStandardSquared is a triple hybridization.
    void SetHybridType(HybridizationType hybrid) override;
    /// @brief Verify if the parameters are consistent with the problem type and hybridization type. Throws an exception if not.
    void CheckSetupConsistency() override;
    TPZMultiphysicsCompMesh *CreateApproximationSpace() override;

protected:
    void GroupAndCondenseElements(TPZMultiphysicsCompMesh *mesh) override;

private:
    MaterialIds fIds;
    int fLowOrder = -1;
    bool fGeometryIsPrepared = false;
    bool fCreated = false;

    void PrepareGeometry();
    /// @brief Divide the geometric elements into triangles oriented towards the macroelement center. The macroelement center must be marked with fIds.fWrap.
    void CreateJohnsonMercier(TPZGeoMesh *gmesh);
    
    void AddWrapElements(TPZGeoMesh *gmesh);
    void AddWrapElementsTripleHybrid(TPZGeoMesh *gmesh);
    TPZCompMesh *CreateTensorSpace(TPZGeoMesh *gmesh, int porder, int porderlow);
    TPZCompMesh *CreateDisplacementSpace(TPZGeoMesh *gmesh, int porder, int porderlow);
    TPZMultiphysicsCompMesh *GenerateMultiphysicsMesh(TPZGeoMesh *gmesh, int porder, int porderlow);
    void AddInterfaceElements(TPZMultiphysicsCompMesh *mesh, int porderlow);
    std::set<int64_t> InactiveInterfaceConnects(TPZMultiphysicsCompMesh *mesh);
    void GroupElements(TPZMultiphysicsCompMesh *mesh, int porderlow);
};
#endif
