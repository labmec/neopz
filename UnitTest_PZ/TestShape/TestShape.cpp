//
// Created by Gustavo on 04/03/2020.
//

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "TPZNullMaterial.h"
#include "TPZMaterialDataT.h"
#include "pzgeotriangle.h"
#include "pzgeoquad.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "pzgeopyramid.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.testshape");
#endif

#include <catch2/catch.hpp>
template<class TGeo>
void AddSampleElement(TPZGeoMesh& gmesh);

template<class TGeo>
void CheckDivergenceOnInternalConnect();

TEST_CASE("internal_connect_divergence_test", "[shape_test]") {
  CheckDivergenceOnInternalConnect<pzgeom::TPZGeoTriangle>();
  CheckDivergenceOnInternalConnect<pzgeom::TPZGeoQuad>();
  CheckDivergenceOnInternalConnect<pzgeom::TPZGeoTetrahedra>();
  CheckDivergenceOnInternalConnect<pzgeom::TPZGeoPrism>();
  CheckDivergenceOnInternalConnect<pzgeom::TPZGeoCube>();
  // Pyramid test fails since its space is not complete
  // CheckDivergenceOnInternalConnect<pzgeom::TPZGeoPyramid>();
}

template<class TGeo>
void AddSampleElement(TPZGeoMesh& gmesh) {
    std::cout << "Creating example mesh.\n";

    const int matId = 1;
    TPZManVector<REAL, 3> lowerCorner(3, 0);
    TPZManVector<REAL, 3> size(3, 1);
    TGeo::InsertExampleElement(gmesh, matId, lowerCorner, size);
    gmesh.BuildConnectivity();
    gmesh.SetDimension(gmesh.Element(0)->Dimension());
}

template<class TGeo>
void CheckDivergenceOnInternalConnect() {

    std::cout << "Starting divergence test on " << TGeo::TypeName() << " element.\n";
    // Create element mesh
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    AddSampleElement<TGeo>(*gmesh);
    int dim = gmesh->Dimension();
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(2);
    cmesh->SetDimModel(dim);

    auto* mat = new TPZNullMaterial<STATE>(1);
    mat->SetDimension(dim);
    cmesh->InsertMaterialObject(mat);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();
    cmesh->InitializeBlock();

    const auto fluxEl = dynamic_cast<TPZInterpolatedElement*>(cmesh->Element(0));
    if (!fluxEl) DebugStop();
    if (fluxEl->Material()->Id() != 1) DebugStop();

    // Get flux geometric element
    const auto gel = fluxEl->Reference();
    if (!gel) DebugStop();

    // Initialize material requirements
    TPZMaterialDataT<STATE> elData;
    fluxEl->InitMaterialData(elData);

    // Get last connect, which is the one that contains internal shape functions
    TPZConnect& con = fluxEl->Connect(fluxEl->NConnects() - 1);

    // Create integration rule on volume side
    const int pOrderIntRule = fluxEl->EffectiveSideOrder(gel->NSides() - 1) * 2;
    TPZIntPoints* intRule = gel->CreateSideIntegrationRule(gel->NSides() - 1, pOrderIntRule);

    // Get shape function indexes to be iterated
    const int nInternalPhi = con.NShape();
    const int firstInternalPhi = fluxEl->NShapeF() - nInternalPhi;

    // Create matrix to store the integration results
    TPZFMatrix<REAL> integrationResult(nInternalPhi, 1, 0);

    // Start divergence integration
    const int npts = intRule->NPoints();
    TPZManVector<REAL, 3> xi(dim, 0);
    REAL w;
    for (auto ipt = 0; ipt < npts; ipt++) {
        intRule->Point(ipt, xi, w);

        fluxEl->ComputeRequiredData(elData, xi);
//        elData.ComputeFunctionDivergence();

        for (int iPhi = 0; iPhi < nInternalPhi; iPhi++) {
            integrationResult(iPhi, 0) += w * fabs(elData.detjac) * elData.divphi(iPhi + firstInternalPhi, 0);
        }
    }

    // Test if obtained results are equal (numerically) to zero
    for (int64_t i = 0; i < integrationResult.Rows(); i++) {
        const bool cond = IsZero(integrationResult(i, 0));
        if(!cond){
          std::cerr << "Divergence test failed for shape function " +
                           std::to_string(i + firstInternalPhi) + " on " +
                           TGeo::TypeName() + " element\n";
        }
        REQUIRE(IsZero(integrationResult(i, 0)));
    }
    std::cout << "Divergence test passed!\n\n";
}
