#include "TPZSymTensorApproxCreator.h"
#include "Elasticity/TPZMixedSymElasticityND.h"
#include "TPZGenGrid2D.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZLinearAnalysis.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeoelbc.h"
#include "pzcondensedcompel.h"
#include <catch2/catch_test_macros.hpp>
#include <memory>
#include <cmath>
#include <exception>

namespace {
std::unique_ptr<TPZGeoMesh> Geometry(bool triangles, int domain = 1, bool neumannRight = false) {
    auto mesh = std::make_unique<TPZGeoMesh>();
    TPZGenGrid2D grid({2,1},{0.,0.,0.},{2.,1.,0.});
    grid.SetElementType(triangles ? MMeshType::ETriangular : MMeshType::EQuadrilateral);
    grid.SetRefpatternElements(true);
    grid.Read(mesh.get(),domain);
    for(int side = 4; side < 8; side++)
        grid.SetBC(mesh.get(),side,neumannRight && side == 5 ? -2 : -1);
    return mesh;
}
void AddMaterial(TPZSymTensorApproxCreator &creator, bool neumannRight = false) {
    auto mat = new TPZMixedSymElasticityND(creator.GetMaterialIds().fDomain,2);
    mat->SetElasticity(1.,0.3);
    mat->SetPlaneStress();
    creator.InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<STATE> val2(2,0.);
    auto bc = mat->CreateBC(mat,-1,0,val1,val2);
    bc->SetForcingFunctionBC([](const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &du) {
        u.Resize(2); u[0] = x[0]; u[1] = x[1];
        du.Redim(2,2); du(0,0) = 1.; du(1,1) = 1.;
    },1);
    creator.InsertMaterialObject(bc);
    if(neumannRight) {
        auto neumann = mat->CreateBC(mat,-2,1,val1,val2);
        creator.InsertMaterialObject(neumann);
    }
}
void DeleteMeshes(TPZMultiphysicsCompMesh *mesh) {
    auto atomic = mesh->MeshVector();
    delete mesh;
    for(auto part : atomic) delete part;
}
}

TEST_CASE("Symmetric tensor default order updates hybridization through the base class", "[symtensor][pre][order]") {
    TPZSymTensorApproxCreator creator(nullptr);
    REQUIRE(creator.GetLowOrder() == 1);
    REQUIRE(creator.GetDefaultOrder() == 2);
    REQUIRE(creator.HybridType() == HybridizationType::EStandardSquared);
    TPZApproxCreator &base = creator;
    base.SetDefaultOrder(2);
    REQUIRE(creator.GetDefaultOrder() == 2);
    REQUIRE(creator.HybridType() == HybridizationType::EStandardSquared);

    creator.SetLowOrder(1);
    REQUIRE(creator.HybridType() == HybridizationType::EStandardSquared);
    base.SetDefaultOrder(1);
    REQUIRE(creator.GetDefaultOrder() == 1);
    REQUIRE(creator.HybridType() == HybridizationType::EStandardSquared);
    base.SetDefaultOrder(3);
    REQUIRE(creator.HybridType() == HybridizationType::EStandardSquared);

    creator.SetHybridType(HybridizationType::ESemiHybrid);
    base.SetDefaultOrder(1);
    REQUIRE(creator.HybridType() == HybridizationType::ESemiHybrid);
    base.SetDefaultOrder(2);
    REQUIRE(creator.HybridType() == HybridizationType::ESemiHybrid);
}

TEST_CASE("Symmetric tensor creator builds and solves macroelements", "[symtensor][pre]") {
    for(bool triangles : {false,true}) {
        TPZManVector<int> orders = {1,2};
        if(triangles == false) orders = {2,3};
        for(int order : orders) {
            for(int low : {1,2}) {
                if(low > order) continue;
                for(bool condense : {true,false}) {
                    CAPTURE(triangles,order,low,condense);
                    auto geo = Geometry(triangles);
                    TPZSymTensorApproxCreator creator(geo.get());
                    creator.SetDefaultOrder(order);
                    creator.SetLowOrder(low);
                    creator.SetShouldCondense(condense);
                    AddMaterial(creator);
                    auto mesh = creator.CreateApproximationSpace();
                    REQUIRE(mesh->MeshVector().size() == 2);
                    REQUIRE(mesh->MeshVector()[0]->Dimension() == 2);
                    REQUIRE(mesh->MeshVector()[1]->Dimension() == 2);
                    REQUIRE(mesh->NEquations() > 0);
                    REQUIRE_THROWS_AS(creator.CreateApproximationSpace(),std::logic_error);
                    int nvol = 0, nwrap = 0, nlow = 0, ngroups = 0;
                    for(auto gel : geo->ElementVec()) {
                        if(!gel || gel->HasSubElement()) continue;
                        nvol += gel->MaterialId() == 1;
                        nwrap += gel->MaterialId() == 5;
                        nlow += gel->MaterialId() == 11;
                    }
                    REQUIRE(nvol == (triangles ? 12 : 8));
                    REQUIRE(nwrap == 3*nvol);
                    // Dirichlet boundaries end at the low-order traction;
                    // only internal macro boundaries need low-order displacement.
                    REQUIRE(nlow == (triangles ? 3 : 1));
                    for(auto cel : mesh->ElementVec()) {
                        if(dynamic_cast<TPZCondensedCompElT<STATE> *>(cel)) ngroups++;
                    }
                    REQUIRE(ngroups == (condense ? (triangles ? 4 : 2) : 0));
                    // Constant isotropic strain is representable even with constant traces:
                    // all boundary averages are exact, and the stress is constant.
                    TPZLinearAnalysis analysis(mesh);
                    TPZFStructMatrix<STATE> structure(mesh);
                    analysis.SetStructuralMatrix(structure);
                    TPZStepSolver<STATE> solver; solver.SetDirect(ELU);
                    analysis.SetSolver(solver);
                    analysis.Run();
                    REQUIRE(std::isfinite(Norm(mesh->Solution())));
                    // Check the computed stress at each atomic triangle center.
                    mesh->LoadReferences();
                    for(auto cel : mesh->MeshVector()[0]->ElementVec()) {
                        if(!cel || cel->Reference()->Dimension() != 2) continue;
                        TPZManVector<REAL,3> qsi(2,0.);
                        cel->Reference()->CenterPoint(cel->Reference()->NSides()-1,qsi);
                        TPZManVector<STATE,3> stress(3,0.);
                        cel->Solution(qsi,0,stress);
                        REQUIRE(stress.size() == 3);
                        if(stress[1]-1./0.7 >= 1.e-5) {
                            std::cout << "I should stop\n";
                        }
                        REQUIRE(std::abs(stress[0]-1./0.7) < 1.e-4);
                        REQUIRE(std::abs(stress[1]-1./0.7) < 1.e-4);
                        REQUIRE(std::abs(stress[2]) < 1.e-4);
                    }
                    DeleteMeshes(mesh);
                }
            }
        }
    }
}

TEST_CASE("Symmetric tensor creator distinguishes Dirichlet and Neumann hybridization", "[symtensor][pre]") {
    for(bool triangles : {false,true}) {
        CAPTURE(triangles);
        auto geo = Geometry(triangles,1,true);
        TPZSymTensorApproxCreator creator(geo.get());
        creator.SetShouldCondense(false);
        AddMaterial(creator,true);
        auto mesh = creator.CreateApproximationSpace();
        const auto &ids = creator.GetMaterialIds();
        int nlow = 0, ndirichlet = 0, nneumann = 0;
        for(auto gel : geo->ElementVec()) {
            if(!gel || gel->HasSubElement()) continue;
            nlow += gel->MaterialId() == ids.fDisplacementLow;
            if(gel->MaterialId() != -1 && gel->MaterialId() != -2) continue;
            TPZGeoElSide side(gel);
            REQUIRE(side.HasNeighbour(ids.fTractionLow));
            const bool hasLowDisplacement = bool(side.HasNeighbour(ids.fDisplacementLow));
            if(gel->MaterialId() == -1) {
                ndirichlet++;
                REQUIRE_FALSE(hasLowDisplacement);
            } else {
                nneumann++;
                REQUIRE(hasLowDisplacement);
            }
        }
        REQUIRE(ndirichlet == 5);
        REQUIRE(nneumann == 1);
        REQUIRE(nlow == (triangles ? 3 : 1) + nneumann);
        DeleteMeshes(mesh);
    }
}

TEST_CASE("Symmetric tensor creator validates parameters and material ids", "[symtensor][pre]") {
    auto geo = Geometry(false,21);
    TPZSymTensorApproxCreator creator(geo.get());
    REQUIRE(creator.HybridData().fWrapMatId == creator.GetMaterialIds().fWrap);
    REQUIRE(creator.HybridData().fLagrangeData.size() == 3);
    REQUIRE(creator.HybridData().fLagrangeData[0].fLagrangeMatId == creator.GetMaterialIds().fDisplacement);
    REQUIRE(creator.HybridData().fLagrangeData[1].fLagrangeMatId == creator.GetMaterialIds().fTractionLow);
    REQUIRE(creator.HybridData().fLagrangeData[2].fLagrangeMatId == creator.GetMaterialIds().fDisplacementLow);
    auto ids = creator.GetMaterialIds(); ids.fDomain = 21;
    creator.SetMaterialIds(ids);
    REQUIRE_THROWS_AS(creator.CheckSetupConsistency(),std::invalid_argument);
    AddMaterial(creator);
    creator.SetDefaultOrder(0);
    REQUIRE_THROWS_AS(creator.CheckSetupConsistency(),std::invalid_argument);
    creator.SetDefaultOrder(2); creator.SetLowOrder(3);
    REQUIRE_THROWS_AS(creator.CheckSetupConsistency(),std::invalid_argument);
    REQUIRE_THROWS_AS(creator.SetLowOrder(-1),std::bad_exception);
    creator.SetLowOrder(1);
    REQUIRE_THROWS_AS(creator.SetLowOrder(0),std::invalid_argument);
    REQUIRE(creator.GetLowOrder() == 1);
    REQUIRE_NOTHROW(creator.CheckSetupConsistency());
    REQUIRE_THROWS_AS(creator.SetProbType(ProblemType::EDarcy),std::invalid_argument);
    REQUIRE_THROWS_AS(creator.SetHybridType(HybridizationType::EStandard),std::invalid_argument);
    ids.fWrap = -1; creator.SetMaterialIds(ids);
    REQUIRE_THROWS_AS(creator.CheckSetupConsistency(),std::invalid_argument);
    ids.fWrap = 105;
    ids.fMacroBoundary = 104;
    ids.fInterfacePositive = 106;
    ids.fInterfaceNegative = 107;
    ids.fDisplacement = 108;
    ids.fInterfacePositiveLow = 109;
    ids.fInterfaceNegativeLow = 110;
    ids.fDisplacementLow = 111;
    ids.fTractionLow = 112;
    creator.SetMaterialIds(ids);
    REQUIRE(creator.HybridData().fWrapMatId == 105);
    REQUIRE(creator.HybridData().fLagrangeData.size() == 3);
    REQUIRE(creator.HybridData().fLagrangeData[0].fLagrangeMatId == ids.fDisplacement);
    REQUIRE(creator.HybridData().fLagrangeData[1].fLagrangeMatId == ids.fTractionLow);
    REQUIRE(creator.HybridData().fLagrangeData[2].fLagrangeMatId == ids.fDisplacementLow);
    creator.SetLowOrder(1);
    creator.SetShouldCondense(false);
    auto mesh = creator.CreateApproximationSpace();
    REQUIRE(mesh->FindMaterial(21) != nullptr);
    REQUIRE(mesh->FindMaterial(creator.HybridData().fWrapMatId) != nullptr);
    DeleteMeshes(mesh);
}

TEST_CASE("Symmetric tensor creator accepts a prepared polygon fan", "[symtensor][pre]") {
    for(int low : {1,2}) {
        TPZGeoMesh geo;
        geo.SetDimension(2);
        constexpr int faces = 6;
        geo.NodeVec().Resize(faces+1);
        TPZManVector<REAL,3> center(3,0.);
        geo.NodeVec()[0].Initialize(center,geo);
        for(int i = 0; i < faces; i++) {
            const REAL angle = 2.*std::acos(-1.)*i/faces;
            TPZManVector<REAL,3> point = {std::cos(angle),std::sin(angle),0.};
            geo.NodeVec()[i+1].Initialize(point,geo);
        }
        for(int i = 0; i < faces; i++) {
            TPZManVector<int64_t,3> nodes = {i+1,(i+1)%faces+1,0};
            int64_t index;
            geo.CreateGeoElement(ETriangle,nodes,1,index);
        }
        geo.BuildConnectivity();
        for(int i = 0; i < faces; i++) {
            TPZGeoElBC(geo.Element(i),3,4);
            TPZGeoElBC(geo.Element(i),3,-1);
        }
        TPZSymTensorApproxCreator creator(&geo);
        creator.SetGeometryIsPrepared(true);
        creator.SetDefaultOrder(2);
        creator.SetLowOrder(low);
        AddMaterial(creator);
        auto mesh = creator.CreateApproximationSpace();
        int ngroups = 0;
        for(auto cel : mesh->ElementVec())
            if(dynamic_cast<TPZCondensedCompElT<STATE> *>(cel)) ngroups++;
        REQUIRE(ngroups == 1);
        // A single condensed macroelement with Dirichlet data on every
        // boundary has no global displacement unknowns.
        REQUIRE(mesh->NEquations() == 0);
        for(auto gel : geo.ElementVec()) {
            if(!gel || gel->HasSubElement()) continue;
            REQUIRE(gel->MaterialId() != creator.GetMaterialIds().fDisplacementLow);
        }
        DeleteMeshes(mesh);
    }
}
