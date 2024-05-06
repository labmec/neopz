//
//  TestTopology.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/6/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/catch_approx.hpp>

#include "pzmanvector.h"
#include "pzvec_extras.h"
#include "pztrnsform.h"
#include "TPZGenGrid2D.h"
#include "tpzautopointer.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZBndCondT.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "tpzpermutation.h"
#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "tpzintpoints.h"
#include "pzintel.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "TPZExtendGridDimension.h"
#include "TPZLinearAnalysis.h"
#include "TPZRefPatternTools.h"
#include "pzshtmat.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"
#include "pzshapepiramHdiv.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMatL2Product.h"
#include "TPZGeoMeshTools.h"

enum class ESpace
{
    H1,
    HCurl,
    HDiv,
    HDivConst,
    L2
};

static const std::map<ESpace, const char *> names({{ESpace::H1, "H1"},
                                                   {ESpace::HCurl, "HCurl"},
                                                   {ESpace::HDiv, "HDiv"},
                                                   {ESpace::HDivConst, "HDivConst"},
                                                   {ESpace::L2, "L2"}});

using namespace pzshape;

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.testhdivconstant");
#endif

static int tetraedra_2[6][4] =
    {
        {1, 2, 5, 4},
        {4, 7, 3, 2},
        {0, 1, 2, 4},
        {0, 2, 3, 4},
        {4, 5, 6, 2},
        {4, 6, 7, 2}};

static bool MyDoubleComparer(REAL a, REAL b)
{
    if (IsZero(a - b))
    {
        return true;
    }
    else
    {
        return false;
    }
}

static void GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem)
{
    gmesh->NodeVec().Resize((nelem + 1) * (nelem + 1) * (nelem + 1));
    for (int64_t i = 0; i <= nelem; i++)
    {
        for (int64_t j = 0; j <= nelem; j++)
        {
            for (int64_t k = 0; k <= nelem; k++)
            {
                TPZManVector<REAL, 3> x(3);
                x[0] = k * 1. / nelem;
                x[1] = j * 1. / nelem;
                x[2] = i * 1. / nelem;
                gmesh->NodeVec()[i * (nelem + 1) * (nelem + 1) + j * (nelem + 1) + k].Initialize(x, *gmesh);
            }
        }
    }
}

static const int gfluxorder = 3;
static TPZAutoPointer<TPZGeoMesh> GenerateMesh(MElementType type, int nelem = 3, int ndiv = 0);
static TPZAutoPointer<TPZGeoMesh> TetrahedralMeshCubo(int64_t nelem, int MaterialId);

static void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);

TPZAutoPointer<TPZCompMesh> CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                                        TPZMaterial *mat, const int matId,
                                        const int k, ESpace space);

int CalcRank(const TPZFMatrix<STATE> &S, const STATE tol);

template <class TSHAPE>
static void IntegralNormal();

/** @brief Create three different shape data objects starting with the same connect orders.
 * We decrease the facet connects order by one for the second shape and increase the volume connect order by two for the third shape.
 * Then, we check if the number of facet functions for shapes 1 and 3 are equal and if the number of volume functions for shapes 1 and 3 are the same.
    @param [in] kFacet polynomial order of the facet connects*/
template <class TSHAPE>
static void CheckConnectOrders(int kFacet);

template <class TSHAPE>
void TestNewHDiv(int kFacet);

/** @brief Checks if the rank of the matrix composed by the L2 product of Hdiv shape functions is equal to the number of shape functions
    Denoting by
    - phi_i : basis functions of the refered space
    This function creates the matrix
    M = phi_i * phi_j

    And then compare if rank(M) = length(phi).
    @param [in] kFacet polynomial order of the facet connects*/
template <ESpace Space, int dim>
void CheckL2ProductRank(int kFacet);

// Tests for the 'voidflux' class.
TEST_CASE("integral_normal", "[hdivconstant_tests]")
{
    std::cout << "Initializing vector_direction check\n";
    //    IntegralNormal<pzshape::TPZShapePrism>();
    IntegralNormal<pzshape::TPZShapeCube>();
    IntegralNormal<pzshape::TPZShapeTetra>();
    IntegralNormal<pzshape::TPZShapeTriang>();
    IntegralNormal<pzshape::TPZShapeQuad>();
    //    IntegralNormal<pzshape::TPZShapePrism>();
    std::cout << "Leaving integral_normal check\n";
}

TEMPLATE_TEST_CASE("Connect order compatibility", "[hdivconstant_tests]",
                   (pzshape::TPZShapeTriang),
                   (pzshape::TPZShapeQuad),
                   (pzshape::TPZShapeTetra),
                   (pzshape::TPZShapeCube))
{
    int kFacet = GENERATE(2, 3, 4);
    SECTION("kFacet=" + std::to_string(kFacet))
    {
        CheckConnectOrders<TestType>(kFacet);
    }
}

TEMPLATE_TEST_CASE("Test New HDiv Shape", "[hdivconstant_tests]",
                   (pzshape::TPZShapeTriang),
                   (pzshape::TPZShapeQuad),
                   (pzshape::TPZShapeTetra),
                   (pzshape::TPZShapeCube))
{
    int kFacet = GENERATE(1, 2, 3, 4);
    SECTION("kFacet=" + std::to_string(kFacet))
    {
        TestNewHDiv<TestType>(kFacet);
    }
}

TEMPLATE_TEST_CASE("Linear Dependence", "[hdivconstant_tests]",
                   (typename std::integral_constant<int, 2>),
                   (typename std::integral_constant<int, 3>))
{
    // SVD requires LAPACK
#ifndef PZ_USING_LAPACK
    return;
#endif
    constexpr int dim = TestType::value;

    ESpace space = GENERATE(ESpace::HDiv,
                            ESpace::HDivConst);

    SECTION(names.at(space))
    {
        int kFacet = GENERATE(1, 2, 3);
        SECTION("kFacet=" + std::to_string(kFacet))
        {
            switch (space)
            {
            case ESpace::HDiv:
            {
                CheckL2ProductRank<ESpace::HDiv, dim>(kFacet);
                break;
            }
            case ESpace::HDivConst:
            {
                CheckL2ProductRank<ESpace::HDivConst, dim>(kFacet);
                break;
            }
            default:
                break;
            }
        }
    }
}

static TPZAutoPointer<TPZGeoMesh> GenerateMesh(MElementType eltype, int nelem, int ndiv)
{
    int dimmodel = 2;
    TPZManVector<int, 3> nx(2, nelem);
    TPZManVector<REAL, 3> x0(3, 0.), x1(3, 1.);
    x1[2] = -1.;
    TPZGenGrid2D grid(nx, x0, x1);
    if (eltype == ETriangle || eltype == EPrisma)
    {
        grid.SetElementType(MMeshType::ETriangular);
    }
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    grid.Read(gmesh.operator->());
    grid.SetBC(gmesh, 4, -1);
    grid.SetBC(gmesh, 5, -1);
    grid.SetBC(gmesh, 6, -1);
    grid.SetBC(gmesh, 7, -1);

    if (eltype == ETriangle || eltype == EPrisma || eltype == ECube || eltype == EQuadrilateral)
    {
        for (int D = 0; D < ndiv; D++)
        {
            int nels = gmesh->NElements();
            for (int elem = 0; elem < nels; elem++)
            {
                TPZVec<TPZGeoEl *> filhos;
                TPZGeoEl *gel = gmesh->ElementVec()[elem];
                gel->Divide(filhos);
            }
        }

        { // queria tanto ver a malha 2d
            std::ofstream Dummyfile("GeometricMesh2d.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh.operator->(), Dummyfile, true);
        }

        if (eltype == EPrisma || eltype == ECube)
        {
            REAL thickness = 1.; // 2.;
            TPZExtendGridDimension extend(gmesh, thickness);
            int numlayers = nelem;
            int bctop = -2;
            int bcbottom = -3; // normal negativa
            gmesh = extend.ExtendedMesh(numlayers, bcbottom, bctop);
            gmesh->SetDimension(3);
            dimmodel = 3;
        }
    }
    else if (eltype == ETetraedro)
    {
        // aqui
        dimmodel = 3;
        // gmesh = CreateOneCuboWithTetraedrons(ndiv); // AQUIDOUGLAS
        ndiv = 1;
        const int64_t NumberOfEl = ndiv;
        const int matid = 1;
        gmesh = TetrahedralMeshCubo(NumberOfEl, matid);
        gmesh->SetDimension(3);
        std::ofstream arg("gmesh.txt");
        gmesh->Print(arg);
    }
    else
    {
        // Elemento nao contemplado
        DebugStop();
    }

    int Axis;
    REAL theta, dump = 0.0;

    theta = 48.0;
    Axis = 1;
    RotateGeomesh(gmesh.operator->(), theta * dump, Axis);

    theta = -45.0;
    Axis = 2;
    RotateGeomesh(gmesh.operator->(), theta * dump, Axis);

    theta = 120.0;
    Axis = 3;
    RotateGeomesh(gmesh.operator->(), theta * dump, Axis);

    {
        //  Print Geometrical Base Mesh
        std::ofstream Dummyfile2("GeometricMesh3d.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, Dummyfile2, true);
    }

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Malha Geo FINAl \n\n";
        gmesh->Print(sout);
        // mphysics->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    return gmesh;
}

TPZAutoPointer<TPZGeoMesh> TetrahedralMeshCubo(int64_t nelem, int MaterialId)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh, nelem);

    for (int64_t i = 0; i < nelem; i++)
    {
        for (int64_t j = 0; j < nelem; j++)
        {
            for (int64_t k = 0; k < nelem; k++)
            {
                TPZManVector<int64_t, 8> nodes(8, 0);
                nodes[0] = k * (nelem + 1) * (nelem + 1) + j * (nelem + 1) + i;
                nodes[1] = k * (nelem + 1) * (nelem + 1) + j * (nelem + 1) + i + 1;
                nodes[2] = k * (nelem + 1) * (nelem + 1) + (j + 1) * (nelem + 1) + i + 1;
                nodes[3] = k * (nelem + 1) * (nelem + 1) + (j + 1) * (nelem + 1) + i;
                nodes[4] = (k + 1) * (nelem + 1) * (nelem + 1) + j * (nelem + 1) + i;
                nodes[5] = (k + 1) * (nelem + 1) * (nelem + 1) + j * (nelem + 1) + i + 1;
                nodes[6] = (k + 1) * (nelem + 1) * (nelem + 1) + (j + 1) * (nelem + 1) + i + 1;
                nodes[7] = (k + 1) * (nelem + 1) * (nelem + 1) + (j + 1) * (nelem + 1) + i;
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Tetrahedral nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                for (int el = 0; el < 6; el++)
                {
                    TPZManVector<int64_t, 4> elnodes(4);
                    int64_t index;
                    for (int il = 0; il < 4; il++)
                    {
                        elnodes[il] = nodes[tetraedra_2[el][il]];
                    }
                    gmesh->CreateGeoElement(ETetraedro, elnodes, MaterialId, index, 0);
                }
            }
        }
    }
    gmesh->BuildConnectivity();

    // Boundary Conditions
    const int numelements = gmesh->NElements();
    const int bczMinus = -3, bczplus = -2, bcids = -1;
    //    const int bczMinus = -1, bczplus = -1, bcids = -1;

    for (int el = 0; el < numelements; el++)
    {
        TPZManVector<TPZGeoNode, 4> Nodefinder(4);
        TPZManVector<REAL, 3> nodecoord(3);
        TPZGeoEl *tetra = gmesh->ElementVec()[el];
        // na face x = 0
        TPZVec<int64_t> ncoordVec(0);
        int64_t sizeOfVec = 0;
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[0], 0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec - 1] = pos;
            }
        }
        if (ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide, bcids);
        }

        ncoordVec.clear();
        sizeOfVec = 0;
        // na face x = 1
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[0], 1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec - 1] = pos;
            }
        }
        if (ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide, bcids);
        }

        ncoordVec.clear();
        sizeOfVec = 0;
        // na face y = 0
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[1], 0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec - 1] = pos;
            }
        }
        if (ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide, bcids);
        }

        ncoordVec.clear();
        sizeOfVec = 0;
        // na face y = 1
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[1], 1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec - 1] = pos;
            }
        }
        if (ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide, bcids);
        }

        ncoordVec.clear();
        sizeOfVec = 0;
        // na face z = 0
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[2], 0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec - 1] = pos;
            }
        }
        if (ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide, bczMinus);
        }

        ncoordVec.clear();
        sizeOfVec = 0;
        // na face z = 1
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[2], 1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec - 1] = pos;
            }
        }
        if (ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide, bczplus);
        }
    }

    return gmesh;
}

void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis)
{
    REAL theta = (M_PI / 180.0) * CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<REAL> RotationMatrix(3, 3, 0.0);

    switch (Axis)
    {
    case 1:
    {
        RotationMatrix(0, 0) = 1.0;
        RotationMatrix(1, 1) = +cos(theta);
        RotationMatrix(1, 2) = -sin(theta);
        RotationMatrix(2, 1) = +sin(theta);
        RotationMatrix(2, 2) = +cos(theta);
    }
    break;
    case 2:
    {
        RotationMatrix(0, 0) = +cos(theta);
        RotationMatrix(0, 2) = +sin(theta);
        RotationMatrix(1, 1) = 1.0;
        RotationMatrix(2, 0) = -sin(theta);
        RotationMatrix(2, 2) = +cos(theta);
    }
    break;
    case 3:
    {
        RotationMatrix(0, 0) = +cos(theta);
        RotationMatrix(0, 1) = -sin(theta);
        RotationMatrix(1, 0) = +sin(theta);
        RotationMatrix(1, 1) = +cos(theta);
        RotationMatrix(2, 2) = 1.0;
    }
    break;
    default:
    {
        RotationMatrix(0, 0) = +cos(theta);
        RotationMatrix(0, 1) = -sin(theta);
        RotationMatrix(1, 0) = +sin(theta);
        RotationMatrix(1, 1) = +cos(theta);
        RotationMatrix(2, 2) = 1.0;
    }
    break;
    }

    TPZVec<REAL> iCoords(3, 0.0);
    TPZVec<REAL> iCoordsRotated(3, 0.0);

    // RotationMatrix.Print("Rotation = ");

    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0, 0) * iCoords[0] + RotationMatrix(0, 1) * iCoords[1] + RotationMatrix(0, 2) * iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1, 0) * iCoords[0] + RotationMatrix(1, 1) * iCoords[1] + RotationMatrix(1, 2) * iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2, 0) * iCoords[0] + RotationMatrix(2, 1) * iCoords[1] + RotationMatrix(2, 2) * iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        gmesh->NodeVec()[inode] = GeoNode;
    }
}

#include "TPZShapeHDivConstant.h"
#include "TPZShapeHDiv.h"
#include "TPZShapeNewHDiv.h"

template <class TSHAPE>
void IntegralNormal()
{
    auto gmesh = GenerateMesh(TSHAPE::Type(), 1, 0);
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->Type() != TSHAPE::Type())
            continue;
        TPZManVector<int64_t> ids(gel->NCornerNodes());
        int nnodes = gel->NCornerNodes();
        for (int no = 0; no < nnodes; no++)
            ids[no] = gel->Node(no).Id();
        TPZShapeHDivConstant<TSHAPE> shape;
        TPZShapeData data;
        TPZManVector<int> orders(TSHAPE::NFacets + 1, 1);
        TPZManVector<int> sideorient(TSHAPE::NFacets, 1);
        shape.Initialize(ids, orders, sideorient, data);
        // comment this out because it will break hcurl shape funcs
        // why was it here?
        // data.fSideTransformationId.Resize(TSHAPE::NFacets, 0);
        data.fHDiv.fSideOrient.Resize(TSHAPE::NFacets, 1);
        int nshape = shape.NHDivShapeF(data);
        TPZManVector<REAL> pt(TSHAPE::Dimension, 0.);
        TPZFNMatrix<8, REAL> phi(TSHAPE::Dimension, nshape), divphi(nshape, 1);
        shape.Shape(pt, data, phi, divphi);
        // integrate the divergence. It should be one
        divphi *= TSHAPE::RefElVolume();
        phi.Print("phi  = ", std::cout);
        divphi.Print("divphi  = ", std::cout);
    }
}

template <class TSHAPE>
void CheckConnectOrders(int kFacet)
{
    TPZMaterialDataT<REAL> data1, data2, data3;
    TPZManVector<int64_t, 27> ids(TSHAPE::NCornerNodes, 0);
    for (int i = 0; i < TSHAPE::NCornerNodes; i++)
    {
        ids[i] = i;
    }
    TPZManVector<int, 27> sideorient(TSHAPE::NFacets, 1);
    TPZManVector<int, 27> orders1(TSHAPE::NFacets + 1, kFacet);
    TPZManVector<int, 27> orders2(orders1);
    for (int i = 0; i < TSHAPE::NFacets; i++)
    {
        orders2[i] -= 1;
    }
    TPZManVector<int, 27> orders3(orders1);
    orders3[TSHAPE::NFacets] += 2;
    TPZShapeHDivConstant<TSHAPE> hdivorig;
    hdivorig.Initialize(ids, orders1, sideorient, data1);
    hdivorig.Initialize(ids, orders2, sideorient, data2);
    hdivorig.Initialize(ids, orders3, sideorient, data3);
    int nshape1 = hdivorig.NHDivShapeF(data1);
    int nshape2 = hdivorig.NHDivShapeF(data2);
    int nshape3 = hdivorig.NHDivShapeF(data3);
    int firstshape1 = 0, firstshape2 = 0, firstshape3 = 0;
    for (int i = 0; i < TSHAPE::NFacets; i++)
    {
        firstshape1 += hdivorig.NConnectShapeF(i, data1);
        firstshape2 += hdivorig.NConnectShapeF(i, data2);
        firstshape3 += hdivorig.NConnectShapeF(i, data3);
    }
    int nvolshape1 = hdivorig.NConnectShapeF(TSHAPE::NFacets, data1);
    int nvolshape2 = hdivorig.NConnectShapeF(TSHAPE::NFacets, data2);
    int nvolshape3 = hdivorig.NConnectShapeF(TSHAPE::NFacets, data3);

    auto hdivOrder1 = data1.fHDiv.fConnectOrders;
    auto hdivOrder2 = data2.fHDiv.fConnectOrders;
    auto hdivOrder3 = data3.fHDiv.fConnectOrders;

    auto hdivNshape1 = data1.fHDiv.fNumConnectShape;
    auto hdivNshape2 = data2.fHDiv.fNumConnectShape;
    auto hdivNshape3 = data3.fHDiv.fNumConnectShape;

    CAPTURE(orders1,orders2,orders3,hdivOrder1,hdivOrder2,hdivOrder3,hdivNshape1,hdivNshape2,hdivNshape3);
    
    REQUIRE(nvolshape1 == nvolshape2);
    REQUIRE(firstshape1 + nvolshape1 == nshape1);
    REQUIRE(firstshape2 + nvolshape2 == nshape2);
    REQUIRE(firstshape3 + nvolshape3 == nshape3);
    REQUIRE(firstshape1 == firstshape3);
    
    for (int i = 0; i < TSHAPE::NFacets+1; i++)
    {
        REQUIRE(hdivOrder1[i] == orders1[i]);
        REQUIRE(hdivOrder2[i] == orders2[i]);
        REQUIRE(hdivOrder3[i] == orders3[i]);
        REQUIRE(hdivNshape1[i] == hdivorig.NConnectShapeF(i, data1));
        REQUIRE(hdivNshape2[i] == hdivorig.NConnectShapeF(i, data2));
        REQUIRE(hdivNshape3[i] == hdivorig.NConnectShapeF(i, data3));
    }

    // Evaluate the shape functions at the integration point.
    // The volume shape functions must be equal for approximation spaces 1 and 2 at every integration point.
    // The facet shape functions must be equal for approximation spaces 1 and 3 at every integration point.
    int dim = TSHAPE::Dimension;
    TPZFMatrix<REAL> phi1(dim, nshape1, 0.), phi2(dim, nshape2, 0.), phi3(dim, nshape3, 0.);
    TPZFNMatrix<60, REAL> divphi1(nshape1, 1), divphi2(nshape2, 1), divphi3(nshape3, 1);
    typename TSHAPE::IntruleType intrule(3);
    int nintpoints = intrule.NPoints();
    TPZManVector<REAL, 3> point(dim, 0.);
    for (int ip = 0; ip < nintpoints; ip++)
    {
        REAL weight;
        intrule.Point(ip, point, weight);
        hdivorig.Shape(point, data1, phi1, divphi1);
        hdivorig.Shape(point, data2, phi2, divphi2);
        hdivorig.Shape(point, data3, phi3, divphi3);

        for (int i = 0; i < nvolshape1; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                REAL val1 = phi1(j, i + firstshape1);
                REAL val2 = phi2(j, i + firstshape2);
                if (abs(val1-val2) > 1.e-10)
                {
                    std::cout << "val1 " << val1 << " val2 " << val2 << std::endl;
                }
                REQUIRE((val1 - val2) == Catch::Approx(0.));
            }
        }

        for (int i = 0; i < firstshape1; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                REAL val1 = phi1(j, i);
                REAL val3 = phi3(j, i);
                if (abs(val1-val3) > 1.e-10)
                {
                    std::cout << "val1 " << val1 << " val3 " << val3 << std::endl;
                }
                REQUIRE((val1 - val3) == Catch::Approx(0.));
            }
        }
    }
}

template <class TSHAPE>
void TestNewHDiv(int kFacet)
{
    TPZMaterialDataT<REAL> dataHDiv, dataNewHDiv, dataHDivConst;
    TPZManVector<int64_t, 27> ids(TSHAPE::NCornerNodes, 0);
    for (int i = 0; i < TSHAPE::NCornerNodes; i++)
    {
        ids[i] = i;
    }
    TPZManVector<int, 27> sideorient(TSHAPE::NFacets, 1);
    TPZManVector<int, 27> orders(TSHAPE::NFacets + 1, kFacet);
    TPZShapeHDiv<TSHAPE> hdiv_std;
    TPZShapeNewHDiv<TSHAPE> hdiv_new;
    TPZShapeHDivConstant<TSHAPE> hdiv_const;
    hdiv_std.Initialize(ids, orders, sideorient, dataHDiv);
    hdiv_new.Initialize(ids, orders, sideorient, dataNewHDiv);
    hdiv_const.Initialize(ids, orders, sideorient, dataHDivConst);
    int nshape_std = hdiv_std.NShapeF(dataHDiv);
    int nshape_new = hdiv_new.NShapeF(dataNewHDiv);
    int nshape_const = hdiv_const.NHDivShapeF(dataHDivConst);
    int nfacet_std = 0, nfacet_new = 0, nfacet_const = 0;
    for (int i = 0; i < TSHAPE::NFacets; i++)
    {
        nfacet_std += hdiv_std.NConnectShapeF(i, dataHDiv);
        nfacet_new += hdiv_new.NConnectShapeF(i, dataNewHDiv);
        nfacet_const += hdiv_const.NConnectShapeF(i, dataHDivConst);
    }
    int nvol_std = hdiv_std.NConnectShapeF(TSHAPE::NFacets, dataHDiv);
    int nvol_new = hdiv_new.NConnectShapeF(TSHAPE::NFacets, dataNewHDiv);
    int nvol_const = hdiv_const.NConnectShapeF(TSHAPE::NFacets, dataHDivConst);

    auto hdiv_std_order = dataHDiv.fHDiv.fConnectOrders;
    auto hdiv_new_order = dataNewHDiv.fHDiv.fConnectOrders;
    auto hdiv_const_order = dataHDivConst.fHDiv.fConnectOrders;

    auto hdiv_std_shape = dataHDiv.fHDiv.fNumConnectShape;
    auto hdiv_new_shape = dataNewHDiv.fHDiv.fNumConnectShape;
    auto hdiv_const_shape = dataHDivConst.fHDiv.fNumConnectShape;

    CAPTURE(orders,hdiv_std_order,hdiv_new_order, hdiv_const_order, hdiv_std_shape, hdiv_new_shape, hdiv_const_shape);
    
    REQUIRE(nvol_std == nvol_new);
    REQUIRE(nfacet_new == nfacet_const);
    REQUIRE(nfacet_std + nvol_std == nshape_std);
    REQUIRE(nfacet_new + nvol_new == nshape_new);
    
    for (int i = 0; i < TSHAPE::NFacets+1; i++)
    {
        REQUIRE(hdiv_std_order[i] == orders[i]);
        REQUIRE(hdiv_new_order[i] == orders[i]);
        REQUIRE(hdiv_std_shape[i] == hdiv_std.NConnectShapeF(i, dataHDiv));
        REQUIRE(hdiv_new_shape[i] == hdiv_new.NConnectShapeF(i, dataNewHDiv));
        REQUIRE(hdiv_const_shape[i] == hdiv_const.NConnectShapeF(i, dataHDivConst));
    }

    // Evaluate the shape functions at the integration point.
    // The volume shape functions must be equal for hdiv_std and hdiv_new at every integration point.
    // The facet shape functions must be equal for hdiv_new and hdiv_const at every integration point.
    constexpr int dim = TSHAPE::Dimension;
    TPZFMatrix<REAL> phi_std(dim, nshape_std, 0.), phi_new(dim, nshape_new, 0.), phi_const(dim, nshape_const, 0.);
    TPZFNMatrix<60, REAL> div_std(nshape_std, 1), div_new(nshape_new, 1), div_const(nshape_const, 1);
    typename TSHAPE::IntruleType intrule(3);
    int nintpoints = intrule.NPoints();
    TPZManVector<REAL, 3> point(dim, 0.);
    for (int ip = 0; ip < nintpoints; ip++)
    {
        REAL weight;
        intrule.Point(ip, point, weight);
        hdiv_std.Shape(point, dataHDiv, phi_std, div_std);
        hdiv_new.Shape(point, dataNewHDiv, phi_new, div_new);
        hdiv_const.Shape(point, dataHDivConst, phi_const, div_const);

        for (int i = 0; i < nvol_std; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                REAL val_std = phi_std(j, i + nfacet_std);
                REAL val_new = phi_new(j, i + nfacet_new);
                if (abs(val_std-val_new) > 1.e-10)
                {
                    std::cout << "val_std " << val_std << " val_new " << val_new << std::endl;
                }
                REQUIRE((val_std - val_new) == Catch::Approx(0.));
            }
        }

        for (int i = 0; i < nfacet_const; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                REAL val_const = phi_const(j, i);
                REAL val_new = phi_new(j, i);
                if (abs(val_const - val_new) > 1.e-10)
                {
                    std::cout << "val_const " << val_const << " val_new " << val_new << std::endl;
                }
                REQUIRE((val_const - val_new) == Catch::Approx(0.));
            }
        }
    }
}

template <ESpace space, int dim>
void CheckL2ProductRank(int kFacet)
{

    /******************************************
    CHECK FUNCTION DECLARATION FOR MORE DETAILS
    *******************************************/
    // select element type
    auto elType = GENERATE(MMeshType::ETriangular,
                           MMeshType::EQuadrilateral,
                           MMeshType::ETetrahedral,
                           MMeshType::EHexahedral,
                           MMeshType::EPrismatic);

    SECTION(MMeshType_Name(elType))
    {
        if (elType == MMeshType::EPrismatic && space == ESpace::HDivConst)
        {
            return; // TODOFIX
        }
        const auto elDim = MMeshType_Dimension(elType);
        if (elType == MMeshType::EPrismatic && space == ESpace::HCurl)
        {
            return; // TODOFIX
        }

        // if dimension does not correspond, skip
        if (elDim != dim)
            return;

        // creating mesh
        TPZAutoPointer<TPZGeoMesh> gmesh{nullptr};
        TPZAutoPointer<TPZCompMesh> cmesh{nullptr};

        constexpr int matId = 1;

        auto *material = new TPZMatL2Product(matId, dim);

        constexpr bool singleElement{true};

        if constexpr (singleElement)
        {
            constexpr bool createBoundEls{false};
            gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(elType, matId, createBoundEls);
        }
        else
        {
            gmesh = CreateGMesh(dim, elType, matId);
        }

        auto EnrichMesh = [](TPZAutoPointer<TPZCompMesh> cmesh, const int k)
        {
            for (auto cel : cmesh->ElementVec())
            {
                if (cel->Dimension() != dim)
                    continue;
                auto intel =
                    dynamic_cast<TPZInterpolatedElement *>(cel);
                const auto nsides = intel->Reference()->NSides();
                intel->SetSideOrder(nsides - 1, k + 1);
                intel->AdjustIntegrationRule();
            }
            cmesh->ExpandSolution();
        };

        cmesh = CreateCMesh(gmesh, material, matId, kFacet, space);

        int enrichLvl = GENERATE(-1, 0, 1);
        // SECTION("Enrichment Level: " + std::to_string(enrichLvl+1));
        EnrichMesh(cmesh, kFacet + enrichLvl);

        // for debuggin
        const auto spaceName = names.at(space);
        const auto elName = MMeshType_Name(elType);
        CAPTURE(spaceName, elName);

        // stores the matrices
        TPZAutoPointer<TPZMatrix<STATE>> matPtr = nullptr;
        {
            constexpr int nThreads{0};
            TPZLinearAnalysis an(cmesh, RenumType::ENone);

            TPZFStructMatrix<STATE> strmtrx(cmesh);
            strmtrx.SetNumThreads(nThreads);
            an.SetStructuralMatrix(strmtrx);

            TPZStepSolver<STATE> step;
            step.SetDirect(ELU);

            an.SetSolver(step);

            an.Assemble();

            matPtr = an.MatrixSolver<STATE>().Matrix();
        }

        TPZFMatrix<STATE> mat(*matPtr);

        const int nShapeHDiv = mat.Rows();

        TPZFMatrix<STATE> S;
        {
            TPZFMatrix<STATE> Udummy, VTdummy;
            mat.SVD(Udummy, S, VTdummy, 'N', 'N');
        }
        static constexpr auto tol = std::numeric_limits<STATE>::epsilon() * 10000;

        TPZManVector<int, 27> connectorder(cmesh->NConnects(), 0);
        for (int i = 0; i < cmesh->NConnects(); i++)
        {
            connectorder[i] = cmesh->ConnectVec()[i].Order();
        }

        TPZManVector<int, 27> nshape(cmesh->NConnects(), 0);
        for (int i = 0; i < cmesh->NConnects(); i++)
        {
            nshape[i] = cmesh->ConnectVec()[i].NShape();
        }

        // rank of matrix (phi_i * phi_j)
        const int rank = CalcRank(S, tol);

        CAPTURE(kFacet, nShapeHDiv, rank);
        CAPTURE(connectorder, nshape);
        REQUIRE(rank == nShapeHDiv);
    }
}

TPZAutoPointer<TPZCompMesh> CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                                        TPZMaterial *mat, const int matId,
                                        const int k, ESpace space)
{

    TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
    const int nel = cmesh->NElements();
    cmesh->SetDefaultOrder(k);
    cmesh->SetDimModel(gmesh->Dimension());
    cmesh->InsertMaterialObject(mat);
    switch (space)
    {
    case ESpace::H1:
        cmesh->SetAllCreateFunctionsContinuous();
        break;
    case ESpace::HCurl:
        cmesh->SetAllCreateFunctionsHCurl();
        break;
    case ESpace::HDiv:
        cmesh->SetAllCreateFunctionsHDiv();
        break;
    case ESpace::HDivConst:
        cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
        cmesh->SetAllCreateFunctionsHDiv();
        break;
    case ESpace::L2:
        if (k == 0)
            cmesh->SetAllCreateFunctionsDiscontinuous();
        else
            cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
        break;
    }   
    cmesh->AutoBuild();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

int CalcRank(const TPZFMatrix<STATE> &S, const STATE tol)
{
    int rank = 0;
    const int dimMat = S.Rows();
    for (int i = 0; i < dimMat; i++)
    {
        rank += S.GetVal(i, 0) > tol ? 1 : 0;
    }
    return rank;
};

template void IntegralNormal<pzshape::TPZShapeQuad>();
