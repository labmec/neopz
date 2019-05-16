
#include <pzgeoelrefless.h>
#include <tpzarc3d.h>
#include <TPZVTKGeoMesh.h>
#include <pzanalysis.h>
#include <pzgeotetrahedra.h>
#include <TPZGeoCube.h>
#include <TPZTriangleSphere.h>
#include <TPZQuadSphere.h>
#include "tpzgeoblend.h"
#include "pzgeoquad.h"
#include "pzgmesh.h"
#include "pzcmesh.h"

namespace blendtest{
    enum meshType {
        EQuad = 0, ETriang = 1, ETetra = 2, EHexa = 3
    };
    void CreateGeoMesh2D(TPZGeoMesh *&gmesh, meshType elType, int nDiv, bool printGMesh, std::string prefix);
    void CreateGeoMesh3D(TPZGeoMesh *&gmesh, meshType elType, int nDiv, bool printGMesh, std::string prefix);
}


using namespace blendtest;

/**
 * Routine for testing the TPZGeoBlend refactoring
 * @return
 */
int main()
{
    InitializePZLOG();
    bool printGMesh = true;
    bool newBlend = true;
    int nDiv = 4;

//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(ETriangle);
//    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
//    gRefDBase.InitializeUniformRefPattern(ETetraedro);
//    gRefDBase.InitializeUniformRefPattern(EPiramide);
//    gRefDBase.InitializeUniformRefPattern(EPrisma);
//    gRefDBase.InitializeUniformRefPattern(ECube);

    TPZGeoMesh * gmesh = nullptr;
    meshType elType = ETetra;

    switch(elType){
        case EQuad:
            pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>::fUseNewX = newBlend;
            break;
        case ETriang:
            pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>::fUseNewX = newBlend;
            break;
        case ETetra:
            pzgeom::TPZGeoBlend<pzgeom::TPZGeoTetrahedra>::fUseNewX = newBlend;
            break;
        case EHexa:
            pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube>::fUseNewX = newBlend;
            break;
        default:
            DebugStop();
    }
    std::string prefix = newBlend? "new" : "old";
    if(elType == EQuad || elType == ETriang){
        CreateGeoMesh2D(gmesh, elType, nDiv, printGMesh, prefix);
    }else{
        CreateGeoMesh3D(gmesh, elType, nDiv, printGMesh, prefix);
    }
    delete gmesh;

//    pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>::fUseNewX = true;
//    CreateGeoMesh2D(gmesh, printGMesh, nDiv);
//    delete gmesh;


}

namespace blendtest {
    void CreateGeoMesh2D(TPZGeoMesh *&gmesh, meshType elType, int nDiv, bool printGMesh, std::string prefix) {
        gmesh = new TPZGeoMesh();

        int64_t nNodes= 5;
        gmesh->NodeVec().Resize(nNodes);

        TPZVec<REAL> coord(3, 0.);
        for (int64_t i = 0; i < 4; i++) {
            coord[0] = ((i + 1) / 2) % 2;
            coord[1] = i / 2;
            gmesh->NodeVec()[i].SetCoord(coord);
            gmesh->NodeVec()[i].SetNodeId(i);
        }
        coord[0] = 0.5 + M_SQRT2 / 2;
        coord[1] = 0.5;

        gmesh->NodeVec()[4].SetCoord(coord);
        gmesh->NodeVec()[4].SetNodeId(4);

        int matIdQuad = 1, matIdArc = 2;
        TPZVec<int64_t> nodesIdVec(1);
        int64_t elId = 0;
        TPZGeoEl *volEl = nullptr;
        TPZGeoElRefPattern<pzgeom::TPZArc3D> *arc = nullptr;
        switch (elType) {
            case EQuad: {
                nodesIdVec.Resize(4);
                nodesIdVec[0] = 0;
                nodesIdVec[1] = 1;
                nodesIdVec[2] = 2;
                nodesIdVec[3] = 3;

                TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > *quadEl =
                        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(nodesIdVec, matIdQuad,
                                                                                         *gmesh);
                volEl = quadEl;

            }
                break;
            case ETriang: {
                nodesIdVec.Resize(3);
                nodesIdVec[0] = 0;
                nodesIdVec[1] = 1;
                nodesIdVec[2] = 2;

                TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> > *triangEl =
                        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >(nodesIdVec, matIdQuad,
                                                                                             *gmesh);
                nodesIdVec[0] = 0;
                nodesIdVec[1] = 3;
                nodesIdVec[2] = 2;
                triangEl =
                        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >(nodesIdVec, matIdQuad,
                                                                                             *gmesh);
                volEl = triangEl;

            }
                break;
            default:
                DebugStop();
        }
        nodesIdVec.Resize(3);
        {
            nodesIdVec[0] = 1;
            nodesIdVec[1] = 2;
            nodesIdVec[2] = 4;
            arc =
                    new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
        }
        gmesh->BuildConnectivity();
        TPZManVector<REAL, 3> qsi(2, 0);
        qsi[0] = 0.25;
        qsi[1] = 0.;
        TPZManVector<REAL, 3> result(3, 0);
        volEl->X(qsi, result);

        {
            TPZVec<TPZGeoEl *> sons;
            for (int iDiv = 0; iDiv < nDiv; iDiv++) {
                int nel = gmesh->NElements();
                for (int iel = 0; iel < nel; iel++) {
                    TPZGeoEl *geo = gmesh->Element(iel);
                    TPZGeoEl *current = geo;
                    int nFather = 0;
                    while (current->Father()) {
                        current = current->Father();
                        nFather++;
                    }
                    if (nFather <= iDiv && !geo->HasSubElement()) {
                        geo->Divide(sons);
                        nel = gmesh->NElements();
                    }
                }
            }
        }

        if (printGMesh) {
            std::string meshFileName = prefix + "gmesh";
            const size_t strlen = meshFileName.length();
            meshFileName.append(".vtk");
            std::ofstream outVTK(meshFileName.c_str());
            meshFileName.replace(strlen, 4, ".txt");
            std::ofstream outTXT(meshFileName.c_str());

            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
            gmesh->Print(outTXT);
            outTXT.close();
            outVTK.close();
        }
    }

    void CreateGeoMesh3D(TPZGeoMesh *&gmesh, meshType elType, int nDiv, bool printGMesh, std::string prefix){
        gmesh = new TPZGeoMesh();

        const int64_t nNodesHexa = 8;
        const int64_t nNodes= nNodesHexa;//+nNodesArc;
        gmesh->NodeVec().Resize(nNodes);

        const REAL radius = std::sqrt(3);
        const REAL cubeSide= 2;
        TPZVec<REAL> coord(3, 0.);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///                                                                                                         ////
        ////the face containing the nodes 0,1,2 and 3 will be joined with a half-sphere                             ////
        ///                                                                                                         ////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for (int64_t i = 0; i < nNodesHexa; i++) {
            int iAux = i%4;
            bool xFactor = (bool) !( (iAux % 2) == (iAux / 2) );
            bool yFactor = (bool) ( i / 4 );
            bool zFactor = (bool) ( ( iAux /2 ) % 2 );
            coord[0] = 0.5 * cubeSide * (2 * (int)xFactor - 1) ;
            coord[1] = 0.5 * cubeSide  * (1 - 2 * (int)yFactor);
            coord[2] = 0.5 * cubeSide * (2 * (int)zFactor - 1) ;
            gmesh->NodeVec()[i].SetCoord(coord);
            gmesh->NodeVec()[i].SetNodeId(i);
//            std::cout<<"x : "<<xFactor<<"\ty : "<<yFactor<<"\tz : "<<zFactor<<std::endl;
//            std::cout<<"x : "<<coord[0]<<"\ty : "<<coord[1]<<"\tz : "<<coord[2]<<std::endl;
        }


        int matIdVol = 1, matIdSphere = 2;
        TPZVec<int64_t> nodesIdVec(1);
        int64_t elId = 0;
        TPZGeoEl *volEl = nullptr;
        TPZGeoElRefPattern<pzgeom::TPZQuadSphere<>> *sphere = nullptr;

        nodesIdVec.Resize(8);
        TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > *hexaEl = nullptr;

        nodesIdVec[0] = 0;
        nodesIdVec[1] = 1;
        nodesIdVec[2] = 2;
        nodesIdVec[3] = 3;
        nodesIdVec[4] = 4;
        nodesIdVec[5] = 5;
        nodesIdVec[6] = 6;
        nodesIdVec[7] = 7;

        hexaEl =
                new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >(nodesIdVec, matIdVol,
                                                                                 *gmesh);
        volEl = hexaEl;
        TPZAutoPointer<TPZRefPattern> refp;
        switch (elType){
            case ETetra:{
                char buf[] =
                        "12 5 "
                        "-99 CubeIntoAll "
                        "-1 -1  -1 "
                        " 1 -1  -1 "
                        " 1  1  -1 "
                        "-1  1  -1 "
                        "-1 -1   0 "
                        " 1 -1   0 "
                        " 1  1   0 "
                        "-1  1   0 "
                        "-1 -1   1 "
                        " 1 -1   1 "
                        " 1  1   1 "
                        "-1  1   1 "
                        "7 8  0 1 2 3 8 9 10 11 " //master cube
                        "7 8  4 5 6 7 8 9 10 11 " //upper half of the cube
                        "6 6  0 5 4 3 6 7 "//prism
                        "5 5  0 3 6 5 2 " //pyramid
                        "4 4  1 2 5 0 "; //tetra

                std::istringstream str(buf);
                refp = new TPZRefPattern(str);
                refp->GenerateSideRefPatterns();
                gRefDBase.InsertRefPattern(refp);
                const int firstFace = pztopology::TPZCube::NumSides(0) + pztopology::TPZCube::NumSides(1);
                const int nFaces = pztopology::TPZCube::NumSides(2);
                hexaEl->SetRefPattern(refp);
                for(int face = firstFace; face < firstFace + nFaces; face++){
                    const int nNodes = pztopology::TPZCube::NSideNodes(face);
                    nodesIdVec.Resize(nNodes);
                    for (int node = 0; node < nNodes; node++){
                       nodesIdVec[node] = pztopology::TPZCube::SideNodeLocId(face,node);
                    }
                    sphere =
                            new TPZGeoElRefPattern<pzgeom::TPZQuadSphere<>>(nodesIdVec, matIdSphere, *gmesh);
                    coord[0] = coord[1] = coord[2] = 0;
                    sphere->Geom().SetData(radius, coord);
                    auto sphereRefp = refp->SideRefPattern(face);
                    sphere->SetRefPattern(sphereRefp);
                }
                gmesh->BuildConnectivity();

                TPZVec<TPZGeoEl *> sons;
                const int nel = gmesh->NElements();
                TPZGeoEl *hexaEl = gmesh->Element(0);
                hexaEl->Divide(sons);
                for (int iel = 1; iel < nel; iel++) {
                    TPZGeoEl *geo = gmesh->Element(iel);
                    auto geoElSide = geo->Neighbour(geo->NSides()-1);
                    if(geoElSide.NSubElements() > 1)
                    {
                        geo->Divide(sons);
                    }
                }
                gmesh->ResetConnectivities();
                gmesh->BuildConnectivity();
            }
                break;
            case EHexa://Do nothing
                break;
        }
        if (printGMesh) {
            std::string meshFileName = prefix + "gmesh_partial";
            const size_t strlen = meshFileName.length();
            meshFileName.append(".vtk");
            std::ofstream outVTK(meshFileName.c_str());
            meshFileName.replace(strlen, 4, ".txt");
            std::ofstream outTXT(meshFileName.c_str());

            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
            gmesh->Print(outTXT);
            outTXT.close();
            outVTK.close();
        }
//        if(volEl){
//            TPZManVector<REAL, 3> qsi(3, 0);
//            qsi[0] = 0.25;
//            qsi[1] = 0.25;
//            qsi[2] = 0.25;
//            TPZManVector<REAL, 3> result(3, 0);
//            volEl->X(qsi, result);
//            qsi[0] = 0.3333333333333333332;
//            qsi[1] = 0.3333333333333333332;
//            qsi[2] = 0.3333333333333333332;
//            volEl->X(qsi, result);
//        }

        {
            TPZVec<TPZGeoEl *> sons;
            for (int iDiv = 0; iDiv < nDiv; iDiv++) {
                int nel = gmesh->NElements();
                for (int iel = 0; iel < nel; iel++) {
                    TPZGeoEl *geo = gmesh->Element(iel);
                    if (!geo->HasSubElement()) {
                        geo->Divide(sons);
                    }
                }
            }
        }

        if (printGMesh) {
            std::string meshFileName = prefix + "gmesh";
            const size_t strlen = meshFileName.length();
            meshFileName.append(".vtk");
            std::ofstream outVTK(meshFileName.c_str());
            meshFileName.replace(strlen, 4, ".txt");
            std::ofstream outTXT(meshFileName.c_str());

            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
            gmesh->Print(outTXT);
            outTXT.close();
            outVTK.close();
        }
    }
}