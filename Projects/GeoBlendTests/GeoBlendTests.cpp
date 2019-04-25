
#include <pzgeoelrefless.h>
#include <tpzarc3d.h>
#include <TPZVTKGeoMesh.h>
#include <pzanalysis.h>
#include <pzgeotetrahedra.h>
#include <TPZTriangleSphere.h>
#include "tpzgeoblend.h"
#include "pzgeoquad.h"
#include "pzgmesh.h"
#include "pzcmesh.h"

namespace blendtest{
    enum meshType {
        EQuad = 0, ETriang = 1, ETetra = 2
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

        const int64_t nNodesTetra = 4;
//        const int64_t nNodesArc = 3;
        const int64_t nNodes= nNodesTetra;//+nNodesArc;
        gmesh->NodeVec().Resize(nNodes);

        const REAL radius = 0.5;
        TPZVec<REAL> coord(3, 0.);

        for (int64_t i = 0; i < nNodesTetra; i++) {
            coord[0] = ((int)(i==1)) * radius;
            coord[1] = ((int)(i==2)) * radius;
            coord[2] = ((int)(i==3)) * radius;
            gmesh->NodeVec()[i].SetCoord(coord);
            gmesh->NodeVec()[i].SetNodeId(i);
        }

//        for (int64_t i = 0; i < nNodesArc; i++) {
//            coord[0] = ((i+1)%2) * radius * M_SQRT1_2;
//            coord[1] = (1-i/2) * radius * M_SQRT1_2;
//            coord[2] = ((int)(!!i))* radius * M_SQRT1_2;
//            gmesh->NodeVec()[i+nNodesTetra].SetCoord(coord);
//            gmesh->NodeVec()[i+nNodesTetra].SetNodeId(i+nNodesTetra);
//            std::cout<<"x:  "<<coord[0]<<"\ty:  "<<coord[1]<<"\tz:  "<<coord[2]<<std::endl;
//        }

        int matIdVol = 1, matIdSphere = 2;
        TPZVec<int64_t> nodesIdVec(1);
        int64_t elId = 0;
        TPZGeoEl *volEl = nullptr;
        TPZGeoElRefPattern<pzgeom::TPZTriangleSphere<>> *sphere = nullptr;
        switch (elType) {
            case ETetra: {
                nodesIdVec.Resize(4);
                nodesIdVec[0] = 0;
                nodesIdVec[1] = 1;
                nodesIdVec[2] = 2;
                nodesIdVec[3] = 3;

                TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTetrahedra> > *tetraEl =
                        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTetrahedra> >(nodesIdVec, matIdVol,
                                                                                         *gmesh);
                volEl = tetraEl;

            }
                break;
            default:
                DebugStop();
        }
        nodesIdVec.Resize(3);
        {
            nodesIdVec[0] = 1;
            nodesIdVec[1] = 2;
            nodesIdVec[2] = 3;
            sphere =
                    new TPZGeoElRefPattern<pzgeom::TPZTriangleSphere<>>(nodesIdVec, matIdSphere, *gmesh);
            coord[0]=coord[1]=coord[2] = 0;
            sphere->Geom().SetData(radius,coord);

        }
        gmesh->BuildConnectivity();
        TPZManVector<REAL, 3> qsi(3, 0);
        qsi[0] = 0.25;
        qsi[1] = 0.25;
        qsi[2] = 0.25;
        TPZManVector<REAL, 3> result(3, 0);
        volEl->X(qsi, result);
        qsi[0] = 0.3333333333333333332;
        qsi[1] = 0.3333333333333333332;
        qsi[2] = 0.3333333333333333332;
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
}