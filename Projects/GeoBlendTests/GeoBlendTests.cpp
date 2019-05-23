
#include <pzgeoelrefless.h>
#include <tpzarc3d.h>
#include <TPZVTKGeoMesh.h>
#include <pzanalysis.h>
#include <pzgeotetrahedra.h>
#include <TPZGeoCube.h>
#include <pzgeopyramid.h>
#include <TPZTriangleSphere.h>
#include <TPZQuadSphere.h>
#include "tpzgeoblend.h"
#include "pzgeoquad.h"
#include "pzgmesh.h"
#include "pzcmesh.h"

namespace blendtest{
    void CreateGeoMesh2D(TPZGeoMesh *&gmesh, int nDiv, bool printGMesh, std::string prefix);
    void CreateGeoMesh3D(TPZGeoMesh *&gmesh, int nDiv, bool printGMesh, std::string prefix);
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
    int nDiv = 3;
    bool run3d = true;

//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(ETriangle);
//    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
//    gRefDBase.InitializeUniformRefPattern(ETetraedro);
//    gRefDBase.InitializeUniformRefPattern(EPiramide);
//    gRefDBase.InitializeUniformRefPattern(EPrisma);
//    gRefDBase.InitializeUniformRefPattern(ECube);


    pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>::fUseNewX = newBlend;
    pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>::fUseNewX = newBlend;
    pzgeom::TPZGeoBlend<pzgeom::TPZGeoTetrahedra>::fUseNewX = newBlend;
    pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube>::fUseNewX = newBlend;
    pzgeom::TPZGeoBlend<pzgeom::TPZGeoPrism>::fUseNewX = newBlend;
    pzgeom::TPZGeoBlend<pzgeom::TPZGeoPyramid>::fUseNewX = newBlend;

    std::string prefix = newBlend? "new" : "old";
    if(run3d){
        TPZGeoMesh * gmesh = nullptr;
        CreateGeoMesh3D(gmesh, nDiv, printGMesh, prefix);
        delete gmesh;
    }
    else{
        TPZGeoMesh * gmesh = nullptr;
        CreateGeoMesh2D(gmesh, nDiv, printGMesh, prefix);
        delete gmesh;
    }

//    pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>::fUseNewX = true;
//    CreateGeoMesh2D(gmesh, printGMesh, nDiv);
//    delete gmesh;


}

namespace blendtest {
    void CreateGeoMesh2D(TPZGeoMesh *&gmesh, int nDiv, bool printGMesh, std::string prefix) {
        gmesh = new TPZGeoMesh();

        int64_t nNodes= 8;
        gmesh->NodeVec().Resize(nNodes);

        TPZVec<REAL> coord(3, 0.);
        for (int64_t i = 0; i < 4; i++) {
            coord[0] = i/2 == i % 2 ? -1 : 1;
            coord[1] = -1 + 2 * (i/2);
//            std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
            gmesh->NodeVec()[i].SetCoord(coord);
            gmesh->NodeVec()[i].SetNodeId(i);
        }//quad nodes
        for (int64_t i = 0; i < 4; i++) {
            const int xOff = i%2;
            const int yOff = 1 - xOff;
            const int xSign = 1- 2*(i/2);
            const int ySign = -1 * xSign;
            coord[0] = xSign * xOff*M_SQRT2;
            coord[1] = ySign * yOff*M_SQRT2;
//            std::cout<<"xOff:\t"<<xOff<<"\txSign:"<<xSign<<std::endl;
//            std::cout<<"yOff:\t"<<yOff<<"\tySign:"<<ySign<<std::endl;
//            std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
            gmesh->NodeVec()[4+i].SetCoord(coord);
            gmesh->NodeVec()[4+i].SetNodeId(4+i);
        }//quad nodes

        int matIdQuad = 1, matIdArc = 2;
        TPZVec<int64_t> nodesIdVec(1);
        int64_t elId = 0;

        nodesIdVec.Resize(4);
        TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > *quadEl = nullptr;
        {
            nodesIdVec[0] = 0;
            nodesIdVec[1] = 1;
            nodesIdVec[2] = 2;
            nodesIdVec[3] = 3;
            quadEl = new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(nodesIdVec, matIdQuad,
                                                                                      *gmesh);
        }


        TPZAutoPointer<TPZRefPattern> refp;
        char buf[] =
                "6 4 "
                "-98 QuadIntoAll "
                "-1 -1   0 "//0
                " 1 -1   0 "//1
                " 1  1   0 "//2
                "-1  1   0 "//3
                " 0 -1   0 "//4
                " 0  1   0 "//5
                "3 4  0 1 2 3 " //master quad
                "3 4  0 4 5 3 " //left half of the cube
                "2 3  4 1 2 "//triang1
                "2 3  4 2 5 "; //triang2

        std::istringstream str(buf);
        refp = new TPZRefPattern(str);
        refp->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp);
        const int firstEdge = pztopology::TPZQuadrilateral::NumSides(0);
        const int nEdges = pztopology::TPZQuadrilateral::NumSides(1);
        quadEl->SetRefPattern(refp);
        TPZGeoElRefPattern<pzgeom::TPZArc3D> *arc = nullptr;
        for(int edge = firstEdge; edge < firstEdge+ nEdges; edge++){
            const int nNodes = pztopology::TPZQuadrilateral::NSideNodes(edge);
            nodesIdVec.Resize(nNodes+1);
            for (int node = 0; node < nNodes; node++){
                nodesIdVec[node] = pztopology::TPZQuadrilateral::SideNodeLocId(edge,node);
                std::cout<<nodesIdVec[node]<<"\t";
            }
            nodesIdVec[nNodes] = edge;
            std::cout<<nodesIdVec[nNodes]<<std::endl;
            arc = new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
            auto arcRefp = refp->SideRefPattern(edge);
            arc->SetRefPattern(arcRefp);
        }
        gmesh->BuildConnectivity();

        TPZVec<TPZGeoEl *> sons;
        const int nel = gmesh->NElements();
        quadEl->Divide(sons);
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

    void CreateGeoMesh3D(TPZGeoMesh *&gmesh, int nDiv, bool printGMesh, std::string prefix){
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

        nodesIdVec.Resize(8);

        nodesIdVec[0] = 0;
        nodesIdVec[1] = 1;
        nodesIdVec[2] = 2;
        nodesIdVec[3] = 3;
        nodesIdVec[4] = 4;
        nodesIdVec[5] = 5;
        nodesIdVec[6] = 6;
        nodesIdVec[7] = 7;
        TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > *hexaEl =
                new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >(nodesIdVec, matIdVol, *gmesh);
        TPZAutoPointer<TPZRefPattern> refp;
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
        TPZGeoElRefPattern<pzgeom::TPZQuadSphere<>> *sphere = nullptr;
        for(int face = firstFace; face < firstFace + nFaces; face++){
            const int nNodes = pztopology::TPZCube::NSideNodes(face);
            nodesIdVec.Resize(nNodes);
            for (int node = 0; node < nNodes; node++){
               nodesIdVec[node] = pztopology::TPZCube::SideNodeLocId(face,node);
            }
            sphere = new TPZGeoElRefPattern<pzgeom::TPZQuadSphere<>>(nodesIdVec, matIdSphere, *gmesh);
            coord[0] = coord[1] = coord[2] = 0;
            sphere->Geom().SetData(radius, coord);
            auto sphereRefp = refp->SideRefPattern(face);
            sphere->SetRefPattern(sphereRefp);
        }
        gmesh->BuildConnectivity();

        {
            TPZVec<TPZGeoEl *> sons;
            const int nel = gmesh->NElements();
            hexaEl->Divide(sons);
            for (int iel = 1; iel < nel; iel++) {
                TPZGeoEl *geo = gmesh->Element(iel);
                auto geoElSide = geo->Neighbour(geo->NSides()-1);
                if(geoElSide.NSubElements() > 1)
                {
                    geo->Divide(sons);
                }
            }
        }

        //Create GeoBlend elements
        {
            TPZGeoMesh *newgmesh = new TPZGeoMesh();
            newgmesh->NodeVec().Resize(gmesh->NNodes());
            for (int64_t i = 0; i < gmesh->NNodes(); i++) {
                coord[0] = gmesh->NodeVec()[i].Coord(0);
                coord[1] = gmesh->NodeVec()[i].Coord(1);
                coord[2] = gmesh->NodeVec()[i].Coord(2);
                newgmesh->NodeVec()[i].SetCoord(coord);
                newgmesh->NodeVec()[i].SetNodeId(i);
                newgmesh->SetMaxNodeId(i);
            }
            //TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >
            const int nel = gmesh->NElements();
            newgmesh->ElementVec().Resize(nel);
            for (int iel = 0; iel < nel; iel++) {
                newgmesh->ElementVec()[iel] = nullptr;
                TPZGeoEl *geo = gmesh->Element(iel);
                nodesIdVec.resize(geo->NNodes());
                for(int i = 0; i < nodesIdVec.size(); i++){
                    nodesIdVec[i] = geo->NodeIndex(i);
                }
                const int matId = geo->MaterialId();
                const int elId = geo->Index();
//                geo->Ind
                switch(geo->Type()){
                    case ETetraedro:
                        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTetrahedra>>(elId,nodesIdVec,matId,*newgmesh);
                        break;
                    case EPiramide:
                        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoPyramid>>(elId,nodesIdVec,matId,*newgmesh);
                        break;
                    case EPrisma:
                        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoPrism>>(elId,nodesIdVec,matId,*newgmesh);
                        break;
                    case ECube:
                        if(geo->HasSubElement()) continue;//the original cube should not be copied
                        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube>>(elId,nodesIdVec,matId,*newgmesh);
                        break;
                    default:
                        geo->Clone(*newgmesh);
                        break; //the element should not be deleted
                }
                //TPZGeoElRefPattern(int64_t id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh);
            }
            delete gmesh;
            gmesh = newgmesh;
        }
        gmesh->ResetConnectivities();
        gmesh->BuildConnectivity();
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

        {
            TPZVec<TPZGeoEl *> sons;
            for (int iDiv = 0; iDiv < nDiv; iDiv++) {
                const int nel = gmesh->NElements();
                for (int iel = 0; iel < nel; iel++) {
                    TPZGeoEl *geo = gmesh->ElementVec()[iel];
                    if (geo && !geo->HasSubElement()) {
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