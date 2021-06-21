//
// Created by Francisco Teixeira Orlandini on 19/12/19.
//

#include "TPZGenGrid3D.h"

#include "pzgmesh.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGenGrid2D.h"//MMeshType


TPZGenGrid3D::TPZGenGrid3D(const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX, const TPZVec<int> &nelDiv, const MMeshType meshType)
                           :fGmesh{nullptr}, fMinX{minX}, fMaxX{maxX},fNelDiv{nelDiv},fMeshType{meshType}{
    const REAL tol{ZeroTolerance()};
    if(fMinX.NElements() != 3 || fMaxX.NElements() != 3 || fNelDiv.NElements() != 3){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"minX, maxX and nelDiv must have size = 3!\n";
        PZError<<"size(minX) = "<<fMinX.NElements()<<"\n"
               <<"size(maxX) = "<<fMaxX.NElements()<<"\n"
               <<"size(nelDiv) = "<<fNelDiv.NElements()<<"\n";
        DebugStop();
    }
    if(MMeshType_Dimension(fMeshType) != fDim){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"Element type "<<fMeshType<<" is not supported. Aborting...\n";
        DebugStop();
    }
    else if(fNelDiv[0] < 1 || fNelDiv[1] < 1 || fNelDiv[2] < 1){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"The read number of grid divisions is not allowed. The parameters are:\n";
        PZError<<"nel x: "<<fNelDiv[0]<<"\n"<<"nel y: "<<fNelDiv[1]<<"\n"<<"nel z: "<<fNelDiv[2]<<"\n";
        DebugStop();
    }else if(fMaxX[0] < tol || fMaxX[1] < tol || fMaxX[2] < tol){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"Dimensions of grid not allowed. The parameters are:\n";
        PZError<<"max x: "<<fMaxX[0]<<"\n"<<"max y: "<<fMaxX[1]<<"\n"<<"max z: "<<fMaxX[2]<<"\n";
        DebugStop();
    }
    else if(fMaxX[0] < fMinX[0] || fMaxX[1] < fMinX[1] || fMaxX[2] < fMinX[2]){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"Dimensions of grid not allowed. The parameters are:\n";
        PZError<<"min x: "<<fMinX[0]<<"\n"<<"min y: "<<fMinX[1]<<"\n"<<"min z: "<<fMinX[2]<<"\n";
        PZError<<"max x: "<<fMaxX[0]<<"\n"<<"max y: "<<fMaxX[1]<<"\n"<<"max z: "<<fMaxX[2]<<"\n";
        DebugStop();
    }
}

TPZGeoMesh* TPZGenGrid3D::BuildVolumetricElements(const int matIdDomain){
    if(fGmesh){
        PZError<<__PRETTY_FUNCTION__<<" error.\n";
        PZError<<"The geometric mesh is not empty (nullptr). This method should be called only once,";
        PZError<<" before any other TPZGenGrid3D methods are called.\n";
        PZError<<"Aborting...\n";
        DebugStop();
        return nullptr;
    }
    fGmesh = new TPZGeoMesh();
    fGmesh->SetDimension(fDim);
    //create nodes
    const int& nelx = fNelDiv[0];
    const int& nely = fNelDiv[1];
    const int& nelz = fNelDiv[2];
    [&](){
        TPZManVector<REAL,3> nodeCoords(3,0);
        int64_t newindex = -1;
        for(auto iZ = 0; iZ <= nelz; iZ++){
            for(auto iY = 0; iY <= nely; iY++){
                for(auto iX = 0; iX <= nelx; iX++){
                    nodeCoords[0] = fMinX[0] + ((fMaxX[0]-fMinX[0]) * iX)/nelx;
                    nodeCoords[1] = fMinX[1] + ((fMaxX[1]-fMinX[1]) * iY)/nely;
                    nodeCoords[2] = fMinX[2] + ((fMaxX[2]-fMinX[2]) * iZ)/nelz;
                    newindex = fGmesh->NodeVec().AllocateNewElement();
                    fGmesh->NodeVec()[newindex].Initialize(nodeCoords, *fGmesh);
                }
            }
        }
    }();

    //create volumetric elements
    [&](){
        const int nNodes = [&](){
            switch(fMeshType){
                case MMeshType::ETetrahedral: return 4;
                case MMeshType::EPrismatic: return 6;
                case MMeshType::EHexahedral:
                case MMeshType::EPyramidal:
                case MMeshType::EHexaPyrMixed: return 8;
                default:DebugStop(); return -1;
            }
        }();
        TPZManVector<int64_t,8> nodesIdVec(nNodes,-1);
        int interiorNodeCount = (nelx+1)*(nely+1)*(nelz+1);
        for(auto iZ = 0; iZ < nelz; iZ++){
            for(auto iY = 0; iY < nely; iY++){
                for(auto iX = 0; iX < nelx; iX++){
                    // const auto count = iX + iY + iZ;
                    const int permut = (iX + iY + iZ)%2; // Alternates 0 and 1 when moving on either of 3 directions
                    const auto firstNodeId = iZ * (nelx+1) * (nely+1) + iY * (nelx+1) + iX;/*lower left node*/


                    // Lablelling node indices by cube node for readability
                    const int64_t cubenode[8] = {
                                                    firstNodeId,
                                                    firstNodeId + 1, /* +1 moves up a layer in x direction*/
                                                    firstNodeId + 1 + (nelx+1), /* +(nelx+1) moves up a layer in y direction*/
                                                    firstNodeId + (nelx+1),
                                                    firstNodeId                 + (nelx+1) * (nely+1), /* +(nelx+1)*(nely+1) moves up a layer in z direction*/
                                                    firstNodeId + 1             + (nelx+1) * (nely+1),
                                                    firstNodeId + 1 + (nelx+1)  + (nelx+1) * (nely+1),
                                                    firstNodeId + (nelx+1)      + (nelx+1) * (nely+1)
                                                };

                    switch (fMeshType) {
                        case MMeshType::ETetrahedral:
                            nodesIdVec[0] = cubenode[(0+permut)%4];     // (i + permut)%4 rotates the cube 90 degrees around the zeta axis if permut == 1
                            nodesIdVec[1] = cubenode[(1+permut)%4];
                            nodesIdVec[2] = cubenode[(3+permut)%4];
                            nodesIdVec[3] = cubenode[(0+permut)%4 + 4]; // +4 gives the top layer cubenode
                            new TPZGeoElRefPattern<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*fGmesh);
                            nodesIdVec[0] = cubenode[(1+permut)%4 + 4];
                            nodesIdVec[1] = cubenode[(0+permut)%4 + 4];
                            nodesIdVec[2] = cubenode[(2+permut)%4 + 4];
                            nodesIdVec[3] = cubenode[(1+permut)%4];
                            new TPZGeoElRefPattern<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*fGmesh);
                            nodesIdVec[0] = cubenode[(2+permut)%4];
                            nodesIdVec[1] = cubenode[(3+permut)%4];
                            nodesIdVec[2] = cubenode[(1+permut)%4];
                            nodesIdVec[3] = cubenode[(2+permut)%4 + 4];
                            new TPZGeoElRefPattern<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*fGmesh);
                            nodesIdVec[0] = cubenode[(3+permut)%4 + 4];
                            nodesIdVec[1] = cubenode[(2+permut)%4 + 4];
                            nodesIdVec[2] = cubenode[(0+permut)%4 + 4];
                            nodesIdVec[3] = cubenode[(3+permut)%4];
                            new TPZGeoElRefPattern<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*fGmesh);
                            nodesIdVec[0] = cubenode[(0+permut)%4 + 4];
                            nodesIdVec[1] = cubenode[(1+permut)%4];
                            nodesIdVec[2] = cubenode[(3+permut)%4];
                            nodesIdVec[3] = cubenode[(2+permut)%4 + 4];
                            new TPZGeoElRefPattern<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*fGmesh);
                            break;
                        case MMeshType::EHexahedral:
                        case MMeshType::EHexaPyrMixed:
                        case MMeshType::EPyramidal:
                            nodesIdVec[0] = cubenode[0];
                            nodesIdVec[1] = cubenode[1];
                            nodesIdVec[2] = cubenode[2];
                            nodesIdVec[3] = cubenode[3];
                            nodesIdVec[4] = cubenode[4];
                            nodesIdVec[5] = cubenode[5];
                            nodesIdVec[6] = cubenode[6];
                            nodesIdVec[7] = cubenode[7];
                            new TPZGeoElRefPattern<pzgeom::TPZGeoCube>(nodesIdVec,matIdDomain,*fGmesh);
                            break;
                        case MMeshType::EPrismatic:
                            nodesIdVec[0] = cubenode[0];
                            nodesIdVec[1] = cubenode[1];
                            nodesIdVec[2] = cubenode[2];
                            nodesIdVec[3] = cubenode[4];
                            nodesIdVec[4] = cubenode[5];
                            nodesIdVec[5] = cubenode[6];
                            new TPZGeoElRefPattern<pzgeom::TPZGeoPrism>(nodesIdVec,matIdDomain,*fGmesh);

                            nodesIdVec[0] = cubenode[0];
                            nodesIdVec[1] = cubenode[2];
                            nodesIdVec[2] = cubenode[3];
                            nodesIdVec[3] = cubenode[4];
                            nodesIdVec[4] = cubenode[6];
                            nodesIdVec[5] = cubenode[7];
                            new TPZGeoElRefPattern<pzgeom::TPZGeoPrism>(nodesIdVec,matIdDomain,*fGmesh);
                            break;
                        default:
                            DebugStop();
                    }
                }
            }
        }

    }();
    //for meshes with pyramids, some elements need to be refined.
    if(fMeshType == MMeshType::EPyramidal || fMeshType == MMeshType::EHexaPyrMixed){

        const char buf[] =
                "9     7  "
                "37     HexIntoPyramids  "
                "-1.    -1.    -1.  "
                " 1.    -1.    -1.  "
                " 1.     1.    -1.  "
                "-1.     1.    -1.  "
                "-1.    -1.     1.  "
                " 1.    -1.     1.  "
                " 1.     1.     1.  "
                "-1.     1.     1.  "
                " 0.     0.     0.  "
                "7     8     0     1     2     3     4     5     6     7  "
                "5     5     0     1     2     3     8 "
                "5     5     0     4     5     1     8 "
                "5     5     1     5     6     2     8 "
                "5     5     2     6     7     3     8 "
                "5     5     3     7     4     0     8 "
                "5     5     7     6     5     4     8 "
        ;
        std::istringstream str(buf);
        TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
        TPZAutoPointer<TPZRefPattern> refpatFound = gRefDBase.FindRefPattern(refpat);
        if(!refpatFound)
        {
            gRefDBase.InsertRefPattern(refpat);
        }
        else
        {
            refpatFound->SetName(refpat->Name());
        }
        refpat->InsertPermuted();

        fGmesh->BuildConnectivity();
        const int nElem = fGmesh->NElements();
        TPZManVector<TPZGeoEl *,6> sons(6,nullptr);
        for(auto iel = 0; iel < nElem; iel++){
            if((iel + iel/nelx + iel/(nelx*nely))% 2 == 1 && fMeshType == MMeshType::EHexaPyrMixed) continue;
            TPZGeoEl *gel = fGmesh->Element(iel);
            gel->SetRefPattern(refpat);
            gel->Divide(sons);
            gel->RemoveConnectivities();
            const int index = fGmesh->ElementIndex(gel);
            fGmesh->ElementVec()[index] = nullptr;
            delete gel;
            for(auto &son : sons){
                son->SetFather(nullptr);
            }
        }
    }
    fGmesh->BuildConnectivity();
    return fGmesh;
}

TPZGeoMesh *TPZGenGrid3D::BuildBoundaryElements(const int matIdZmin, const int matIdXmin,
                                                const int matIdYmin, const int matIdXmax,
                                                const int matIdYmax, const int matIdZmax) {
    if(!fGmesh){
        PZError<<__PRETTY_FUNCTION__<<" error.\n";
        PZError<<"The geometric mesh is empty (nullptr). Currently this class only supports building";
        PZError<<" boundary elements after the volumetric elements have been created. have you called TPZGenGrid3D::BuildVolumetricElements?\n";
        PZError<<"Aborting...\n";
        DebugStop();
        return nullptr;
    }

    //create boundary elements
    const int& nelx = fNelDiv[0];
    const int& nely = fNelDiv[1];
    const int& nelz = fNelDiv[2];
    TPZManVector<int64_t,8> nodesIdVec(0);
    switch (fMeshType) {
        case MMeshType::EHexahedral:
        case MMeshType::EHexaPyrMixed:
            nodesIdVec.Resize(4);
            break;
        case MMeshType::ETetrahedral:
        case MMeshType::EPrismatic:
            nodesIdVec.Resize(3);
            break;
        default:
            DebugStop();
    }
    //top/bottom
    int aux = 0;
    for(auto iZ = 0; iZ < nelz + 1; iZ+=nelz){
        const auto matIdBoundary = iZ == 0 ? matIdZmin : matIdZmax;
        for(auto iY = 0; iY < nely; iY++){
            for(auto iX = 0; iX < nelx; iX++){
                const auto firstNodeId = iZ * (nelx+1) * (nely+1) + iY * (nelx+1) + iX;/*lower left node*/
                const auto count = iX+iY+iZ;
                switch (fMeshType) {
                    case MMeshType::ETetrahedral:
                        aux = count % 2 ? 1 : 0;
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + 1;
                        nodesIdVec[2] = firstNodeId + (nelx+1) + aux;
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        aux = count % 2 ? 0 : 1;
                        nodesIdVec[0] = firstNodeId + aux;
                        nodesIdVec[1] = firstNodeId + (nelx+1);
                        nodesIdVec[2] = firstNodeId + (nelx+1) + 1;
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    case MMeshType::EHexahedral:
                    case MMeshType::EHexaPyrMixed:
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + 1;
                        nodesIdVec[2] = firstNodeId + (nelx+1) + 1;
                        nodesIdVec[3] = firstNodeId + (nelx+1);
                        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    case MMeshType::EPrismatic:
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + 1;
                        nodesIdVec[2] = firstNodeId + (nelx+1) + 1;
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + (nelx+1);
                        nodesIdVec[2] = firstNodeId + (nelx+1) + 1;
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    default:
                        DebugStop();
                }
            }
        }
    }
    switch (fMeshType) {
        case MMeshType::ETetrahedral:
            nodesIdVec.Resize(3);
            break;
        case MMeshType::EHexaPyrMixed:
        case MMeshType::EHexahedral:
        case MMeshType::EPrismatic:
            nodesIdVec.Resize(4);
            break;
        default:
            DebugStop();
    }
    //left/right
    for(auto iZ = 0; iZ < nelz; iZ++){
        for(auto iY = 0; iY < nely+1; iY+=nely){
            const auto matIdBoundary = iY == 0 ? matIdYmin : matIdYmax;
            for(auto iX = 0; iX < nelx; iX++){
                const auto firstNodeId = iZ * (nelx+1) * (nely+1) + iY * (nelx+1) + iX;/*lower left node*/
                const auto count = iX+iY+iZ;
                switch (fMeshType) {
                    case MMeshType::ETetrahedral:
                        aux = count % 2 ? 1 : 0;
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + 1;
                        nodesIdVec[2] = firstNodeId + (nelx+1) * (nely+1) + aux;
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        aux = count % 2 ? 0 : 1;
                        nodesIdVec[0] = firstNodeId + aux;
                        nodesIdVec[1] = firstNodeId + (nelx+1) * (nely+1)+1;
                        nodesIdVec[2] = firstNodeId + (nelx+1) * (nely+1);
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    case MMeshType::EHexahedral:
                    case MMeshType::EHexaPyrMixed:
                    case MMeshType::EPrismatic:
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + 1;
                        nodesIdVec[2] = firstNodeId + (nelx+1) * (nely+1) + 1 ;
                        nodesIdVec[3] = firstNodeId + (nelx+1) * (nely+1);
                        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    default:
                        DebugStop();
                }
            }
        }
    }
    //front/back
    for(auto iZ = 0; iZ < nelz; iZ++){
        for(auto iY = 0; iY < nely; iY++){
            for(auto iX = 0; iX < nelx+1; iX+=nelx){
                const auto matIdBoundary = iX == 0 ? matIdXmin : matIdXmax;
                const auto firstNodeId = iZ * (nelx+1) * (nely+1) + iY * (nelx+1) + iX;/*lower left node*/
                const auto count = iX+iY+iZ;
                switch (fMeshType) {
                    case MMeshType::ETetrahedral:
                        aux = count % 2 ? nelx+1 : 0;
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + (nelx+1);
                        nodesIdVec[2] = firstNodeId + (nelx+1) * (nely+1) + aux;
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        aux = count % 2 ? 0 : nelx+1;
                        nodesIdVec[0] = firstNodeId + aux;
                        nodesIdVec[1] = firstNodeId + (nelx+1) * (nely+1) + (nelx+1);
                        nodesIdVec[2] = firstNodeId + (nelx+1) * (nely+1);
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    case MMeshType::EHexahedral:
                    case MMeshType::EHexaPyrMixed:
                    case MMeshType::EPrismatic:
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + (nelx+1);
                        nodesIdVec[2] = firstNodeId + (nelx+1) * (nely+1) + (nelx+1);
                        nodesIdVec[3] = firstNodeId + (nelx+1) * (nely+1);
                        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    default:
                        DebugStop();
                }
            }
        }
    }

    fGmesh->BuildConnectivity();
    return fGmesh;
}
