//
// Created by Francisco Teixeira Orlandini on 19/12/19.
//

#include "TPZGenGrid3D.h"

#include "pzgmesh.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "tpzgeoelrefpattern.h"


TPZGenGrid3D::TPZGenGrid3D(const REAL minX, const REAL minY, const REAL minZ,
                           const REAL maxX, const REAL maxY, const REAL maxZ,
                           const int nelx, const int nely, const int nelz, const MElementType elType)
                           :fGmesh{nullptr}, fMinX{minX}, fMinY{minY}, fMinZ{minZ},
                           fMaxX{maxX}, fMaxY{maxY}, fMaxZ{maxZ},
                           fNelx{nelx},fNely{nely},fNelz{nelz},fEltype{elType}{
    const REAL tol{ZeroTolerance()};
    if(fEltype != ETetraedro && fEltype != ECube && fEltype != EPrisma){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"Element type "<<fEltype<<" is not supported. Aborting...\n";
        DebugStop();
    }
    else if(fNelx < 1 || fNely < 1 || fNelz < 1){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"The read number of grid divisions is not allowed. The parameters are:\n";
        PZError<<"nel x: "<<fNelx<<"\n"<<"nel y: "<<fNely<<"\n"<<"nel z: "<<fNelz<<"\n";
        DebugStop();
    }else if(fMaxX < tol || fMaxY < tol || fMaxZ < tol){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"Dimensions of grid not allowed. The parameters are:\n";
        PZError<<"max x: "<<fMaxX<<"\n"<<"max y: "<<fMaxY<<"\n"<<"max z: "<<fMaxZ<<"\n";
        DebugStop();
    }
    else if(fMaxX < fMinX || fMaxY < fMinY || fMaxZ < fMinZ){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"Dimensions of grid not allowed. The parameters are:\n";
        PZError<<"min x: "<<fMinX<<"\n"<<"min y: "<<fMinY<<"\n"<<"min z: "<<fMinZ<<"\n";
        PZError<<"max x: "<<fMaxX<<"\n"<<"max y: "<<fMaxY<<"\n"<<"max z: "<<fMaxZ<<"\n";
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
    [&](){
        TPZManVector<REAL,3> nodeCoords(3,0);
        int64_t newindex = -1;
        for(auto iZ = 0; iZ <= fNelz; iZ++){
            for(auto iY = 0; iY <= fNely; iY++){
                for(auto iX = 0; iX <= fNelx; iX++){
                    nodeCoords[0] = fMinX + (fMaxX * iX)/fNelx;
                    nodeCoords[1] = fMinY + (fMaxY * iY)/fNely;
                    nodeCoords[2] = fMinZ + (fMaxZ * iZ)/fNelz;
                    newindex = fGmesh->NodeVec().AllocateNewElement();
                    fGmesh->NodeVec()[newindex].Initialize(nodeCoords, *fGmesh);
                }
            }
        }
    }();

    //create volumetric elements
    [&](){
        TPZManVector<int64_t,8> nodesIdVec(MElementType_NNodes(fEltype),-1);
        for(auto iZ = 0; iZ < fNelz; iZ++){
            for(auto iY = 0; iY < fNely; iY++){
                for(auto iX = 0; iX < fNelx; iX++){
                    const auto count = iX + iY + iZ;
                    const auto firstNodeId = iZ * (fNelx+1) * (fNely+1) + iY * (fNelx+1) + iX;/*lower left node*/
                    switch (fEltype) {
                        case ETetraedro:
                            nodesIdVec[0] = firstNodeId;
                            nodesIdVec[1] = firstNodeId + 1;
                            nodesIdVec[2] = firstNodeId + (fNelx+1) + count % 2;
                            nodesIdVec[3] = firstNodeId + (fNelx+1) * (fNely+1) + count % 2;
                            new TPZGeoElRefPattern<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*fGmesh);
                            nodesIdVec[0] = firstNodeId + 1 - count % 2;
                            nodesIdVec[1] = firstNodeId + (fNelx+1) * (fNely+1) + 1;
                            nodesIdVec[2] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1) + 1 -  count % 2;
                            nodesIdVec[3] = firstNodeId + (fNelx+1) * (fNely+1);
                            new TPZGeoElRefPattern<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*fGmesh);
                            nodesIdVec[0] = firstNodeId + 1 - count % 2;
                            nodesIdVec[1] = firstNodeId + (fNelx+1) + 1;
                            nodesIdVec[2] = firstNodeId + (fNelx+1);
                            nodesIdVec[3] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1) + 1 - count % 2;
                            new TPZGeoElRefPattern<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*fGmesh);
                            nodesIdVec[0] = firstNodeId + (fNelx+1) + count % 2;
                            nodesIdVec[1] = firstNodeId + (fNelx+1) * (fNely+1)+ count % 2;
                            nodesIdVec[2] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1);
                            nodesIdVec[3] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1) + 1;
                            new TPZGeoElRefPattern<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*fGmesh);
                            nodesIdVec[0] = firstNodeId + 1 - count % 2;
                            nodesIdVec[1] = firstNodeId + (fNelx+1) + count % 2;
                            nodesIdVec[2] = firstNodeId + (fNelx+1) * (fNely+1) + count % 2;
                            nodesIdVec[3] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1) + 1 - count % 2;
                            new TPZGeoElRefPattern<pzgeom::TPZGeoTetrahedra>(nodesIdVec,matIdDomain,*fGmesh);
                            break;
                        case ECube:
                            nodesIdVec[0] = firstNodeId;
                            nodesIdVec[1] = firstNodeId + 1;
                            nodesIdVec[2] = firstNodeId + (fNelx+1) + 1;
                            nodesIdVec[3] = firstNodeId + (fNelx+1);
                            nodesIdVec[4] = firstNodeId + (fNelx+1) * (fNely+1);
                            nodesIdVec[5] = firstNodeId + (fNelx+1) * (fNely+1) + 1;
                            nodesIdVec[6] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1) + 1;
                            nodesIdVec[7] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1);
                            new TPZGeoElRefPattern<pzgeom::TPZGeoCube>(nodesIdVec,matIdDomain,*fGmesh);
                            break;
                        case EPrisma:
                            nodesIdVec[0] = firstNodeId;
                            nodesIdVec[1] = firstNodeId + 1;
                            nodesIdVec[2] = firstNodeId + (fNelx+1) + 1;
                            nodesIdVec[3] = firstNodeId + (fNelx+1) * (fNely+1);
                            nodesIdVec[4] = firstNodeId + (fNelx+1) * (fNely+1) + 1;
                            nodesIdVec[5] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1) + 1;
                            new TPZGeoElRefPattern<pzgeom::TPZGeoPrism>(nodesIdVec,matIdDomain,*fGmesh);

                            nodesIdVec[0] = firstNodeId;
                            nodesIdVec[1] = firstNodeId + (fNelx+1);
                            nodesIdVec[2] = firstNodeId + (fNelx+1) + 1;
                            nodesIdVec[3] = firstNodeId + (fNelx+1) * (fNely+1);
                            nodesIdVec[4] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1);
                            nodesIdVec[5] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1) + 1;
                            new TPZGeoElRefPattern<pzgeom::TPZGeoPrism>(nodesIdVec,matIdDomain,*fGmesh);
                            break;
                        default:
                            DebugStop();
                    }
                }
            }
        }

    }();
    fGmesh->BuildConnectivity();
    return fGmesh;
}

TPZGeoMesh *TPZGenGrid3D::BuildBoundaryElements(const int matIdXmin, const int matIdXmax,
                                                const int matIdYmin, const int matIdYmax,
                                                const int matIdZmin, const int matIdZmax) {
    if(!fGmesh){
        PZError<<__PRETTY_FUNCTION__<<" error.\n";
        PZError<<"The geometric mesh is empty (nullptr). Currently this class only supports building";
        PZError<<" boundary elements after the volumetric elements have been created. have you called TPZGenGrid3D::BuildVolumetricElements?\n";
        PZError<<"Aborting...\n";
        DebugStop();
        return nullptr;
    }

    //create boundary elements

    TPZManVector<int64_t,8> nodesIdVec(MElementType_NNodes(fEltype),-1);
    switch (fEltype) {
        case ECube:
            nodesIdVec.Resize(4);
            break;
        case ETetraedro:
        case EPrisma:
            nodesIdVec.Resize(3);
            break;
        default:
            DebugStop();
    }
    //top/bottom
    int aux = 0;
    for(auto iZ = 0; iZ < fNelz + 1; iZ+=fNelz){
        const auto matIdBoundary = iZ == 0 ? matIdZmin : matIdZmax;
        for(auto iY = 0; iY < fNely; iY++){
            for(auto iX = 0; iX < fNelx; iX++){
                const auto firstNodeId = iZ * (fNelx+1) * (fNely+1) + iY * (fNelx+1) + iX;/*lower left node*/
                const auto count = iX+iY+iZ;
                switch (fEltype) {
                    case ETetraedro:
                        aux = count % 2 ? 1 : 0;
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + 1;
                        nodesIdVec[2] = firstNodeId + (fNelx+1) + aux;
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        aux = count % 2 ? 0 : 1;
                        nodesIdVec[0] = firstNodeId + aux;
                        nodesIdVec[1] = firstNodeId + (fNelx+1);
                        nodesIdVec[2] = firstNodeId + (fNelx+1) + 1;
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    case ECube:
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + 1;
                        nodesIdVec[2] = firstNodeId + (fNelx+1) + 1;
                        nodesIdVec[3] = firstNodeId + (fNelx+1);
                        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    case EPrisma:
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + 1;
                        nodesIdVec[2] = firstNodeId + (fNelx+1) + 1;
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + (fNelx+1);
                        nodesIdVec[2] = firstNodeId + (fNelx+1) + 1;
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    default:
                        DebugStop();
                }
            }
        }
    }
    switch (fEltype) {
        case ETetraedro:
            nodesIdVec.Resize(3);
            break;
        case ECube:
        case EPrisma:
            nodesIdVec.Resize(4);
            break;
        default:
            DebugStop();
    }
    //left/right
    for(auto iZ = 0; iZ < fNelz; iZ++){
        for(auto iY = 0; iY < fNely+1; iY+=fNely){
            const auto matIdBoundary = iY == 0 ? matIdYmin : matIdYmax;
            for(auto iX = 0; iX < fNelx; iX++){
                const auto firstNodeId = iZ * (fNelx+1) * (fNely+1) + iY * (fNelx+1) + iX;/*lower left node*/
                const auto count = iX+iY+iZ;
                switch (fEltype) {
                    case ETetraedro:
                        aux = count % 2 ? 1 : 0;
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + 1;
                        nodesIdVec[2] = firstNodeId + (fNelx+1) * (fNely+1) + aux;
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        aux = count % 2 ? 0 : 1;
                        nodesIdVec[0] = firstNodeId + aux;
                        nodesIdVec[1] = firstNodeId + (fNelx+1) * (fNely+1)+1;
                        nodesIdVec[2] = firstNodeId + (fNelx+1) * (fNely+1);
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    case ECube:
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + 1;
                        nodesIdVec[2] = firstNodeId + (fNelx+1) * (fNely+1) + 1 ;
                        nodesIdVec[3] = firstNodeId + (fNelx+1) * (fNely+1);
                        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    case EPrisma:
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + 1;
                        nodesIdVec[2] = firstNodeId + (fNelx+1) * (fNely+1) + 1;
                        nodesIdVec[3] = firstNodeId + (fNelx+1) * (fNely+1);
                        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    default:
                        DebugStop();
                }
            }
        }
    }
    //front/back
    for(auto iZ = 0; iZ < fNelz; iZ++){
        for(auto iY = 0; iY < fNely; iY++){
            for(auto iX = 0; iX < fNelx+1; iX+=fNelx){
                const auto matIdBoundary = iX == 0 ? matIdXmin : matIdXmax;
                const auto firstNodeId = iZ * (fNelx+1) * (fNely+1) + iY * (fNelx+1) + iX;/*lower left node*/
                const auto count = iX+iY+iZ;
                switch (fEltype) {
                    case ETetraedro:
                        aux = count % 2 ? fNelx+1 : 0;
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + (fNelx+1);
                        nodesIdVec[2] = firstNodeId + (fNelx+1) * (fNely+1) + aux;
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        aux = count % 2 ? 0 : fNelx+1;
                        nodesIdVec[0] = firstNodeId + aux;
                        nodesIdVec[1] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1);
                        nodesIdVec[2] = firstNodeId + (fNelx+1) * (fNely+1);
                        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    case ECube:
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + (fNelx+1);
                        nodesIdVec[2] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1);
                        nodesIdVec[3] = firstNodeId + (fNelx+1) * (fNely+1);
                        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,matIdBoundary,*fGmesh);
                        break;
                    case EPrisma:
                        nodesIdVec[0] = firstNodeId;
                        nodesIdVec[1] = firstNodeId + (fNelx+1);
                        nodesIdVec[2] = firstNodeId + (fNelx+1) * (fNely+1) + (fNelx+1);
                        nodesIdVec[3] = firstNodeId + (fNelx+1) * (fNely+1);
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
