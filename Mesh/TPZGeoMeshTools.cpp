#include "TPZGeoMeshTools.h"
#include "TPZRefPattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"
#include "tpzpyramid.h"
#include "pzgeoelbc.h"

void TPZGeoMeshTools::DividePyramidsIntoTetra(TPZGeoMesh *gmesh) {
    const char buf[] =
            "3     5  "
            "37     PyramidsIntoTetra  "
            "-1.    -1.    0.  "
            " 1.    -1.    0.  "
            " 1.     1.    0.  "
            "-1.     1.    0.  "
            " 0.     0.    1.  "
            "7     5     0     1     2     3     4 "
            "4     4     0     1     2     4 "
            "4     4     0     2     3     4 "
    ;
    std::istringstream str(buf);
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
    TPZAutoPointer<TPZRefPattern> refpatFound = gRefDBase.FindRefPattern(refpat);
    if(!refpatFound){
        gRefDBase.InsertRefPattern(refpat);
    }
    else{
        refpatFound->SetName(refpat->Name());
    }
    refpat->InsertPermuted();

    TPZManVector<TPZGeoEl *, 2> sons;
    for(auto &gel : gmesh->ElementVec()){
        if(gel->Type() == EPiramide || gel->NSubElements() == 0){
            gel->SetRefPattern(refpat);
            gel->Divide(sons);
        }
    }
    gmesh->BuildConnectivity();
}

TPZGeoMesh *
TPZGeoMeshTools::CreateGeoMesh1D(const REAL minX, const REAL maxX, const int nEls,
                const TPZVec<int> &matids, bool createBoundEls){
#ifdef PZDEBUG
    if(minX > maxX){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Invalid values of minX and maxX\n";
        PZError<<"minX: "<<minX<<" maxX: "<<maxX;
        DebugStop();
    }
    if(nEls < 1){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Invalid value of nEls\n";
        PZError<<"nEls: "<<nEls;
        DebugStop();
    }
#endif
    if((createBoundEls && matids.size()!=3)||
       (!createBoundEls && matids.size()!=1)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Wrong size of matids vec!\n";
        DebugStop();
    }
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    const int matId{matids[0]};
    const auto nnodes = nEls + 1;
    const auto elSize = (maxX-minX)/nEls;
    gmesh->NodeVec().Resize(nnodes);
    for (auto i = 0 ; i < nnodes; i++)
      {
        const REAL pos = minX + i*elSize;
        TPZManVector <REAL,3> coord= {pos,0.,0.};
        coord[0] = pos;
        gmesh->NodeVec()[i].SetCoord(coord);
        gmesh->NodeVec()[i].SetNodeId(i);
      }
    TPZManVector <int64_t,2> nodeVec(2);
    int64_t id;
    for (auto iel = 0; iel < nEls; iel++)
      {
        nodeVec={iel,iel+1};
        gmesh->CreateGeoElement(EOned, nodeVec, matId, id);
        gmesh->ElementVec()[id];
      }
    if(createBoundEls){
        const int bc0{matids[1]},bc1{matids[2]};
        TPZManVector<int64_t,1> nodeVec = {0};
        gmesh->CreateGeoElement(EPoint, nodeVec, bc0, id);
        nodeVec = {nEls};
        gmesh->CreateGeoElement(EPoint, nodeVec, bc1, id);
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

TPZGeoMesh *
TPZGeoMeshTools::CreateGeoMeshOnGrid(int dim, const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX, const TPZVec<int> &matids,
                                     const TPZVec<int> nDivs, MMeshType meshType, bool createBoundEls) {
#ifdef PZDEBUG
    const REAL tol{ZeroTolerance()};
    if(dim !=3 && dim!= 2){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"Dimension = "<<dim<<" is not supported. Aborting...\n";
        DebugStop();
    }
    if(MMeshType_Dimension(meshType) != dim){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"Element type "<<meshType<<" is not supported. Aborting...\n";
        DebugStop();
    }
    if(minX.NElements() != 3 || maxX.NElements() != 3){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"minX and maxX must have size = 3!\n";
        PZError<<"size(minX) = "<<minX.NElements()<<"\n"
               <<"size(maxX) = "<<maxX.NElements()<<"\n";
        DebugStop();
    }
    if(nDivs.NElements() != dim){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"size(nDivs) != dim !\n";
        PZError<<"dim = "<<minX.NElements()<<"\n"
               <<"size(nDivs) = "<<nDivs.NElements()<<"\n";
        DebugStop();
    }
    if(matids.size() != dim*2 + 1 && createBoundEls){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"The number of material ids must be equal to the number of boundaries + 1\n"
               <<"# of matids: "<<matids.size()<<"\n"
               <<"# of boundaries + 1: "<<dim*2+1<<"\n";
        DebugStop();
    } else if(matids.size() != 1 && !createBoundEls){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"When no boundary is created, the number of material ids must be equal to one\n"
               <<"# of matids: "<<matids.size()<<"\n";
        DebugStop();
    }
    if(dim == 2){
        if(nDivs[0] < 1 || nDivs[1] < 1){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"The read number of grid divisions is not allowed. The parameters are:\n";
            PZError<<"nel x: "<<nDivs[0]<<"\n"<<"nel y: "<<nDivs[1]<<"\n";
            DebugStop();
        }else if(maxX[0] < tol || maxX[1] < tol){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"Dimensions of grid not allowed. The parameters are:\n";
            PZError<<"max x: "<<maxX[0]<<"\n"<<"max y: "<<maxX[1]<<"\n";
            DebugStop();
        }
        else if(maxX[0] < minX[0] || maxX[1] < minX[1]){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"Dimensions of grid not allowed. The parameters are:\n";
            PZError<<"min x: "<<minX[0]<<"\n"<<"min y: "<<minX[1]<<"\n";
            PZError<<"max x: "<<maxX[0]<<"\n"<<"max y: "<<maxX[1]<<"\n";
            DebugStop();
        }
    }
    else{//dim == 3
        if(nDivs[0] < 1 || nDivs[1] < 1 || nDivs[2] < 1){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"The read number of grid divisions is not allowed. The parameters are:\n";
            PZError<<"nel x: "<<nDivs[0]<<"\n"<<"nel y: "<<nDivs[1]<<"\n"<<"nel z: "<<nDivs[2]<<"\n";
            DebugStop();
        }else if(maxX[0] < tol || maxX[1] < tol || maxX[2] < tol){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"Dimensions of grid not allowed. The parameters are:\n";
            PZError<<"max x: "<<maxX[0]<<"\n"<<"max y: "<<maxX[1]<<"\n"<<"max z: "<<maxX[2]<<"\n";
            DebugStop();
        }
        else if(maxX[0] < minX[0] || maxX[1] < minX[1] || maxX[2] < minX[2]){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"Dimensions of grid not allowed. The parameters are:\n";
            PZError<<"min x: "<<minX[0]<<"\n"<<"min y: "<<minX[1]<<"\n"<<"min z: "<<minX[2]<<"\n";
            PZError<<"max x: "<<maxX[0]<<"\n"<<"max y: "<<maxX[1]<<"\n"<<"max z: "<<maxX[2]<<"\n";
            DebugStop();
        }
    }
#endif
    return [&](){
        switch(dim){
            case 2:{
                TPZGeoMesh *gmesh = new TPZGeoMesh();
                TPZGenGrid2D gengrid(nDivs, minX, maxX);
                gengrid.SetElementType(meshType);
                gengrid.Read(gmesh, matids[0]);
                if(createBoundEls){
                    gengrid.SetBC(gmesh, 4, matids[1]);
                    gengrid.SetBC(gmesh, 5, matids[2]);
                    gengrid.SetBC(gmesh, 6, matids[3]);
                    gengrid.SetBC(gmesh, 7, matids[4]);
                }
                gmesh->BuildConnectivity();
                return gmesh;
            }
            case 3:{
                TPZGeoMesh *gmesh = nullptr;
                TPZGenGrid3D genGrid3D(minX,maxX,nDivs,meshType);
                gmesh = genGrid3D.BuildVolumetricElements(matids[0]);
                if(createBoundEls){
                    gmesh = genGrid3D.BuildBoundaryElements(matids[1],matids[2],matids[3],matids[4],matids[5],matids[6]);
                }

                gmesh->BuildConnectivity();
                return gmesh;
            }
            default:{
                PZError<<__PRETTY_FUNCTION__;
                PZError<<"Invalid dimension, returning nullptr\n";
                return (TPZGeoMesh*)nullptr;
            }
        }
    }();

}

TPZGeoMesh *TPZGeoMeshTools::CreateGeoMeshSingleEl(const MMeshType meshType,
                                                   const int matid,
                                                   const bool createBoundEls,
                                                   const int matidbc)
{
  auto gmesh = [=]() ->TPZGeoMesh* {
    switch (meshType) {
    case MMeshType::ETriangular:
      return CreateGeoMeshSingleElT<pztopology::TPZTriangle>(
          matid, createBoundEls, matidbc);
    case MMeshType::EQuadrilateral:
      return CreateGeoMeshSingleElT<pztopology::TPZQuadrilateral>(
          matid, createBoundEls, matidbc);
    case MMeshType::ETetrahedral:
      return CreateGeoMeshSingleElT<pztopology::TPZTetrahedron>(
          matid, createBoundEls, matidbc);
    case MMeshType::EPyramidal:
      return CreateGeoMeshSingleElT<pztopology::TPZPyramid>(
          matid, createBoundEls, matidbc);
    case MMeshType::EPrismatic:
      return CreateGeoMeshSingleElT<pztopology::TPZPrism>(
          matid, createBoundEls, matidbc);
    case MMeshType::EHexahedral:
      return CreateGeoMeshSingleElT<pztopology::TPZCube>(
          matid, createBoundEls, matidbc);
    case MMeshType::EHexaPyrMixed:
      return nullptr;
    case MMeshType::ENoType:
      return nullptr;
    }
  }();
  return gmesh;
}

template <class TGEOM>
TPZGeoMesh *TPZGeoMeshTools::CreateGeoMeshSingleElT(const int matid,
                                                    const bool createBoundEls,
                                                    const int matidbc)
{
  constexpr auto nNodes{TGEOM::NCornerNodes};
  constexpr auto dim{TGEOM::Dimension};
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  gmesh->SetDimension(dim);
  // Auxiliary vector to store coordinates:
  TPZManVector<REAL, 3> coords(3, 0.);
  for (int iNode = 0; iNode < nNodes; iNode++) {
    TGEOM::ParametricDomainNodeCoord(iNode, coords);
    coords.Resize(3);
    auto newindex = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[newindex].Initialize(coords, *gmesh);
  }
  TPZManVector<int64_t, TGEOM::NCornerNodes> nodesIdVec(TGEOM::NCornerNodes,
                                                        -1);
  for (int i = 0; i < TGEOM::NCornerNodes; i++)
    nodesIdVec[i] = i;
  int64_t index{-1};
  TPZGeoEl *gel =
      gmesh->CreateGeoElement(TGEOM::Type(), nodesIdVec, matid, index, 0);
  gmesh->BuildConnectivity();

  if (createBoundEls) {
    const auto fside = gel->FirstSide(dim - 1);
    const auto lastside = fside + gel->NSides(dim - 1);
    for (int iside = fside; iside < lastside; iside++)
      TPZGeoElBC(gel, iside, matidbc);
  }
  return gmesh;
}