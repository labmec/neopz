#include "TPZGeoMeshTools.h"
#include "TPZRefPattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"
#include "TPZAcademicGeoMesh.h"

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
TPZGeoMeshTools::CreateGeoMeshOnGrid(int dim, const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX, const TPZVec<int> &matids,
                                     const TPZVec<int> nDivs, MMeshType meshType, bool createBoundEls) {
#ifdef PZDEBUG
    const REAL tol{ZeroTolerance()};
    if(dim !=3 && dim!= 2 && dim != 1){
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
    if(dim == 1)
    {
        if(nDivs[0] < 1){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"The read number of grid divisions is not allowed. The parameters are:\n";
            PZError<<"nel x: "<<nDivs[0]<<"\n"<<"nel y: "<<nDivs[1]<<"\n";
            DebugStop();
        }else if(maxX[0] < tol){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"Dimensions of grid not allowed. The parameters are:\n";
            PZError<<"max x: "<<maxX[0]<<"\n"<<"max y: "<<maxX[1]<<"\n";
            DebugStop();
        }
        else if(maxX[0] < minX[0]){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"Dimensions of grid not allowed. The parameters are:\n";
            PZError<<"min x: "<<minX[0]<<"\n"<<"min y: "<<minX[1]<<"\n";
            PZError<<"max x: "<<maxX[0]<<"\n"<<"max y: "<<maxX[1]<<"\n";
            DebugStop();
        }
    } else if(dim == 2){
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
            case 1:{
                int nelem = 1 << (nDivs[0]-1);
                TPZAcademicGeoMesh academic(nelem,TPZAcademicGeoMesh::MMeshType::EOned);
                return academic.CreateGeoMesh();
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
            default:
                return (TPZGeoMesh*)nullptr;
        }
    }();


}

    /**
     * This method creates a TPZGeoMesh of dim 3(2) based on a cuboid(rectangle) aligned with the x,y and z axis (x and y).
     * @param meshsize : size of the domain
     * @param matids Material identifiers. If creating boundaries, it must have size dim*2 + 1.
     * Their order follows the side numbering of the hexahedral element (quadrilateral element)
     * @param nDivs Number of elements in each direction
     * @param meshType What type of elements to create
     * @param createBoundEls Whether to create the boundary elements
     * @return Pointer to the TPZGeoMesh that has been created
     */
TPZGeoMesh * TPZGeoMeshTools::CreateGeoMeshOnGrid(REAL meshsize,
            const TPZVec<int> &matids, const TPZVec<int> nDivs, MMeshType meshType, bool createBoundEls)
{
    TPZManVector<REAL,3> x0(3,0.),x1(3,meshsize);
    int dim = 1;
    switch (meshType) {
        case MMeshType::EOneDimensional:
            dim = 1;
            break;
        case MMeshType::EQuadrilateral:
        case MMeshType::ETriangular:
            dim = 2;
            break;
        default:
            dim = 3;
            break;
    }
    for(int i=dim; i<3; i++) x1[i] = 0.;
    TPZManVector<int> divloc(nDivs);
    divloc.Resize(dim);
    return CreateGeoMeshOnGrid(dim, x0, x1, matids, divloc, meshType, createBoundEls);
}
