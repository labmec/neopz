
#include <pzgeoelrefless.h>
#include <tpzarc3d.h>
#include <TPZVTKGeoMesh.h>
#include <pzanalysis.h>
#include <pzgeotetrahedra.h>
#include <TPZGeoCube.h>
#include <pzgeopyramid.h>
#include <TPZTriangleSphere.h>
#include <TPZQuadSphere.h>
#include <tpzquadraticquad.h>
#include <tpzquadratictrig.h>
#include <tpzquadraticline.h>
#include "tpzgeoblend.h"
#include "pzgeoquad.h"
#include "TPZGeoLinear.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzgeoelmapped.h"

namespace blendtest{
    void CreateGeoMesh2D(TPZGeoMesh *&gmesh, int nDiv, bool printGMesh, std::string prefix);
    template <class TGeo>
    void CreateGeoMesh2DNew(TPZGeoMesh *&gmesh, int nDiv, bool printGMesh, std::string prefix);
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
    bool run3d = false;

    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    if(run3d){
        gRefDBase.InitializeUniformRefPattern(ETetraedro);
        gRefDBase.InitializeUniformRefPattern(EPiramide);
        gRefDBase.InitializeUniformRefPattern(EPrisma);
        gRefDBase.InitializeUniformRefPattern(ECube);
    }


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
//        CreateGeoMesh2D(gmesh, nDiv, printGMesh, prefix);
//        delete gmesh;

        gmesh = nullptr;
//        CreateGeoMesh2D(gmesh, nDiv, printGMesh, prefix);
        CreateGeoMesh2DNew<pzgeom::TPZGeoTriangle>(gmesh, nDiv, printGMesh, prefix);
        delete gmesh;

        gmesh = nullptr;
        CreateGeoMesh2DNew<pzgeom::TPZGeoQuad>(gmesh, nDiv, printGMesh, prefix);
        delete gmesh;
    }

//    pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>::fUseNewX = true;
//    CreateGeoMesh2D(gmesh, printGMesh, nDiv);
//    delete gmesh;


}

namespace blendtest {
    template <class TGeo>
    void CreateGeoMesh2DNew(TPZGeoMesh *&gmesh, int nDiv, bool printGMesh, std::string prefix) {
        auto elType = TGeo::Type();
        switch(elType){
            case ETriangle:
            case EQuadrilateral:
                break;
            default:
                DebugStop();
        }
        const int64_t nNodesVolEl = TGeo::NNodes;
        const int64_t dim = TGeo::Dimension;
        gmesh = new TPZGeoMesh();


        TPZManVector<REAL,3> coordsOffset(3,0);
        coordsOffset[0] = 2;

        TPZManVector<int64_t,4> nodesIdArcsVec(nNodesVolEl,-1);
        {
            TPZVec<REAL> coord(3, 0.);
            const REAL xIni = 1;
            const REAL yIni = 0;

            ///CREATE NODES FOR QUADRATIC ELEMENT
            for (int64_t i = 0; i < nNodesVolEl; i++) {
                const REAL theta = ((REAL)i/nNodesVolEl) *2*M_PI;
                coord[0] = cos(theta) * xIni - sin(theta) * yIni;
                coord[1] = sin(theta) * xIni + cos(theta) * yIni;
                coord[2] = 0;
//            std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
                const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
                gmesh->NodeVec()[newindex].Initialize(coord,*gmesh);
            }//quadraticElCornerNodes

            REAL thetaOff = (1./nNodesVolEl) *M_PI;
            for (int64_t i = 0; i < nNodesVolEl; i++) {
                const REAL theta = ((REAL)i/nNodesVolEl) *2*M_PI + thetaOff;
                coord[0] = cos(theta) * xIni - sin(theta) * yIni;
                coord[1] = sin(theta) * xIni + cos(theta) * yIni;
                coord[2] = 0;
//            std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
                const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
                gmesh->NodeVec()[newindex].Initialize(coord,*gmesh);
            }//quadraticElMidSideNodes

            ///CREATE NODES FOR BLEND ELEMENT
            for (int64_t i = 0; i < nNodesVolEl; i++) {
                const REAL theta = ((REAL)i/nNodesVolEl) *2*M_PI;
                coord[0] = cos(theta) * xIni - sin(theta) * yIni + coordsOffset[0];
                coord[1] = sin(theta) * xIni + cos(theta) * yIni + coordsOffset[1];
                coord[2] = 0 + coordsOffset[2];
//            std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
                const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
                gmesh->NodeVec()[newindex].Initialize(coord,*gmesh);
            }//blendElCornerNodes

            for (int64_t i = 0; i < nNodesVolEl; i++) {
                const REAL theta = ((REAL)i/nNodesVolEl) *2*M_PI + thetaOff;
                coord[0] = cos(theta) * xIni - sin(theta) * yIni + coordsOffset[0];
                coord[1] = sin(theta) * xIni + cos(theta) * yIni + coordsOffset[1];
                coord[2] = 0 + coordsOffset[2];
//            std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
                const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
                gmesh->NodeVec()[newindex].Initialize(coord,*gmesh);
                nodesIdArcsVec[i] = gmesh->NodeVec()[newindex].Id();
            }//blendElMidSideNodes
        }
        const int matIdVol = 1, matIdSphere = 2;
        TPZManVector<int64_t,8> nodesIdVec(1);
        int64_t elId = 0;

        ///////CREATE MAX DIM ELEMENTS
        nodesIdVec.Resize(nNodesVolEl * 2);
        for(int i = 0; i < nodesIdVec.size(); i++ ) nodesIdVec[i] = i;

        TPZGeoEl *quadraticEl = nullptr, *blendEl = nullptr;
        switch(elType){
            case ETriangle:
                quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticTrig>(elId,nodesIdVec,matIdVol,*gmesh);
                nodesIdVec.Resize(nNodesVolEl);
                for(int i = 0; i < nodesIdVec.size(); i++ ) nodesIdVec[i] = i + 2 * nNodesVolEl;
                blendEl = new TPZGeoElRefLess<pzgeom::TPZGeoBlend<TGeo>>(elId,nodesIdVec,matIdVol,*gmesh);
                break;
            case EQuadrilateral:
                quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticQuad>(elId,nodesIdVec,matIdVol,*gmesh);
                nodesIdVec.Resize(nNodesVolEl);
                for(int i = 0; i < nodesIdVec.size(); i++ ) nodesIdVec[i] = i + 2 * nNodesVolEl;
                blendEl = new TPZGeoElRefLess<pzgeom::TPZGeoBlend<TGeo>>(elId,nodesIdVec,matIdVol,*gmesh);
                break;
            default:
                DebugStop();
                break; //the element should not be deleted
        }

        if (printGMesh) {
            std::string meshFileName = prefix + "gmesh_" + TGeo::TypeName();
            const size_t strlen = meshFileName.length();
            meshFileName.replace(strlen, 4, ".txt");
            std::ofstream outTXT(meshFileName.c_str());
            gmesh->Print(outTXT);
            outTXT.close();
        }
        ///////CREATE MAX DIM - 1 ELEMENTS FOR BLENDING
        {
            TPZStack<int> containedSides;
            TGeo::LowerDimensionSides(TGeo::NSides - 1, containedSides, TGeo::Dimension - 1);
            for (int64_t subSideIndex = 0;
                 subSideIndex < containedSides.NElements(); subSideIndex++){
                const int64_t subSide = containedSides[subSideIndex];
                const int64_t nNodesSide = TGeo::NSideNodes(subSide);
                nodesIdVec.Resize(nNodesSide+1);//TODO: depends on the type of the side
                for (int64_t node = 0; node < nNodesSide; node++){
                    const int64_t nodeId = TGeo::SideNodeLocId(subSide,node);
                    nodesIdVec[node] = blendEl->Node(nodeId).Id();
//                std::cout<<nodesIdVec[node]<<"\t";
                }
                nodesIdVec[nNodesSide] = nodesIdArcsVec[subSideIndex];//TODO: depends on the type of the side
//            std::cout<<nodesIdVec[nNodes]<<std::endl;
                new TPZGeoElRefLess<pzgeom::TPZQuadraticLine>(nodesIdVec, matIdSphere, *gmesh);//TODO: depends on the type of the side
            }
        }
        gmesh->BuildConnectivity();


        if (printGMesh) {
            std::string meshFileName = prefix + "gmesh_" + TGeo::TypeName();
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
        //X COMPARE
        {
            TPZManVector<REAL,3> xi;
            REAL notUsedHere = -1;
            const REAL tol = 1e-8;
            const int pOrder = 10;
            uint64_t errors = 0;
            auto intRule = blendEl->CreateSideIntegrationRule(blendEl->NSides()-1, pOrder);
            xi.Resize(dim,0);
            for(int iPt = 0; iPt < intRule->NPoints(); iPt++){
                bool hasAnErrorOccurred = false;
                intRule->Point(iPt,xi,notUsedHere);
                TPZManVector<REAL,3> xBlend(3);
                blendEl->X(xi, xBlend);
                TPZManVector<REAL,3> xQuad(3);
                quadraticEl->X(xi, xQuad);
                std::ostringstream xBlendM, xQuadM;
                xBlendM<<"x_blend:"<<std::endl;
                xQuadM<<"x_quad:"<<std::endl;
                const auto VAL_WIDTH = 15;
                REAL diff = 0;
                for(int i = 0; i < xBlend.size(); i++){
                    xBlendM<<std::setw(VAL_WIDTH) << std::right<<xBlend[i] - coordsOffset[i]<<"\t";
                    xQuadM <<std::setw(VAL_WIDTH) << std::right<<xQuad[i] <<"\t";
                    diff += (xBlend[i] - coordsOffset[i] - xQuad[i])*(xBlend[i] - coordsOffset[i] - xQuad[i]);
                }
                xBlendM<<std::endl;
                xQuadM<<std::endl;
                if(diff > tol){
                    if(!hasAnErrorOccurred) errors++;
                    hasAnErrorOccurred = true;
                }
                if(hasAnErrorOccurred){
                    std::cout<<std::flush;
                    std::cout<<xBlendM.str()<<std::endl;
                    std::cout<<xQuadM.str()<<std::endl;
                    std::cout<<"diff :"<<diff<<std::endl;
                }
            }
            std::cout<<"============================"<<std::endl;
            std::cout<<"Element: "<<MElementType_Name(TGeo::Type());
            std::cout<<"\tNumber of points: "<<intRule->NPoints()<<"\tErrors: "<<errors<<std::endl;
        }
    }

    void CreateGeoMesh2D(TPZGeoMesh *&gmesh, int nDiv, bool printGMesh, std::string prefix) {
        gmesh = new TPZGeoMesh();

        TPZVec<REAL> coord(3, 0.);
        for (int64_t i = 0; i < 4; i++) {
            coord[0] = i/2 == i % 2 ? -1 : 1;
            coord[1] = -1 + 2 * (i/2);
//            std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
            const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newindex].Initialize(coord,*gmesh);
        }//quad nodes

        TPZManVector<int64_t,4> nodesIdArcsVec(4);
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
            const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newindex].Initialize(coord,*gmesh);
            nodesIdArcsVec[i] = newindex;
        }//quad nodes

        int matIdQuad = 1, matIdArc = 2;
        TPZManVector<int64_t,4> nodesIdVec(1);
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
//                std::cout<<nodesIdVec[node]<<"\t";
            }
            nodesIdVec[nNodes] = nodesIdArcsVec[edge-firstEdge];
//            std::cout<<nodesIdVec[nNodes]<<std::endl;
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

        //Create GeoBlend elements
        {
            TPZGeoMesh *newgmesh = new TPZGeoMesh();
            for (int64_t i = 0; i < gmesh->NNodes(); i++) {
                coord[0] = gmesh->NodeVec()[i].Coord(0);
                coord[1] = gmesh->NodeVec()[i].Coord(1);
                coord[2] = gmesh->NodeVec()[i].Coord(2);
                const int64_t newindex = newgmesh->NodeVec().AllocateNewElement();
                newgmesh->NodeVec()[newindex].Initialize(coord,*newgmesh);
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
                    case ETriangle:
                        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(elId,nodesIdVec,matId,*newgmesh);
                        break;
                    case EQuadrilateral:
                        if(geo->HasSubElement()) continue;//the original cube should not be copied
                        new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>>(elId,nodesIdVec,matId,*newgmesh);
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
            std::string meshFileName = prefix + "gmesh2D_partial";
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

        #ifdef _AUTODIFF
        {
            TPZManVector<REAL,3> xiReal;
            TPZManVector<Fad<REAL>,3> xiFad;
            REAL weight = -1;//useless
            const REAL tol = 1e-8;
            const int pOrder = 3;
            const int nel = gmesh->NElements();
            for (int iel = 0; iel < nel; iel++) {
                TPZGeoEl *geo = gmesh->ElementVec()[iel];
                if (geo && !geo->HasSubElement()) {
                    uint64_t errors = 0;
                    auto intRule = geo->CreateSideIntegrationRule(geo->NSides()-1, pOrder);
                    xiReal.Resize(geo->Dimension(),0);
                    xiFad.Resize(geo->Dimension(),0);
                    for(int iPt = 0; iPt < intRule->NPoints(); iPt++){
                        bool hasAnErrorOccurred = false;
                        intRule->Point(iPt,xiReal,weight);
                        for(int x = 0; x < geo->Dimension(); x++){
                            xiFad[x] = Fad<REAL>(geo->Dimension(),x,xiReal[x]);
                        }
//                        Fad<REAL> func;
//                        func = 9*xiFad[0];
//                        for(int i =0; i < func.size(); i++){
//                            std::cout<<"dx["<<i<<"]:\t"<<func.dx(i)<<std::endl;
//                        }
                        TPZManVector<Fad<REAL>,3> xFad(3);
                        geo->X(xiFad, xFad);
                        TPZFNMatrix<9,REAL> gradXreal(3,3,0);
                        geo->GradX(xiReal, gradXreal);
                        std::ostringstream xFadM, gradXM;
                        xFadM<<"xFad:"<<std::endl;
                        gradXM<<"grad x:"<<std::endl;
                        const auto VAL_WIDTH = 10;
                        for(int i = 0; i < gradXreal.Rows(); i++){
                            for(int j = 0; j < gradXreal.Cols(); j++){
                                xFadM<<std::setw(VAL_WIDTH) << std::right<<xFad[i].dx(j)<<"\t";
                                gradXM<<std::setw(VAL_WIDTH) << std::right<<gradXreal(i,j)<<"\t";
                                const REAL diff = (xFad[i].dx(j)-gradXreal(i,j))*(xFad[i].dx(j)-gradXreal(i,j));
                                if(diff > tol){
                                    if(!hasAnErrorOccurred) errors++;
                                    hasAnErrorOccurred = true;
                                }
                            }
                            xFadM<<std::endl;
                            gradXM<<std::endl;
                        }
                        if(hasAnErrorOccurred){
                            std::cout<<std::flush;
                            std::cout<<xFadM.str()<<std::endl;
                            std::cout<<gradXM.str()<<std::endl;
                        }
                    }
                    std::cout<<"============================"<<std::endl;
                    std::cout<<"Element: "<<geo->Id()<<"\tType: "<<MElementType_Name(geo->Type());
                    std::cout<<"\tIs blend? : "<<geo->IsGeoBlendEl()<<std::endl;
                    std::cout<<"\tNumber of points: "<<intRule->NPoints()<<"\tErrors: "<<errors<<std::endl;
                }
            }
        }
        #endif


        {
            TPZVec<TPZGeoEl *> sons;
            std::vector<std::string> loading = {"-","/","|","\\"};
            for (int iDiv = 0; iDiv < nDiv; iDiv++) {
                std::cout<<"Performing "<<iDiv+1<<" ref step out of " << nDiv<<std::endl;
                const int nel = gmesh->NElements();
                for (int iel = 0; iel < nel; iel++) {
                    std::cout<<"\b"<<loading[iel%4]<<std::flush;
                    TPZGeoEl *geo = gmesh->ElementVec()[iel];
                    if (geo && !geo->HasSubElement()) {
                        geo->Divide(sons);
                    }
                }
                std::cout<<"\b"<<std::endl;
            }
        }

        if (printGMesh) {
            std::string meshFileName = prefix + "gmesh2D";
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
            const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[newindex].Initialize(coord,*gmesh);
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
            for (int64_t i = 0; i < gmesh->NNodes(); i++) {
                coord[0] = gmesh->NodeVec()[i].Coord(0);
                coord[1] = gmesh->NodeVec()[i].Coord(1);
                coord[2] = gmesh->NodeVec()[i].Coord(2);
                const int64_t newindex = newgmesh->NodeVec().AllocateNewElement();
                newgmesh->NodeVec()[newindex].Initialize(coord,*newgmesh);
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
            std::string meshFileName = prefix + "gmesh3D_partial";
            std::cout<<"Printing "<<meshFileName<<".vtk and .txt"<<std::endl;
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
#ifdef _AUTODIFF
        {
            TPZManVector<REAL,3> xiReal;
            TPZManVector<Fad<REAL>,3> xiFad;
            REAL weight = -1;//useless
            const REAL tol = 1e-8;
            const int pOrder = 3;
            const int nel = gmesh->NElements();
            for (int iel = 0; iel < nel; iel++) {
                TPZGeoEl *geo = gmesh->ElementVec()[iel];
                if (geo && !geo->HasSubElement()) {
                    uint64_t errors = 0;
                    auto intRule = geo->CreateSideIntegrationRule(geo->NSides()-1, pOrder);
                    xiReal.Resize(geo->Dimension(),0);
                    xiFad.Resize(geo->Dimension(),0);
                    for(int iPt = 0; iPt < intRule->NPoints(); iPt++){
                        bool hasAnErrorOccurred = false;
                        intRule->Point(iPt,xiReal,weight);
                        for(int x = 0; x < geo->Dimension(); x++){
                            xiFad[x] = Fad<REAL>(geo->Dimension(),x,xiReal[x]);
                        }
//                        Fad<REAL> func;
//                        func = 9*xiFad[0];
//                        for(int i =0; i < func.size(); i++){
//                            std::cout<<"dx["<<i<<"]:\t"<<func.dx(i)<<std::endl;
//                        }
                        TPZManVector<Fad<REAL>,3> xFad(3);
                        geo->X(xiFad, xFad);
                        TPZFNMatrix<9,REAL> gradXreal(3,3,0);
                        geo->GradX(xiReal, gradXreal);
                        std::ostringstream xFadM, gradXM;
                        xFadM<<"xFad:"<<std::endl;
                        gradXM<<"grad x:"<<std::endl;
                        const auto VAL_WIDTH = 10;
                        for(int i = 0; i < gradXreal.Rows(); i++){
                            for(int j = 0; j < gradXreal.Cols(); j++){
                                xFadM<<std::setw(VAL_WIDTH) << std::right<<xFad[i].dx(j)<<"\t";
                                gradXM<<std::setw(VAL_WIDTH) << std::right<<gradXreal(i,j)<<"\t";
                                const REAL diff = (xFad[i].dx(j)-gradXreal(i,j))*(xFad[i].dx(j)-gradXreal(i,j));
                                if(diff > tol){
                                    if(!hasAnErrorOccurred) errors++;
                                    hasAnErrorOccurred = true;
                                }
                            }
                            xFadM<<std::endl;
                            gradXM<<std::endl;
                        }
                        if(hasAnErrorOccurred){
                            std::cout<<std::flush;
                            std::cout<<xFadM.str()<<std::endl;
                            std::cout<<gradXM.str()<<std::endl;
                        }
                    }
                    std::cout<<"============================"<<std::endl;
                    std::cout<<"Element: "<<geo->Id()<<"\tType: "<<MElementType_Name(geo->Type());
                    std::cout<<"\tIs blend? : "<<geo->IsGeoBlendEl()<<std::endl;
                    std::cout<<"\tNumber of points: "<<intRule->NPoints()<<"\tErrors: "<<errors<<std::endl;
                }
            }
        }
#endif
        {
            TPZVec<TPZGeoEl *> sons;
            std::vector<std::string> loading = {"-","/","|","\\"};
            for (int iDiv = 0; iDiv < nDiv; iDiv++) {
                std::cout<<"Performing "<<iDiv+1<<" ref step out of " << nDiv<<std::endl;
                const int nel = gmesh->NElements();
                for (int iel = 0; iel < nel; iel++) {
                    std::cout<<"\b"<<loading[iel%4]<<std::flush;
                    TPZGeoEl *geo = gmesh->ElementVec()[iel];
                    if (geo && !geo->HasSubElement()) {
                        geo->Divide(sons);
                    }
                }
                std::cout<<"\b"<<std::endl;
            }
        }

        if (printGMesh) {
            std::string meshFileName = prefix + "gmesh3D";
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