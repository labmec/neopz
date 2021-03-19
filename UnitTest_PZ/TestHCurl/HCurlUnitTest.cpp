/**
 * @file HCurlUnitTest.cpp
 * @brief Define a Unit Test using Boost for validation of HCurl-conforming approximation spaces
 *
 */
#include <iostream>

#include "pzlog.h"
#include "TPZTopologyUtils.h"
#include "pzgmesh.h"

#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"

#include "pzgeoelrefless.h"
#ifdef PZ_LOG
static PZLogger logger("pz.mesh.testhcurl");
#endif

#include "fad.h"

#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"
#include "pzcmesh.h"
#include "TPZMatHCurlProjection.cpp"
#include "pzanalysis.h"
#include "pzintel.h"
#include "TPZCompElHCurl.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "TPZVTKGeoMesh.h"
// Using Unit Test of the Boost Library
#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz #define BOOST_TEST_MAIN pz hcurl_tests tests

#include "boost/test/unit_test.hpp"
#include "boost/test/tools/output_test_stream.hpp"


#define NOISY_HCURL //outputs useful debug info
#define NOISY_HCURL_VTK


struct SuiteInitializer{
    SuiteInitializer(){
        boost::unit_test::unit_test_log.set_threshold_level( boost::unit_test::log_warnings );
    }
};
BOOST_FIXTURE_TEST_SUITE(hcurl_tests,SuiteInitializer)
    namespace hcurltest{
        constexpr REAL tol = 1e-10;
        /**
         * This unit test aims to verify the vectors used to build the HCurl approximation space.
         * It was rewritten based on TestFunctionTracesUniformMesh, therefore it could be simplified.
         * @param type element type (triangle, quadrilateral, etc.
         * @param dim dimension of the element
         */
        void TestVectorTracesUniformMesh(MMeshType type, const int dim);
        /**
         * This unit test aims to verify the trace compatibility of the HCurl approximation space in a UNIFORM mesh.
         * @param type element type (triangle, quadrilateral, etc.
         * @param dim dimension of the element
         */
        void TestFunctionTracesUniformMesh(MMeshType type, const int pOrder, const int dim);

        template <class TGEOM>
        void PrintShapeFunctions(const int pOrder);

        namespace auxiliaryfuncs{
            void VectorProduct(const TPZVec<REAL> &, const TPZVec<REAL> &, TPZVec<REAL> &);

            TPZGeoMesh *CreateGeoMesh2D(int nelx, int nely, MMeshType meshType, TPZVec<int> &matIds);

            TPZGeoMesh *CreateGeoMesh3D(int nelx, int nely, int nelz, MMeshType meshType, TPZVec<int> &matIds);
        }
    }


    BOOST_AUTO_TEST_CASE(hcurl_topology_tests) {
        constexpr int dim2D{2},dim3D{3};
        hcurltest::TestVectorTracesUniformMesh(MMeshType::ETriangular,dim2D);
        hcurltest::TestVectorTracesUniformMesh(MMeshType::EQuadrilateral,dim2D);
        hcurltest::TestVectorTracesUniformMesh(MMeshType::ETetrahedral,dim3D);
        hcurltest::TestVectorTracesUniformMesh(MMeshType::EHexahedral,dim3D);
        hcurltest::TestVectorTracesUniformMesh(MMeshType::EPrismatic,dim3D);
    }

    BOOST_AUTO_TEST_CASE(hcurl_mesh_tests) {
//        hcurltest::PrintShapeFunctions<pzgeom::TPZGeoTetrahedra>(4);
        constexpr int dim2D{2},dim3D{3};
        for(auto k = 1; k <= 5; k++) hcurltest::TestFunctionTracesUniformMesh(MMeshType::ETriangular, k,dim2D);
        for(auto k = 1; k <= 5; k++) hcurltest::TestFunctionTracesUniformMesh(MMeshType::EQuadrilateral, k,dim2D);
        for(auto k = 1; k <= 5; k++) hcurltest::TestFunctionTracesUniformMesh(MMeshType::ETetrahedral, k,dim3D);
        for(auto k = 1; k <= 5; k++) hcurltest::TestFunctionTracesUniformMesh(MMeshType::EHexahedral, k,dim3D);
        for(auto k = 1; k <= 5; k++) hcurltest::TestFunctionTracesUniformMesh(MMeshType::EPrismatic, k,dim3D);

    }

    namespace hcurltest{

        void TestVectorTracesUniformMesh(MMeshType type, const int dim){
            static std::string testName = __PRETTY_FUNCTION__;
            std::cout << testName << std::endl;
            std::cout<<"\t"<<type<<std::endl;
            constexpr int pOrder{4};//testing all vectors
            constexpr int ndiv{1};
            TPZManVector<int,2> matIds(2,-1);
            TPZGeoMesh *gmesh = [&]() -> TPZGeoMesh* {
                switch(dim){
                    case 2:
                        return auxiliaryfuncs::CreateGeoMesh2D(ndiv,ndiv,type,matIds);
                    case 3:
                        return auxiliaryfuncs::CreateGeoMesh3D(ndiv,ndiv,ndiv,type,matIds);
                    default:
                        DebugStop();
                }
                return nullptr;
            }();

#ifdef NOISY_HCURL
            std::ofstream outTXT("gmesh_"+MMeshType_Name(type)+"_ndiv_"+std::to_string(ndiv)+ ".txt");
            gmesh->Print(outTXT);
            outTXT.close();
#endif
#ifdef NOISY_HCURL_VTK
            std::ofstream outVTK("gmesh_"+MMeshType_Name(type)+"_ndiv_"+std::to_string(ndiv)+ ".vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK);
            outVTK.close();
#endif
            auto cmesh = new TPZCompMesh(gmesh);
            cmesh->SetDefaultOrder(pOrder);//with p=4 you will have all the vectors
            cmesh->SetDimModel(dim);

            auto mat = new TPZMatHCurlProjection(dim,matIds[0],1);
            cmesh->InsertMaterialObject(mat);
            auto bcond = mat->CreateBC(mat,matIds[1],0,TPZFNMatrix<1,STATE>(1,1,0),TPZFNMatrix<1,STATE>(1,1,0));
            cmesh->InsertMaterialObject(bcond);
            cmesh->SetAllCreateFunctionsHCurl();
            cmesh->AutoBuild();
            cmesh->AdjustBoundaryElements();
            cmesh->CleanUpUnconnectedNodes();

#ifdef NOISY_HCURL
            {
                std::ofstream outTXT("cmesh_"+MMeshType_Name(type)+"_ndiv_"+std::to_string(ndiv)+ ".txt");
                cmesh->Print(outTXT);
                outTXT.close();
            }
#endif

            ///this lambda expression will help calculating the tangential traces
            auto CheckTracesFunc = [](REAL &diffTrace, const TPZVec<REAL> &elVec, const TPZVec<REAL>&neighVec,
                                      const TPZVec<REAL> &vec, const int sideDim, TPZVec<REAL> &elTrace, TPZVec<REAL> &neighTrace){
                REAL averageTrace{0};
                diffTrace = 0;
                switch(sideDim){
                    case 1:{
                        elTrace.Resize(1,0);
                        neighTrace.Resize(1,0);
                        for (auto x = 0; x < 3; x++) elTrace[0] += elVec[x]*vec[x];
                        for (auto x = 0; x < 3; x++) neighTrace[0] += neighVec[x]*vec[x];
                        diffTrace = std::abs(elTrace[0] - neighTrace[0]);
                        averageTrace = (sqrt(elTrace[0]*elTrace[0])+sqrt(neighTrace[0]*neighTrace[0]))/2;
                    }
                        break;
                    case 2:{
                        TPZManVector<REAL,3> temp(3,0);
                        elTrace.Resize(3,0);
                        neighTrace.Resize(3,0);
                        auxiliaryfuncs::VectorProduct(vec,elVec,temp);
                        auxiliaryfuncs::VectorProduct(vec,temp,elTrace);

                        auxiliaryfuncs::VectorProduct(vec,neighVec,temp);
                        auxiliaryfuncs::VectorProduct(vec,temp,neighTrace);
                        for(auto x = 0; x < 3; x++) diffTrace+= (neighTrace[x]-elTrace[x])*(neighTrace[x]-elTrace[x]);
                        diffTrace= sqrt(diffTrace);
                        REAL firstTrace = 0, secondTrace = 0;
                        for(auto x = 0; x < 3; x++) firstTrace+= elTrace[x]*elTrace[x];
                        for(auto x = 0; x < 3; x++) secondTrace+= neighTrace[x]*neighTrace[x];
                        firstTrace = sqrt(firstTrace);
                        secondTrace = sqrt(secondTrace);
                        averageTrace =(firstTrace+secondTrace)/2;
                    }
                        break;
                    default:
                        return false;
                }
                if(averageTrace < tol)  averageTrace = 1;

                return  diffTrace/averageTrace < tol;
            };

            MElementType elType = [&](){
               switch(type){
                   case MMeshType::ETriangular: return ETriangle;
                   case MMeshType::EQuadrilateral: return EQuadrilateral;
                   case MMeshType::ETetrahedral: return ETetraedro;
                   case MMeshType::EHexahedral: return ECube;
                   case MMeshType::EPrismatic: return EPrisma;
                   case MMeshType::EPyramidal: return EPiramide;
                   default: return ENoType;
               }
            }();
            for(auto dummyCel : cmesh->ElementVec()){
                const auto cel = dynamic_cast<TPZInterpolatedElement *>(dummyCel);
                const auto gel = cel->Reference();
                //skips boundary els
                if(!cel || cel->Reference()->Type() != elType) continue;

                TPZMaterialData elData;
                cel->InitMaterialData(elData);
                const int elNNodes = MElementType_NNodes(elType);
                for (auto iCon = 0; iCon <cel->NConnects(); iCon++) {
                    auto &con = cel->Connect(iCon);

                    if(con.NElConnected() < 2) continue;

                    const int iSide = iCon + elNNodes;

                    TPZTransform<> elTransform(gel->SideToSideTransform(iSide, gel->NSides() - 1));
                    TPZGeoElSide gelSide(gel, iSide);
                    const auto sideDim = gelSide.Dimension();
                    //the following vector will be the edge tg vector if 2D, the normal vector if 3D
                    TPZManVector<REAL,3> vec(3,0);
                    switch(sideDim){
                        case 1:{
                            TPZManVector<int, 2> edgeNodes(2, 0);
                            for (auto i = 0; i < 2; i++) edgeNodes[i] = gel->SideNodeIndex(iSide, i);
                            const REAL sign = edgeNodes[0] < edgeNodes[1] ? 1 : -1;
                            REAL edgeLength = 0;
                            TPZVec<REAL> p0(3), p1(3);
                            gmesh->NodeVec()[edgeNodes[0]].GetCoordinates(p0);
                            gmesh->NodeVec()[edgeNodes[1]].GetCoordinates(p1);
                            edgeLength =
                                    std::sqrt((p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
                                              (p0[2] - p1[2]) * (p0[2] - p1[2]));
                            for (auto x = 0; x < 3; x++) vec[x] = sign * (p1[x] - p0[x]) / edgeLength;
                        }
                            break;
                        case 2:{
                            TPZManVector<REAL,3> xCenter(2,0);
                            switch(gelSide.NSideNodes()){
                                case 3:
                                    pztopology::TPZTriangle::CenterPoint(pztopology::TPZTriangle::NSides-1,xCenter);
                                    break;
                                case 4:
                                    pztopology::TPZQuadrilateral::CenterPoint(pztopology::TPZQuadrilateral::NSides-1,xCenter);
                                    break;
                                default:
                                    DebugStop();
                            }
                            gelSide.Normal(xCenter,vec);
                        }
                            break;
                        default:
                            DebugStop();
                    }

                    TPZManVector<REAL,3> centerPt(sideDim,0);
                    gelSide.CenterPoint(centerPt);
                    TPZManVector<REAL,3> ptEl(gel->Dimension(),0);
                    elTransform.Apply(centerPt,ptEl);
                    cel->ComputeRequiredData(elData, ptEl);
                    TPZFNMatrix<60,REAL> elDeformedDirections(elData.fDeformedDirections);

                    //gather all the sides contained in the closure of iSide
                    TPZStack<int> smallSides;
                    gel->LowerDimensionSides(iSide,smallSides);
                    smallSides.Push(iSide);//include the side itself

                    TPZGeoElSide neighGelSide = gelSide.Neighbour();
                    while(neighGelSide != gelSide) {
                        const auto neighCel = dynamic_cast<TPZInterpolatedElement *> (neighGelSide.Element()->Reference());
                        if (!neighCel) {
                            neighGelSide = neighGelSide.Neighbour();
                            continue;
                        }
                        const auto neighGel = neighCel->Reference();
                        const int neighNNodes = MElementType_NNodes(neighGel->Type());
                        const auto neighSide = neighGelSide.Side();
                        const auto neighDim = neighGelSide.Element()->Dimension();
                        TPZTransform<> neighTransform(neighCel->Reference()->SideToSideTransform(neighSide,
                                                                                                 neighGel->NSides() -1));
                        TPZTransform<> localTransf(sideDim);
                        gelSide.SideTransform3(neighGelSide,localTransf);
                        TPZMaterialData neighData;
                        neighCel->InitMaterialData(neighData);

                        TPZManVector <REAL,3> ptsN(neighGelSide.Dimension(),0),ptNeigh(neighDim);
                        localTransf.Apply(centerPt,ptsN);
                        neighTransform.Apply(ptsN,ptNeigh);
                        neighCel->ComputeRequiredData(neighData, ptNeigh);

                        TPZFNMatrix<60,REAL> neighDeformedDirections(neighData.fDeformedDirections);
                        //since this test was based on the basis functions tests, the following std::set
                        //is a lazy way to avoid multiple testing of the same vectors
                        std::set<int64_t> testedVectors;
                        for(auto subSide : smallSides){
                            TPZGeoElSide gelSubSide(gel, subSide);
                            if(gel->SideDimension(subSide) < 1) continue;
                            const int subConnect = subSide - elNNodes;
                            const int neighSubSide = [&](){
                                TPZGeoElSide neighGelSubSide = gelSubSide.Neighbour();
                                while(neighGelSubSide.Element() != neighGelSide.Element()) {
                                    neighGelSubSide = neighGelSubSide.Neighbour();
                                    if(neighGelSubSide.Element() == gel){
                                        DebugStop();
                                    }
                                }
                                return neighGelSubSide.Side();
                            }();
                            const int nShapes = cel->NConnectShapeF(subConnect,pOrder);

                            const int firstElShape = [&](){
                                int firstElShapeTemp = 0;
                                for(auto jCon = 0; jCon < subConnect; jCon++){
                                    firstElShapeTemp += cel->NConnectShapeF(jCon,cel->EffectiveSideOrder(jCon+elNNodes));
                                }
                                return firstElShapeTemp;
                            }();
                    

                            const int firstNeighShape = [&](){
                                int firstNeighShapeTemp = 0;
                                const int neighCon = neighSubSide - neighGelSide.Element()->NNodes();
                                for(auto jCon = 0; jCon < neighCon; jCon++){
                                    firstNeighShapeTemp += neighCel->NConnectShapeF(jCon,neighCel->EffectiveSideOrder(jCon+neighNNodes));
                                }
                                return firstNeighShapeTemp;
                            }();

                            TPZManVector<REAL,3> elVec(3,0), neighVec(3,0);
                            for(auto iShape = 0; iShape < nShapes; iShape ++){
                                const int elPhiIndex = firstElShape+iShape;
                                const int neighPhiIndex = firstNeighShape+iShape;
                                const int elVecIndex = elData.fVecShapeIndex[elPhiIndex].first;
                                const int neighVecIndex = neighData.fVecShapeIndex[neighPhiIndex].first;
                                if (testedVectors.find(elVecIndex) == testedVectors.end()) {
                                    testedVectors.insert(elVecIndex);
                                }else{
                                    continue;
                                }
                                for (auto x = 0; x < 3; x++) elVec[x] = elDeformedDirections(x,elVecIndex);
                                for (auto x = 0; x < 3; x++) neighVec[x] = neighDeformedDirections(x,neighVecIndex);
                                REAL diffTrace{0};
                                TPZManVector<REAL,3> elTrace,neighTrace;
                                const bool checkTraces = CheckTracesFunc(diffTrace,elVec,neighVec,vec,sideDim,elTrace,neighTrace);

                                {
                                    std::ostringstream traceMsg;
                                    if(!checkTraces){
                                        traceMsg <<"el index: "<<gel->Index()<<std::endl;
                                        traceMsg <<"neigh index: "<<neighGel->Index()<<std::endl;
                                        traceMsg <<"el vec index: "<<elVecIndex<<std::endl;
                                        traceMsg <<"neigh vec index: "<<neighVecIndex<<std::endl;
                                        traceMsg <<"el vec (master el): ";
                                        for (auto x = 0; x < gel->Dimension(); x++) traceMsg << elData.fMasterDirections(x,elVecIndex)<<"\t";
                                        traceMsg <<"\n";
                                        traceMsg<<"neigh vec (master el): ";
                                        for (auto x = 0; x < neighDim; x++) traceMsg << neighData.fMasterDirections(x,neighVecIndex)<<"\t";
                                        traceMsg <<"\n";
                                        traceMsg <<"el vec (deformed): ";
                                        for (auto x = 0; x < 3; x++) traceMsg << elVec[x]<<"\t";
                                        traceMsg<<"\nneigh vec (deformed): ";
                                        for (auto x = 0; x < 3; x++) traceMsg << neighVec[x]<<"\t";
                                        traceMsg <<"\n";
                                        traceMsg <<"el trace: ";
                                        for (auto x = 0; x < elTrace.size(); x++) traceMsg << elTrace[x]<<"\t";
                                        traceMsg<<"\nneigh trace: ";
                                        for (auto x = 0; x < neighTrace.size(); x++) traceMsg << neighTrace[x]<<"\t";
                                        traceMsg <<"\n";
                                        if(sideDim == 2) traceMsg<<"normal vector: ";
                                        else traceMsg<<"tangent vector: ";
                                        for (auto x = 0; x < vec.size(); x++) traceMsg << vec[x]<<"\t";
                                        traceMsg <<"\n";
                                        traceMsg<<"el jacobian: ";
                                        elData.jacobian.Print(traceMsg);
                                        traceMsg<<"el axes: ";
                                        elData.axes.Print(traceMsg);
                                        traceMsg<<"neigh jacobian: ";
                                        neighData.jacobian.Print(traceMsg);
                                        traceMsg<<"neigh axes: ";
                                        neighData.axes.Print(traceMsg);
                                    }
                                    BOOST_CHECK_MESSAGE(checkTraces,"\n"+testName+" failed"+
                                                                    "\ntopology: "+MMeshType_Name(type)+"\n"+
                                                                    "side: "+std::to_string(iSide)+"\n"+
                                                                    "diff traces: "+std::to_string(diffTrace)+"\n"+
                                                                    traceMsg.str()
                                    );
                                }
                            }
                        }
                        neighGelSide = neighGelSide.Neighbour();
                    }
                }
            }
            delete cmesh;
            delete gmesh;
        }//hcurltest::TestVectorTracesUniformMesh

        void TestFunctionTracesUniformMesh(MMeshType type, const int pOrder, const int dim){
            static std::string testName = __PRETTY_FUNCTION__;
            std::cout << testName << std::endl;
            std::cout<<"\t"<<MMeshType_Name(type)<<" p = "<<pOrder<<std::endl;
            constexpr int ndiv{1};
            TPZManVector<int,2> matIds(2,-1);
            TPZGeoMesh *gmesh = [&]() -> TPZGeoMesh* {
                switch(dim){
                    case 2:
                        return auxiliaryfuncs::CreateGeoMesh2D(ndiv,ndiv,type,matIds);
                    case 3:
                        return auxiliaryfuncs::CreateGeoMesh3D(ndiv,ndiv,ndiv,type,matIds);
                    default:
                        DebugStop();
                }
                return nullptr;
            }();
//            {
//                TPZVec<TPZGeoEl*> sons;
//                gmesh->Element(0)->Divide(sons);
//            }
#ifdef NOISY_HCURL
            std::ofstream outTXT("gmesh_"+MMeshType_Name(type)+"_ndiv_"+std::to_string(ndiv)+ ".txt");
            gmesh->Print(outTXT);
            outTXT.close();
#endif
#ifdef NOISY_HCURL_VTK
            std::ofstream outVTK("gmesh_"+MMeshType_Name(type)+"_ndiv_"+std::to_string(ndiv)+ ".vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK);
            outVTK.close();
#endif
            auto cmesh = new TPZCompMesh(gmesh);
            cmesh->SetDefaultOrder(pOrder);
            cmesh->SetDimModel(dim);

            auto mat = new TPZMatHCurlProjection(dim,matIds[0],1);
            cmesh->InsertMaterialObject(mat);
            auto bcond = mat->CreateBC(mat,matIds[1],0,TPZFNMatrix<1,STATE>(1,1,0),TPZFNMatrix<1,STATE>(1,1,0));
            cmesh->InsertMaterialObject(bcond);
            cmesh->SetAllCreateFunctionsHCurl();
            cmesh->AutoBuild();
            cmesh->AdjustBoundaryElements();
            cmesh->CleanUpUnconnectedNodes();

#ifdef NOISY_HCURL
            {
                std::ofstream outTXT("cmesh_"+MMeshType_Name(type)+"_ndiv_"+std::to_string(ndiv)+ ".txt");
                cmesh->Print(outTXT);
                outTXT.close();
            }
#endif

            ///this lambda expression will help calculating the tangential traces
            auto CheckTracesFunc = [](REAL &diffTrace, const TPZVec<REAL> &elShapeFunc, const TPZVec<REAL>&neighShapeFunc,
                                      const TPZVec<REAL> &vec, const int sideDim, TPZVec<REAL> &elTrace, TPZVec<REAL> &neighTrace){
                REAL averageTrace{0};
                diffTrace = 0;
                switch(sideDim){
                    case 1:{
                        elTrace.Resize(1,0);
                        neighTrace.Resize(1,0);
                        for (auto x = 0; x < 3; x++) elTrace[0] += elShapeFunc[x]*vec[x];
                        for (auto x = 0; x < 3; x++) neighTrace[0] += neighShapeFunc[x]*vec[x];
                        diffTrace = std::abs(elTrace[0] - neighTrace[0]);
                        averageTrace = (sqrt(elTrace[0]*elTrace[0])+sqrt(neighTrace[0]*neighTrace[0]))/2;
                    }
                        break;
                    case 2:{
                        TPZManVector<REAL,3> temp(3,0);
                        elTrace.Resize(3,0);
                        neighTrace.Resize(3,0);
                        auxiliaryfuncs::VectorProduct(vec,elShapeFunc,temp);
                        auxiliaryfuncs::VectorProduct(temp,vec,elTrace);

                        auxiliaryfuncs::VectorProduct(vec,neighShapeFunc,temp);
                        auxiliaryfuncs::VectorProduct(temp,vec,neighTrace);
                        for(auto x = 0; x < 3; x++) diffTrace+= (neighTrace[x]-elTrace[x])*(neighTrace[x]-elTrace[x]);
                        diffTrace= sqrt(diffTrace);
                        REAL firstTrace = 0, secondTrace = 0;
                        for(auto x = 0; x < 3; x++) firstTrace+= elTrace[x]*elTrace[x];
                        for(auto x = 0; x < 3; x++) secondTrace+= neighTrace[x]*neighTrace[x];
                        firstTrace = sqrt(firstTrace);
                        secondTrace = sqrt(secondTrace);
                        averageTrace =(firstTrace+secondTrace)/2;
                    }
                        break;
                    default:
                        return false;
                }
                if(averageTrace < tol)  averageTrace = 1;

                return  diffTrace/averageTrace < tol;
            };

            MElementType elType = [&](){
                switch(type){
                    case MMeshType::ETriangular: return ETriangle;
                    case MMeshType::EQuadrilateral: return EQuadrilateral;
                    case MMeshType::ETetrahedral: return ETetraedro;
                    case MMeshType::EHexahedral: return ECube;
                    case MMeshType::EPrismatic: return EPrisma;
                    case MMeshType::EPyramidal: return EPiramide;
                    default: return ENoType;
                }
            }();
            for(auto dummyCel : cmesh->ElementVec()){
                const auto cel = dynamic_cast<TPZInterpolatedElement *>(dummyCel);
                if(!cel) continue;
                const auto gel = cel->Reference();
                //skips boundary els
                if(!cel || cel->Reference()->Type() != elType) continue;

                TPZMaterialData elData;
                cel->InitMaterialData(elData);
                const int nState = cel->Material()->NStateVariables();
                const int elNNodes = MElementType_NNodes(elType);
                for (auto iCon = 0; iCon <cel->NConnects(); iCon++) {
                    auto &con = cel->Connect(iCon);
                    const int nShape = con.NShape();
                    {
                        const bool check = con.NDof( *cmesh) == nShape * nState;
                        BOOST_CHECK_MESSAGE(check,"\n"+testName+" failed"+
                                                  "\ntopology: "+MMeshType_Name(type)+"\n"
                        );
                    }

                    if(con.NElConnected() < 2) continue;

                    const int iSide = iCon + elNNodes;
                    const int pOrderIntRule = cel->EffectiveSideOrder(iSide)*2;
                    TPZIntPoints *sideIntRule = gel->CreateSideIntegrationRule(iSide, pOrderIntRule);
                    const int npts = sideIntRule->NPoints();

                    TPZTransform<> elTransform(gel->SideToSideTransform(iSide, gel->NSides() - 1));
                    TPZGeoElSide gelSide(gel, iSide);
                    const auto sideDim = gelSide.Dimension();

                    TPZGeoElSide neighGelSide = gelSide.Neighbour();
                    //the following vector will be the edge tg vector if 2D, the normal vector if 3D
                    TPZManVector<REAL,3> vec(3,0);
                    switch(sideDim){
                        case 1:{
                            TPZManVector<int, 2> edgeNodes(2, 0);
                            for (auto i = 0; i < 2; i++) edgeNodes[i] = gel->SideNodeIndex(iSide, i);
                            const REAL sign = edgeNodes[0] < edgeNodes[1] ? 1 : -1;
                            REAL edgeLength = 0;
                            TPZVec<REAL> p0(3), p1(3);
                            gmesh->NodeVec()[edgeNodes[0]].GetCoordinates(p0);
                            gmesh->NodeVec()[edgeNodes[1]].GetCoordinates(p1);
                            edgeLength =
                                    std::sqrt((p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
                                              (p0[2] - p1[2]) * (p0[2] - p1[2]));
                            for (auto x = 0; x < 3; x++) vec[x] = sign * (p1[x] - p0[x]) / edgeLength;
                        }
                            break;
                        case 2:{
                            TPZManVector<REAL,3> xCenter(2,0);
                            switch(gelSide.NSideNodes()){
                                case 3:
                                    pztopology::TPZTriangle::CenterPoint(pztopology::TPZTriangle::NSides-1,xCenter);
                                    break;
                                case 4:
                                    pztopology::TPZQuadrilateral::CenterPoint(pztopology::TPZQuadrilateral::NSides-1,xCenter);
                                    break;
                                default:
                                    DebugStop();
                            }
                            gelSide.Normal(xCenter,vec);
                        }
                            break;
                        default:
                            DebugStop();
                    }
                    bool firstNeighbour{true};


                    //gather all the sides contained in the closure of iSide
                    TPZStack<int> smallSides;
                    gel->LowerDimensionSides(iSide,smallSides);
                    smallSides.Push(iSide);//include the side itself

                    while(neighGelSide != gelSide) {
                        const auto neighCel = dynamic_cast<TPZInterpolatedElement *> (neighGelSide.Element()->Reference());
                        if (!neighCel) {
                            neighGelSide = neighGelSide.Neighbour();
                            continue;
                        }
                        const auto neighGel = neighCel->Reference();
                        const int neighNNodes = MElementType_NNodes(neighGel->Type());
                        const auto neighSide = neighGelSide.Side();
                        TPZTransform<> neighTransform(neighCel->Reference()->SideToSideTransform(neighSide,
                                                                                                 neighCel->Reference()->NSides() -1));
                        TPZTransform<> localTransf(sideDim);
                        gelSide.SideTransform3(neighGelSide,localTransf);
                        TPZMaterialData neighData;
                        neighCel->InitMaterialData(neighData);

                        const auto neighDim = neighGelSide.Element()->Dimension();
                        TPZManVector <REAL,3> pts(sideDim,0),ptsN(sideDim,0), ptEl(dim,0),ptNeigh(neighDim,0);
                        REAL w;
                        TPZFNMatrix<60,REAL> elShape,neighShape;

                        for(auto subSide : smallSides){
                            TPZGeoElSide gelSubSide(gel, subSide);
                            if(gel->SideDimension(subSide) < 1) continue;
                            const int subConnect = subSide - elNNodes;
                            const int neighSubSide = [&](){
                                TPZGeoElSide neighGelSubSide = gelSubSide.Neighbour();
                                while(neighGelSubSide.Element() != neighGelSide.Element()) {
                                    neighGelSubSide = neighGelSubSide.Neighbour();
                                    if(neighGelSubSide.Element() == gel){
                                        DebugStop();
                                    }
                                }
                                return neighGelSubSide.Side();
                            }();

                            const int firstElShape = [&](){
                                int firstElShapeTemp = 0;
                                for(auto jCon = 0; jCon < subConnect; jCon++){
                                    firstElShapeTemp += cel->NConnectShapeF(jCon,cel->EffectiveSideOrder(jCon+elNNodes));
                                }
                                return firstElShapeTemp;
                            }();
                            const int nShapes = cel->NConnectShapeF(subConnect,pOrder);

                            {
                                const auto elConIndex = cel->SideConnectIndex(cel->NSideConnects(subSide) - 1, subSide);
                                const auto neighConIndex = neighCel->SideConnectIndex(neighCel->NSideConnects(neighSubSide) - 1, neighSubSide);
                                const bool check = elConIndex == neighConIndex;
                                BOOST_CHECK_MESSAGE(check, "\n" + testName + " failed" +
                                                           "\ntopology: " + MMeshType_Name(type) + "\n"
                                );
                            }
                            const int firstNeighShape = [&](){
                                int firstNeighShapeTemp = 0;
                                const int neighCon = neighSubSide - neighGelSide.Element()->NNodes();
                                for(auto jCon = 0; jCon < neighCon; jCon++){
                                    firstNeighShapeTemp += neighCel->NConnectShapeF(jCon,neighCel->EffectiveSideOrder(jCon+neighNNodes));
                                }
                                return firstNeighShapeTemp;
                            }();

                            for (auto ipt = 0; ipt < npts; ipt++) {
                                sideIntRule->Point(ipt, pts, w);

                                elTransform.Apply(pts,ptEl);
                                cel->ComputeRequiredData(elData, ptEl);
                                TPZHCurlAuxClass::ComputeShape(elData.fVecShapeIndex, elData.phi,
                                                               elData.fDeformedDirections,elShape);
                                localTransf.Apply(pts,ptsN);
                                neighTransform.Apply(ptsN,ptNeigh);
                                neighCel->ComputeRequiredData(neighData, ptNeigh);
                                TPZHCurlAuxClass::ComputeShape(neighData.fVecShapeIndex, neighData.phi,
                                                               neighData.fDeformedDirections,neighShape);

                                TPZManVector<REAL,3> elShapeFunc(3,0), neighShapeFunc(3,0);
                                bool anyWrongCheck = false;
                                for(auto iShape = 0; iShape < nShapes; iShape ++){
                                    const int elPhiIndex = firstElShape+iShape;
                                    const int neighPhiIndex = firstNeighShape+iShape;
                                    const int elH1phiIndex = elData.fVecShapeIndex[elPhiIndex].second;
                                    const int neighH1phiIndex = neighData.fVecShapeIndex[neighPhiIndex].second;
                                    const bool checkPhis = std::abs(elData.phi(elH1phiIndex,0) - neighData.phi(neighH1phiIndex,0)) < tol;

                                    BOOST_CHECK_MESSAGE(checkPhis,"\n"+testName+" failed: phis are different!"+
                                                                  "\ntopology: "+MMeshType_Name(type)+"\n"+
                                                                  "side: "+std::to_string(iSide)+"\n"+
                                                                  "p order: "+std::to_string(pOrder)+"\n"+
                                                                  "elem phi index: "+std::to_string(elPhiIndex)+"\n"+
                                                                  "neig phi index: "+std::to_string(neighPhiIndex)+"\n"
                                    );
                                    anyWrongCheck = !checkPhis || anyWrongCheck;
                                    if(anyWrongCheck) {
                                        break;
                                    }
                                    for (auto x = 0; x < 3; x++) elShapeFunc[x] = elShape(elPhiIndex,x);
                                    for (auto x = 0; x < 3; x++) neighShapeFunc[x] = neighShape(neighPhiIndex,x);
                                    REAL diffTrace{0};
                                    TPZManVector<REAL,3> elTrace,neighTrace;
                                    const bool checkTraces = CheckTracesFunc(diffTrace,elShapeFunc,neighShapeFunc,vec,sideDim,elTrace,neighTrace);

                                    {
                                        std::ostringstream traceMsg;
                                        if(!checkTraces){
                                            traceMsg <<"el func: ";
                                            for (auto x = 0; x < 3; x++) traceMsg << elShapeFunc[x]<<"\t";
                                            traceMsg<<"\nneigh func: ";
                                            for (auto x = 0; x < 3; x++) traceMsg << neighShapeFunc[x]<<"\t";
                                            traceMsg <<"\n";
                                            traceMsg <<"el trace: ";
                                            for (auto x = 0; x < elTrace.size(); x++) traceMsg << elTrace[x]<<"\t";
                                            traceMsg<<"\nneigh trace: ";
                                            for (auto x = 0; x < neighTrace.size(); x++) traceMsg << neighTrace[x]<<"\t";
                                            traceMsg <<"\n";
                                        }
                                        BOOST_CHECK_MESSAGE(checkTraces,"\n"+testName+" failed"+
                                                                        "\ntopology: "+MMeshType_Name(type)+"\n"+
                                                                        "side: "+std::to_string(iSide)+"\n"+
                                                                        "p order: "+std::to_string(pOrder)+"\n"+
                                                                        "diff traces: "+std::to_string(diffTrace)+"\n"+
                                                                        traceMsg.str()
                                        );
                                    }
                                    anyWrongCheck = !checkTraces || anyWrongCheck;
                                    if(anyWrongCheck) {
                                        break;
                                    }
                                }

                                //now checking sideshapefunction. this needs to be checked only once
                                if(firstNeighbour){
                                    TPZFNMatrix<30,REAL> sideShapeFuncs,sideShapeCurl;
                                    cel->SideShapeFunction(gelSide.Side(),pts,sideShapeFuncs,sideShapeCurl);

                                    const int firstSideShape = [&](){
                                        int firstSideShapeTemp = 0;
                                        TPZStack<int> subSubSides;
                                        gel->LowerDimensionSides(gelSide.Side(),subSubSides);
                                        for(auto subsubSide : subSubSides){
                                            if(gel->SideDimension(subsubSide) < 1) continue;
                                            if(subsubSide == subSide) break;
                                            firstSideShapeTemp += cel->NConnectShapeF(subsubSide - elNNodes,pOrder);
                                        }
                                        return firstSideShapeTemp;
                                    }();
                                    TPZManVector<REAL,3> sideShapeFunc(3,0);
                                    for(auto iShape = 0; iShape < nShapes; iShape ++){
                                        for (auto x = 0; x < 3; x++) elShapeFunc[x] = elShape(firstElShape + iShape,x);
                                        for (auto x = 0; x < 3; x++) sideShapeFunc[x] = sideShapeFuncs(firstSideShape + iShape,x);
                                        REAL diffTrace{0};
                                        TPZManVector<REAL,3> elTrace, sideTrace;
                                        const bool checkTraces = CheckTracesFunc(diffTrace,elShapeFunc,sideShapeFunc,vec,sideDim,elTrace,sideTrace);
                                        {
                                            std::ostringstream traceMsg;
                                            if(!checkTraces){
                                                traceMsg <<"el func: ";
                                                for (auto x = 0; x < 3; x++) traceMsg << elShapeFunc[x]<<"\t";
                                                traceMsg<<"\nside func: ";
                                                for (auto x = 0; x < 3; x++) traceMsg << sideShapeFunc[x]<<"\t";
                                                traceMsg <<"\n";
                                                traceMsg <<"el trace: ";
                                                for (auto x = 0; x < elTrace.size(); x++) traceMsg << elTrace[x]<<"\t";
                                                traceMsg<<"\nside trace: ";
                                                for (auto x = 0; x < elTrace.size(); x++) traceMsg << sideTrace[x]<<"\t";
                                                traceMsg <<"\n";
                                            }
                                            BOOST_CHECK_MESSAGE(checkTraces,"\n"+testName+" failed"+
                                                                            "\ntopology: "+MMeshType_Name(type)+"\n"+
                                                                            "side: "+std::to_string(iSide)+"\n"+
                                                                            "p order: "+std::to_string(pOrder)+"\n"+
                                                                            "diff traces: "+std::to_string(diffTrace)+"\n"+
                                                                            traceMsg.str()
                                            );
                                        }
                                        anyWrongCheck = !checkTraces || anyWrongCheck;
                                    }
                                    if(anyWrongCheck){
                                        break;
                                    }
                                }
                            }
                        }
                        firstNeighbour = false;
                        neighGelSide = neighGelSide.Neighbour();
                    }
                }
            }
            /***********************************************************************************************************
             *              the following lines might be useful for analysing the basis functions
            ***********************************************************************************************************/
//            TPZAnalysis an(cmesh,false);
//            const int postProcessResolution = 3;
//            const std::string executionInfo = [&](){
//                std::string name("");
//                name.append(MMeshType_Name(type));
//                name.append(std::to_string(pOrder));
//                return name;
//            }();
//
//            const std::string plotfile = "solution"+executionInfo+".vtk";//where to print the vtk files
//            TPZStack<std::string> scalnames, vecnames;
//            vecnames.Push("E");//print the state variable
//            auto sol = an.Solution();
//            sol.Zero();
//            for(int i = 0; i < sol.Rows(); i++){
//                sol(i - 1 < 0 ? 0 : i - 1 , 0) = 0;
//                sol(i,0) = 1;
//                an.LoadSolution(sol);
//                an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
//                an.SetStep(i);
//                an.PostProcess(postProcessResolution);
//            }
            delete cmesh;
            delete gmesh;
        }//hcurltest::TestFunctionTracesUniformMesh


        template <class TGEOM>
        void PrintShapeFunctions(const int pOrder){
            constexpr auto nNodes{TGEOM::NCornerNodes};
            constexpr auto dim{TGEOM::Dimension};
            TPZGeoMesh *gmesh = new TPZGeoMesh();
            gmesh->SetDimension(dim);
            //Auxiliary vector to store coordinates:
            TPZManVector<REAL,3> coords(3, 0.);
            for(int iNode = 0 ; iNode < nNodes; iNode++){
                TGEOM::ParametricDomainNodeCoord(iNode,coords);
                auto newindex = gmesh->NodeVec().AllocateNewElement();
                gmesh->NodeVec()[newindex].Initialize(coords, *gmesh);
            }
            TPZManVector<int64_t,TGEOM::NCornerNodes> nodesIdVec(TGEOM::NCornerNodes,-1);
            for(int i = 0 ; i < TGEOM::NCornerNodes; i++) nodesIdVec[i] = i;
            int64_t index{-1};
            TPZGeoEl *gel = gmesh->CreateGeoElement(TGEOM::Type(), nodesIdVec, 1, index, 0);
            gmesh->BuildConnectivity();
            TPZCompMesh *cmesh = new TPZCompMesh (gmesh);

            cmesh->SetDefaultOrder(pOrder);
            cmesh->SetDimModel(dim);

            auto mat = new TPZMatHCurlProjection(dim,1,1);
            cmesh->InsertMaterialObject(mat);
            cmesh->SetAllCreateFunctionsHCurl();
            cmesh->AutoBuild();
            cmesh->CleanUpUnconnectedNodes();

            TPZAnalysis an(cmesh,false);
            const int postProcessResolution = 3;
            const std::string executionInfo = [&](){
                std::string name("");
                name.append(MElementType_Name(TGEOM::Type()));
                name.append(std::to_string(pOrder));
                return name;
            }();

            const std::string plotfile = "shapeFuncs"+executionInfo+".vtk";//where to print the vtk files
            TPZStack<std::string> scalnames, vecnames;
            vecnames.Push("E");//print the state variable
            auto sol = an.Solution();
            sol.Zero();
            for(int i = 0; i < sol.Rows(); i++){
                sol(i - 1 < 0 ? 0 : i - 1 , 0) = 0;
                sol(i,0) = 1;
                an.LoadSolution(sol);
                an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
                an.SetStep(i);
                an.PostProcess(postProcessResolution);
            }
            delete cmesh;
            delete gmesh;
        }


        namespace auxiliaryfuncs{
            void ComputeDirections (TPZGeoEl *gel, TPZFMatrix<REAL> &deformedDirections,
                    TPZFMatrix<REAL> &masterDirections, const TPZVec<int> &transformationIds) {
                int dim = gel->Dimension();
                TPZFMatrix<REAL> gradX(dim, dim, 0);
                for (auto x = 0; x < dim; x++) gradX(x, x) = 1;
                const int nVec = dim * gel->NSides();
                masterDirections.Resize(dim, nVec);
                deformedDirections.Resize(3, nVec);
                switch (gel->Type()) {
                    case EOned:
                        pztopology::TPZLine::ComputeHCurlDirections(gradX,
                                                                        masterDirections,transformationIds);//these are the vectors on the master element
                        break;
                    case ETriangle:
                        pztopology::TPZTriangle::ComputeHCurlDirections(gradX,
                                                                        masterDirections,transformationIds);//these are the vectors on the master element
                        break;
                    case EQuadrilateral:
                        pztopology::TPZQuadrilateral::ComputeHCurlDirections(gradX,
                                                                             masterDirections,transformationIds);//these are the vectors on the master element
                        break;
                    case ETetraedro:
                        pztopology::TPZTetrahedron::ComputeHCurlDirections(gradX,
                                                                           masterDirections,transformationIds);//these are the vectors on the master element
                        break;
                    case EPiramide:
                        DebugStop();//pztopology::TPZPyramid::ComputeHCurlDirections(gradX, masterDirections,transformationIds);//these are the vectors on the master element
                        break;
                    case EPrisma:
                        pztopology::TPZPrism::ComputeHCurlDirections(gradX, masterDirections,transformationIds);//these are the vectors on the master element
                        break;
                    case ECube:
                        pztopology::TPZCube::ComputeHCurlDirections(gradX,
                                                                    masterDirections,transformationIds);//these are the vectors on the master element
                        break;
                    default:
                        DebugStop();
                }
                TPZVec<REAL> xiCenter(dim, 0);
                TPZFMatrix<REAL> jacobian(dim, dim, 0), axes(dim, 3), jacInv(dim, dim, 0);
                REAL detJac = 0;
                gel->Jacobian(xiCenter, jacobian, axes, detJac, jacInv);

#ifdef NOISY_HCURL
                std::cout << std::endl;
                jacobian.Print("Original Jacobian:");
                axes.Print("Original Axes:");
#endif
                for (auto iVec = 0; iVec < nVec; iVec++) {

                    TPZManVector<REAL, 3> tempDirection(dim, 0);
                    for (auto i = 0; i < dim; i++) {
                        //covariant piola transform: J^{-T}
                        tempDirection[i] = 0;
                        for (auto j = 0; j < dim; j++) tempDirection[i] += jacInv(j, i) * masterDirections(j, iVec);
                    }
                    for (auto i = 0; i < 3; i++) {
                        deformedDirections(i, iVec) = 0;
                        for (auto j = 0; j < dim; j++) deformedDirections(i, iVec) += axes(j, i) * tempDirection[j];
                    }
                }
            }//ComputeDirections

            void VectorProduct(const TPZVec<REAL> &v1, const TPZVec<REAL> &v2,TPZVec<REAL> &result){
                REAL x1=v1[0], y1=v1[1],z1=v1[2];
                REAL x2=v2[0], y2=v2[1],z2=v2[2];
                result.Resize(v1.NElements());
                result[0]=y1*z2-z1*y2;
                result[1]=z1*x2-x1*z2;
                result[2]=x1*y2-y1*x2;
            };//VectorProduct

            TPZGeoMesh *CreateGeoMesh2D(int nelx, int nely, MMeshType meshType, TPZVec<int> &matIds) {
                //Creating geometric mesh, nodes and elements.
                //Including nodes and elements in the mesh object:
                TPZGeoMesh *gmesh = new TPZGeoMesh();
                constexpr int dim{2};
                gmesh->SetDimension(dim);

                //Auxiliary vector to store coordinates:
                TPZVec<REAL> coord1(3, 0.);
                TPZVec<REAL> coord2(3, 0.);
                coord1[0] = 0;coord1[1] = 0;coord1[2] = 0;
                coord2[0] = 1;coord2[1] = 1;coord2[2] = 0;

                TPZManVector<int> nelem(2, 1);
                nelem[0] = nelx;
                nelem[1] = nely;

                TPZGenGrid2D gengrid(nelem, coord1, coord2);
                gengrid.SetElementType(meshType);
                constexpr int matIdDomain = 1, matIdBoundary = 2;
                gengrid.Read(gmesh, matIdDomain);
                gengrid.SetBC(gmesh, 4, matIdBoundary);
                gengrid.SetBC(gmesh, 5, matIdBoundary);
                gengrid.SetBC(gmesh, 6, matIdBoundary);
                gengrid.SetBC(gmesh, 7, matIdBoundary);

                gmesh->BuildConnectivity();

                {
                    TPZCheckGeom check(gmesh);
                    check.CheckUniqueId();
                }
                matIds.Resize(2);
                matIds[0] = matIdDomain;
                matIds[1] = matIdBoundary;
                //Printing geometric mesh:

                //ofstream bf("before.vtk");
                //TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
                return gmesh;
            }//CreateGeoMesh2D

            TPZGeoMesh *CreateGeoMesh3D(int nelx, int nely, int nelz, MMeshType meshType, TPZVec<int> &matIds) {
                //Creating geometric mesh, nodes and elements.
                //Including nodes and elements in the mesh object:
                //create boundary elements
                TPZManVector<REAL,3>minX(3,0);
                TPZManVector<REAL,3>maxX(3,1);
                TPZManVector<int,3> nelDiv(3,-1);
                nelDiv[0] = nelx;
                nelDiv[1] = nely;
                nelDiv[2] = nelz;
                constexpr int matIdDomain{1};
                constexpr int matIdBoundary{2};

                TPZGenGrid3D genGrid3D(minX,maxX,nelDiv,meshType);
                genGrid3D.BuildVolumetricElements(matIdDomain);
                TPZGeoMesh *gmesh = genGrid3D.BuildBoundaryElements(matIdBoundary,matIdBoundary,matIdBoundary,matIdBoundary,matIdBoundary,matIdBoundary);
                gmesh->BuildConnectivity();

                {
                    TPZCheckGeom check(gmesh);
                    check.CheckUniqueId();
                    check.PerformCheck();
                }
                matIds.Resize(2);
                matIds[0] = matIdDomain;
                matIds[1] = matIdBoundary;
                //Printing geometric mesh:

                //ofstream bf("before.vtk");
                //TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
                return gmesh;
            }//CreateGeoMesh3D
        }//namespace
    }//namespace


BOOST_AUTO_TEST_SUITE_END()

#ifdef NOISY_HCURL
#undef NOISY_HCURL
#endif
#ifdef NOISYVTK_HCURL
#undef NOISYVTK_HCURL
#endif

#endif

