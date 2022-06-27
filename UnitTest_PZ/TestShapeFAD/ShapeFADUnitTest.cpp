/**
 * @file BlendUnitTest.cpp
 * @brief Define a Unit Test using Boost for blend elements over all kind of geometries available in the library
 *
 */
#include <iostream>

#include "TPZRefPatternDataBase.h"


#include "TPZVTKGeoMesh.h"


#include "TPZGenSpecialGrid.h"
#include "TPZShapeH1.h"
#include "TPZShapeHDiv.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "TPZMaterialDataT.h"
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.testgeom");
#endif

#include "fad.h"


#include<catch2/catch.hpp>

#define SHAPEFAD_VERBOSE //outputs x and grad comparisons
#define SHAPEFAD_OUTPUT_TXT//prints all elements in .txt format
#define SHAPEFAD_OUTPUT_VTK//prints all elements in .vtk format

namespace shapetest{
const int pOrder = 3;
const REAL tol = 1e-8;
void TestMesh2D(TPZGeoMesh *gmesh, int nDiv);
void TestMesh3D(TPZGeoMesh *gmesh, int nDiv);
//check the fobrenius norm of the difference between two given matrices
bool CheckMatrices(const TPZFMatrix<REAL> &gradx1, std::string name1,
                   const TPZFMatrix<REAL> &gradx2, std::string name2, const REAL &tol);

}

TEST_CASE("shapefad_tests","[shape_tests]") {
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    {
        const int nDiv = 4;
        auto gmesh = TPZGenSpecialGrid::CreateGeoMesh2D_Circle(0);
        shapetest::TestMesh2D(gmesh,nDiv);
        delete gmesh;
    }
    
    gRefDBase.InitializeUniformRefPattern(ETetraedro);
    gRefDBase.InitializeUniformRefPattern(EPiramide);
    gRefDBase.InitializeUniformRefPattern(EPrisma);
    gRefDBase.InitializeUniformRefPattern(ECube);
    {
        const int nDiv = 2;
        auto gmesh = TPZGenSpecialGrid::CreateGeoMesh3D_DividedSphere(0);
        shapetest::TestMesh3D(gmesh, nDiv);
        delete gmesh;
    }
}


namespace shapetest{



void TestMesh2D(TPZGeoMesh *gmesh, int nDiv)
{
    if(gmesh->Dimension() != 2) DebugStop();
    {
        TPZManVector<REAL,3> xiReal;
        TPZManVector<Fad<REAL>,3> xiFad;
        REAL weight = -1;//useless
        const int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *geo = gmesh->ElementVec()[iel];
            if (geo && geo->Dimension()==2 && !geo->HasSubElement()) {
                uint64_t errors = 0;
                int nnodes = geo->NNodes();
                TPZManVector<int64_t,4> nodeids(nnodes);
                for (int in = 0; in<nnodes; in++) {
                    nodeids[in]= geo->NodePtr(in)->Id();
                }
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
                    TPZFNMatrix<9,REAL> gradXfad(gradXreal.Rows(),gradXreal.Cols(),0.0);
                    for(int i = 0; i < gradXreal.Rows(); i++){
                        for(int j = 0; j < gradXreal.Cols(); j++){
                            gradXfad(i,j) = xFad[i].dx(j);
                        }
                    }
                    hasAnErrorOccurred = shapetest::CheckMatrices(gradXreal,"gradXreal",gradXfad,"gradXfad",shapetest::tol);
                    REQUIRE(!hasAnErrorOccurred);
                    if(hasAnErrorOccurred){
                        errors++;
                    }
                    // initialize the value of xi with derivatives in xy
                    TPZManVector<Fad<REAL>,2> qsifad(2);
                    TPZFNMatrix<9,Fad<REAL>> gradXFAD(2,2);
                    geo->ComputeQsiGrad(xiReal, qsifad);
                    geo->GradX(qsifad, gradXFAD);
                    Fad<REAL> detjac = gradXFAD(0,0)*gradXFAD(1,1)-gradXFAD(0,1)*gradXFAD(1,0);
                    TPZFNMatrix<9,REAL> gradx(2,2),gradinv(2,2),gradxLU(2,2);
                    for (int i=0; i<2; i++) {
                        for (int j=0; j<2; j++) {
                            gradx(i,j) = gradXreal(i,j);
                        }
                    }
//                    gradx.Print("gradx = ",std::cout,EMathematicaInput);
                    gradxLU = gradx;
                    gradxLU.Inverse(gradinv, ELU);
//                    gradinv.Print("gradinv = ",std::cout,EMathematicaInput);
//                    std::cout << "qsifad" << qsifad << std::endl;
                    TPZMaterialDataT<REAL> data;
                    TPZManVector<int, 4> orders(geo->NSides()-geo->NCornerNodes(),shapetest::pOrder);
                    TPZFNMatrix<25,REAL> dphix,dphi;
                    TPZFNMatrix<25,Fad<REAL>> phiFAD;
                    TPZManVector<int,4> sideorient(geo->NSides(1),1);
                    if(geo->Type() == ETriangle)
                    {
                        TPZShapeH1<pzshape::TPZShapeTriang> shapeh1;
                        shapeh1.Initialize(nodeids, orders, data);
                        phiFAD.Resize(data.fPhi.Rows(),1);
                        TPZFMatrix<REAL> phi;
                        shapeh1.Shape(xiReal, data, phi, dphi);
                        TPZFMatrix<Fad<REAL>> dphiFAD(2,data.fDPhi.Cols());
                        shapeh1.Shape(qsifad, data, phiFAD, dphiFAD);
                    }
                    if(geo->Type() == EQuadrilateral)
                    {
                        TPZShapeH1<pzshape::TPZShapeQuad> shapeh1;
                        shapeh1.Initialize(nodeids, orders, data);
                        phiFAD.Resize(data.fPhi.Rows(),1);
                        TPZFMatrix<REAL> phi;
                        shapeh1.Shape(xiReal, data, phi, dphi);
//                        dphi.Print("dphi = ",std::cout, EMathematicaInput);
                        TPZFMatrix<Fad<REAL>> dphiFAD(2,data.fDPhi.Cols());
                        shapeh1.Shape(qsifad, data, phiFAD, dphiFAD);
//                        phiFAD.Print("phiFAD",std::cout);
                    }
                    gradinv.MultAdd(dphi, dphi, dphix, 1., 0., 1);
                    TPZFNMatrix<50,REAL> dphixFAD(dphix);
                    for (int i=0; i<dphixFAD.Cols(); i++) {
                        dphixFAD(0,i) = phiFAD(i,0).dx(0);
                        dphixFAD(1,i) = phiFAD(i,0).dx(1);
                    }
                    hasAnErrorOccurred = shapetest::CheckMatrices(dphix,"dphiXreal",dphixFAD,"dphiXfad",shapetest::tol);
                    REQUIRE(!hasAnErrorOccurred);
                    if(hasAnErrorOccurred){
                        errors++;
                    }
                    if(geo->Type() == ETriangle)
                    {
                        TPZShapeHDiv<pzshape::TPZShapeTriang> shapehdiv;
                        shapehdiv.Initialize(nodeids, orders, sideorient, data);
                        phiFAD.Resize(data.fPhi.Rows(),1);
                        TPZFMatrix<REAL> phi, divphi;
                        shapehdiv.Shape(xiReal, data, phi, divphi);
                        TPZFMatrix<Fad<REAL>> dphiFAD(2,data.fDPhi.Cols());
                        shapehdiv.Shape(qsifad, data, phiFAD, dphiFAD);
                        TPZFMatrix<REAL> phiFADREAL(phi);
                        for (int i=0; i<phi.Rows(); i++) {
                            for (int j=0; j<phi.Cols(); j++) {
                                phiFADREAL(i,j) = phiFAD(i,j).val();
                            }
                        }
                        hasAnErrorOccurred = shapetest::CheckMatrices(phi,"hdivphi",phiFADREAL,"hdivphiFAD",shapetest::tol);
                        REQUIRE(!hasAnErrorOccurred);
                        if(hasAnErrorOccurred){
                            errors++;
                        }
                        TPZFMatrix<Fad<REAL>> dirHDiv;
                        gradXFAD.MultAdd(phiFAD, phiFAD, dirHDiv, 1./detjac, 0.);
                        TPZFMatrix<REAL> divmaster(dirHDiv.Cols(),1);
                        for (int i = 0; i<divmaster.Rows(); i++) {
                            divmaster(i,0) = (dirHDiv(0,i).dx(0)+dirHDiv(1,i).dx(1))*detjac.val();
                        }
//                        std::cout << "div phi master " << divphi(0,0) << std::endl;
//                        std::cout << "div computed fad " << divmaster(0,0) << std::endl;
//                        std::cout << "ratio ";
//                        for(int i=0; i<10; i++) std::cout << divphi(i,0)/divmaster(i,0) << " ";
//                        std::cout  << std::endl;
                        hasAnErrorOccurred = shapetest::CheckMatrices(divmaster,"div master",divphi,"div computed FAD",shapetest::tol);
                        REQUIRE(!hasAnErrorOccurred);
                        if(hasAnErrorOccurred){
                            errors++;
                        }
                    }
                    if(geo->Type() == EQuadrilateral)
                    {
                        TPZShapeHDiv<pzshape::TPZShapeQuad> shapehdiv;
                        shapehdiv.Initialize(nodeids, orders, sideorient, data);
                        phiFAD.Resize(data.fPhi.Rows(),1);
                        TPZFMatrix<REAL> phi,divphi;
                        shapehdiv.Shape(xiReal, data, phi, divphi);
                        TPZFMatrix<Fad<REAL>> dphiFAD(2,data.fDPhi.Cols());
                        shapehdiv.Shape(qsifad, data, phiFAD, dphiFAD);
//                        phi.Print("phihdiv ",std::cout);
//                        phiFAD.Print("phiFAD hdiv ",std::cout);
                        TPZFMatrix<REAL> phiFADREAL(phi);
                        for (int i=0; i<phi.Rows(); i++) {
                            for (int j=0; j<phi.Cols(); j++) {
                                phiFADREAL(i,j) = phiFAD(i,j).val();
                            }
                        }
                        hasAnErrorOccurred = shapetest::CheckMatrices(phi,"hdivphi",phiFADREAL,"hdivphiFAD",shapetest::tol);
                        REQUIRE(!hasAnErrorOccurred);
                        if(hasAnErrorOccurred){
                            errors++;
                        }
                        TPZFMatrix<Fad<REAL>> dirHDiv;
                        gradXFAD.MultAdd(phiFAD, phiFAD, dirHDiv, 1./detjac, 0.);
                        TPZFMatrix<REAL> divmaster(dirHDiv.Cols(),1);
                        for (int i = 0; i<divmaster.Rows(); i++) {
                            divmaster(i,0) = (dirHDiv(0,i).dx(0)+dirHDiv(1,i).dx(1))*detjac.val();
                        }
//                        std::cout << "div phi master " << divphi(0,0) << std::endl;
//                        std::cout << "div computed fad " << divmaster(0,0) << std::endl;
//                        std::cout << "ratio ";
//                        for(int i=0; i<10; i++) std::cout << divphi(i,0)/divmaster(i,0) << " ";
//                        std::cout  << std::endl;
                        hasAnErrorOccurred = shapetest::CheckMatrices(divmaster,"div master",divphi,"div computed FAD",shapetest::tol);
                        REQUIRE(!hasAnErrorOccurred);
                        if(hasAnErrorOccurred){
                            errors++;
                        }
                    }

                }
#ifdef SHAPEFAD_VERBOSE
                if(errors > 0 || geo->IsGeoBlendEl()){
                    std::cout << "============================"
                    << std::endl;
                    std::cout
                    << "Element: " << geo->Id()
                    << "\tType: " << MElementType_Name(geo->Type());
                    std::cout << "\tIs blend? : " << geo->IsGeoBlendEl()
                    << std::endl;
                    std::cout
                    << "\tNumber of points: " << intRule->NPoints()
                    << "\tErrors: " << errors << std::endl;
                }
#endif
                delete intRule;
            }
        }
    }
    
    
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
            std::cout<<"\b";
        }
    }
#ifdef SHAPEFAD_OUTPUT_TXT
    {
        const std::string meshFileName =
        "blendmesh2D.txt";
        std::ofstream outTXT(meshFileName.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }
#endif
#ifdef SHAPEFAD_OUTPUT_VTK
    {
        const std::string meshFileName =
        "blendmesh2D.vtk";
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
    }
#endif
}


void TestMesh3D(TPZGeoMesh *gmesh, int nDiv)
{
    {
        TPZManVector<REAL,3> xiReal;
        TPZManVector<Fad<REAL>,3> xiFad;
        REAL weight = -1;//useless
        const int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *geo = gmesh->ElementVec()[iel];
            if (geo && !geo->HasSubElement()) {
                uint64_t errors = 0;
                auto intRule = geo->CreateSideIntegrationRule(geo->NSides()-1, pOrder);
                xiReal.Resize(geo->Dimension(),0);
                xiFad.Resize(geo->Dimension(),0);
                for(int iPt = 0; iPt < intRule->NPoints(); iPt++) {
                    bool hasAnErrorOccurred = false;
                    intRule->Point(iPt, xiReal, weight);
                    for (int x = 0; x < geo->Dimension(); x++) {
                        xiFad[x] = Fad<REAL>(geo->Dimension(), x, xiReal[x]);
                    }
                    
                    TPZManVector<Fad<REAL>, 3> xFad(3);
                    geo->X(xiFad, xFad);
                    TPZFNMatrix<9, REAL> gradXreal(3, 3, 0);
                    geo->GradX(xiReal, gradXreal);
                    TPZFNMatrix<9,REAL> gradXfad(gradXreal.Rows(),gradXreal.Cols(),0.0);
                    for (int i = 0; i < gradXreal.Rows(); i++) {
                        for (int j = 0; j < gradXreal.Cols(); j++) {
                            gradXfad(i, j) = xFad[i].dx(j);
                        }
                    }
                    hasAnErrorOccurred = shapetest::CheckMatrices(gradXreal, "gradXreal", gradXfad, "gradXfad",
                                                                  shapetest::tol);
                    REQUIRE(!hasAnErrorOccurred);
                    if (hasAnErrorOccurred) {
                        errors++;
                    }
                }
#ifdef SHAPEFAD_VERBOSE
                if(errors > 0 || geo->IsGeoBlendEl()){
                    std::cout << "============================"
                    << std::endl;
                    std::cout
                    << "Element: " << geo->Id()
                    << "\tType: " << MElementType_Name(geo->Type());
                    std::cout << "\tIs blend? : " << geo->IsGeoBlendEl()
                    << std::endl;
                    std::cout
                    << "\tNumber of points: " << intRule->NPoints()
                    << "\tErrors: " << errors << std::endl;
                }
#endif
            }
        }
    }
    {
        TPZVec<TPZGeoEl *> sons;
        std::vector<std::string> loading = {"-","/","|","\\"};
        for (int iDiv = 0; iDiv < nDiv; iDiv++) {
#ifdef SHAPEFAD_OUTPUT_TXT
            {
                const std::string meshFileName =
                "blendmesh3D_ref" + std::to_string(iDiv) + ".txt";
                std::ofstream outTXT(meshFileName.c_str());
                gmesh->Print(outTXT);
                outTXT.close();
            }
#endif
            std::cout<<"Performing "<<iDiv+1<<" ref step out of " << nDiv<<std::endl;
            const int nel = gmesh->NElements();
            for (int iel = 0; iel < nel; iel++) {
                std::cout<<"\b"<<loading[iel%4]<<std::flush;
                TPZGeoEl *geo = gmesh->ElementVec()[iel];
                if (geo && !geo->HasSubElement()) {
                    geo->Divide(sons);
                }
            }
            std::cout<<"\b";
        }
    }
#ifdef SHAPEFAD_OUTPUT_TXT
    {
        const std::string meshFileName =
        "blendmesh3D.txt";
        std::ofstream outTXT(meshFileName.c_str());
        gmesh->Print(outTXT);
        outTXT.close();
    }
#endif
#ifdef SHAPEFAD_OUTPUT_VTK
    {
        const std::string meshFileName =
        "blendmesh3D.vtk";
        std::ofstream outVTK(meshFileName.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        outVTK.close();
    }
#endif
}

bool CheckMatrices(const TPZFMatrix<REAL> &gradx1, std::string name1,
                   const TPZFMatrix<REAL> &gradx2, std::string name2, const REAL &tol){
    const auto VAL_WIDTH = 15;
    REAL diff = 0;
    bool hasAnErrorOccurred = false;
    if((gradx1.Rows() != gradx2.Rows()) || (gradx1.Cols() != gradx2.Cols())){
        hasAnErrorOccurred = true;
        return hasAnErrorOccurred;
    }
    for (int i = 0; i < gradx1.Rows(); i++) {
        for (int j = 0; j < gradx1.Cols(); j++) {
            diff += (gradx1.GetVal(i, j) - gradx2.GetVal(i, j)) * (gradx1.GetVal(i, j) - gradx2.GetVal(i, j));
        }
    }
    if (diff > tol) {
        hasAnErrorOccurred = true;
    }
    if (hasAnErrorOccurred) {
#ifdef BLEND_VERBOSE
        std::ostringstream x1m, x2m;
        x1m << name1 << std::endl;
        x2m << name2 << std::endl;
        for (int i = 0; i < gradx1.Rows(); i++) {
            for (int j = 0; j < gradx1.Cols(); j++) {
                x1m << std::setw(VAL_WIDTH) << std::right
                << gradx1.GetVal(i, j) << "\t";
                x2m << std::setw(VAL_WIDTH) << std::right
                << gradx2.GetVal(i, j) << "\t";
            }
            x1m << '\n';
            x2m << '\n';
        }
        std::cout << x1m.str() << '\n';
        std::cout << x2m.str() << '\n';
        std::cout << "diff :" << diff << std::endl;
#endif
    }
    return hasAnErrorOccurred;
}


}
