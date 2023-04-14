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
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"
#include "pzshapetetra.h"
#include "TPZMaterialDataT.h"
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.testgeom");
#endif

#include "fad.h"


#include <catch2/catch_test_macros.hpp>

//#define SHAPEFAD_VERBOSE //outputs x and grad comparisons
//#define SHAPEFAD_OUTPUT_TXT//prints all elements in .txt format
//#define SHAPEFAD_OUTPUT_VTK//prints all elements in .vtk format

namespace shapetest{
const int pOrder = 3;
const REAL tol = 1e-8;

template<class TSHAPE>
void TestMesh(TPZGeoMesh *gmesh, int nDiv);
void TestMesh3D(TPZGeoMesh *gmesh, int nDiv);
template<class TSHAPE>
void ComputePhiFAD(const TPZVec<Fad<REAL>> &pt, TPZVec<int64_t> &nodeids, TPZFMatrix<REAL> &dphi, TPZFMatrix<Fad<REAL>> &phiFAD);

// compute the HDiv functions with standard point e computed with Fad
template<class TSHAPE>
void ComputeHDivFunctions(const TPZVec<Fad<REAL>> &pt, TPZVec<int64_t> &nodeids, TPZFMatrix<REAL> &phiHDiv, TPZFMatrix<REAL> &phiHDivFAD);

// verify if the divergence of the HDiv functions is conforming the Piola transform
template<class TSHAPE>
void VerifyDivergenceCompatibility(const TPZVec<Fad<REAL>> &pt, const TPZFMatrix<Fad<REAL>> &gradX, TPZVec<int64_t> &nodeids, TPZFMatrix<REAL> &divphimaster, TPZFMatrix<REAL> &divphicomputed);

//check the fobrenius norm of the difference between two given matrices
bool CheckMatrices(const TPZFMatrix<REAL> &gradx1, std::string name1,
                   const TPZFMatrix<REAL> &gradx2, std::string name2, const REAL &tol);

template<class T>
T DetJac(const TPZFMatrix<T> &mat);

}

TEST_CASE("shapefad_tests","[shape_tests]") {
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETetraedro);
    gRefDBase.InitializeUniformRefPattern(EPiramide);
    gRefDBase.InitializeUniformRefPattern(EPrisma);
    gRefDBase.InitializeUniformRefPattern(ECube);
    {
        const int nDiv = 0;
        auto gmesh = TPZGenSpecialGrid::CreateGeoMesh3D_DividedSphere(0);        
        shapetest::TestMesh<pzshape::TPZShapeTetra>(gmesh, nDiv);
        shapetest::TestMesh<pzshape::TPZShapeCube>(gmesh, nDiv);
        shapetest::TestMesh<pzshape::TPZShapePrism>(gmesh, nDiv);
        shapetest::TestMesh<pzshape::TPZShapePiram>(gmesh, nDiv);
        delete gmesh;
    }
    {
        const int nDiv = 0;
        auto gmesh = TPZGenSpecialGrid::CreateGeoMesh2D_Circle(0);
        
        shapetest::TestMesh<pzshape::TPZShapeQuad>(gmesh,nDiv);
        shapetest::TestMesh<pzshape::TPZShapeTriang>(gmesh,nDiv);
        delete gmesh;
    }
    
}


namespace shapetest{


template<class TSHAPE>
void TestMesh(TPZGeoMesh *gmesh, int nDiv)
{
    if(gmesh->Dimension() != TSHAPE::Dimension) DebugStop();
    {
        const int dim = TSHAPE::Dimension;
        TPZManVector<REAL,3> xiReal;
        TPZManVector<Fad<REAL>,3> xiFad;
        REAL weight = -1;//useless
        const int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *geo = gmesh->ElementVec()[iel];
            if(!geo) continue;
            if(geo->HasSubElement()) continue;
            if(geo->Type() != TSHAPE::Type()) continue;
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
                TPZManVector<Fad<REAL>,3> qsifad(dim);
                TPZFNMatrix<9,Fad<REAL>> gradXFAD(3,dim);
                geo->ComputeQsiGrad(xiReal, qsifad);
                geo->GradX(qsifad, gradXFAD);
                TPZFNMatrix<9,REAL> gradx(dim,dim),gradinv(dim,dim),gradxLU(dim,dim);
                for (int i=0; i<dim; i++) {
                    for (int j=0; j<dim; j++) {
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
                ComputePhiFAD<TSHAPE>(qsifad, nodeids, dphi, phiFAD);
                gradinv.MultAdd(dphi, dphi, dphix, 1., 0., 1);
                TPZFNMatrix<50,REAL> dphixFAD(dphix);
                for (int i=0; i<dphixFAD.Cols(); i++) {
                    for(int d=0; d<dim; d++)
                    {
                        dphixFAD(d,i) = phiFAD(i,0).dx(d);
                    }                }
                hasAnErrorOccurred = shapetest::CheckMatrices(dphix,"dphiXreal",dphixFAD,"dphiXfad",shapetest::tol);
                REQUIRE(!hasAnErrorOccurred);
                if(hasAnErrorOccurred){
                    errors++;
                }
                
                if(geo->Type() != EPiramide)
                {
                    TPZFMatrix<REAL> phihdiv,phihdivFad;
                    ComputeHDivFunctions<TSHAPE>(qsifad, nodeids, phihdiv, phihdivFad);
                    hasAnErrorOccurred = shapetest::CheckMatrices(phihdiv,"hdivphi",phihdivFad,"hdivphiFAD",shapetest::tol);
                    REQUIRE(!hasAnErrorOccurred);
                    if(hasAnErrorOccurred){
                        errors++;
                    }
                    TPZFMatrix<REAL> divphi,divphiFad;
                    
                    VerifyDivergenceCompatibility<TSHAPE>(qsifad, gradXFAD, nodeids, divphi, divphiFad);

//                    divphi.Print("divphi master",std::cout);
//
//                    divphiFad.Print("divphi deformed", std::cout);
//                    TPZFMatrix<REAL> diff = divphi-divphiFad;
//                    diff.Print("diff ", std::cout);
                    
                    hasAnErrorOccurred = shapetest::CheckMatrices(divphi,"hdivphi",divphiFad,"hdivphiFAD",shapetest::tol);
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
    
    
    {
        TPZGeoMesh copy(*gmesh);
        TPZVec<TPZGeoEl *> sons;
        std::vector<std::string> loading = {"-","/","|","\\"};
        for (int iDiv = 0; iDiv < nDiv; iDiv++) {
            std::cout<<"Performing "<<iDiv+1<<" ref step out of " << nDiv<<std::endl;
            const int nel = copy.NElements();
            for (int iel = 0; iel < nel; iel++) {
                std::cout<<"\b"<<loading[iel%4]<<std::flush;
                TPZGeoEl *geo = copy.ElementVec()[iel];
                if (geo && !geo->HasSubElement()) {
                    geo->Divide(sons);
                }
            }
            std::cout<<"\b";
        }
#ifdef SHAPEFAD_OUTPUT_TXT
        {
            const std::string meshFileName =
            "blendmesh2D.txt";
            std::ofstream outTXT(meshFileName.c_str());
            copy.Print(outTXT);
            outTXT.close();
        }
#endif
#ifdef SHAPEFAD_OUTPUT_VTK
        {
            const std::string meshFileName =
            "blendmesh2D.vtk";
            std::ofstream outVTK(meshFileName.c_str());
            TPZVTKGeoMesh::PrintGMeshVTK(&copy, outVTK, true);
            outVTK.close();
        }
#endif
    }
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

template<class TSHAPE>
void ComputePhiFAD(const TPZVec<Fad<REAL>> &pt, TPZVec<int64_t> &nodeids, TPZFMatrix<REAL> &dphi, TPZFMatrix<Fad<REAL>> &phiFAD)
{
    TPZManVector<int,27> orders(TSHAPE::NSides-TSHAPE::NCornerNodes,pOrder);
    TPZShapeH1<TSHAPE> shapeh1;
    TPZShapeData data;
    shapeh1.Initialize(nodeids, orders, data);
    
    TPZFMatrix<REAL> phi(data.fPhi);
    dphi.Resize(data.fDPhi.Rows(),data.fDPhi.Cols());
    TPZManVector<REAL> xiReal(pt.size());
    for(int i=0; i<xiReal.size(); i++) xiReal[i] = pt[i].val();
    shapeh1.Shape(xiReal, data, phi, dphi);
    phiFAD.Resize(data.fPhi.Rows(),1);
    TPZFMatrix<Fad<REAL>> dphiFAD(pt.size(),data.fDPhi.Cols());
    shapeh1.Shape(pt, data, phiFAD, dphiFAD);
}

// compute the HDiv functions with standard point e computed with Fad
template<class TSHAPE>
void ComputeHDivFunctions(const TPZVec<Fad<REAL>> &qsifad, TPZVec<int64_t> &nodeids, TPZFMatrix<REAL> &phiHDiv, TPZFMatrix<REAL> &phiHDivFADReal)
{
    const int nfaces = TSHAPE::NFacets;
    TPZManVector<int,27> orders(nfaces+1,pOrder);
    const int dim = TSHAPE::Dimension;
    TPZManVector<int,27> sideorient(nfaces,1);
    TPZShapeHDiv<TSHAPE> shapehdiv;
    TPZShapeData data;
    shapehdiv.Initialize(nodeids, orders, sideorient, data);
    TPZFMatrix<REAL> divphi;
    TPZManVector<REAL,3> xiReal(dim);
    for(int i=0; i<dim; i++) xiReal[i]=qsifad[i].val();
    shapehdiv.Shape(xiReal, data, phiHDiv, divphi);
    TPZFMatrix<Fad<REAL>> phiHDivFAD,divphiFAD;
    shapehdiv.Shape(qsifad, data, phiHDivFAD, divphiFAD);
    phiHDivFADReal.Resize(phiHDivFAD.Rows(),phiHDivFAD.Cols());
    for (int i=0; i<phiHDiv.Rows(); i++) {
        for (int j=0; j<phiHDiv.Cols(); j++) {
            phiHDivFADReal(i,j) = phiHDivFAD(i,j).val();
        }
    }
}


template<class T>
T DetJac(const TPZFMatrix<T> &mat)
{
    if(mat.Cols() == 2) return mat(0,0)*mat(1,1)-mat(0,1)*mat(1,0);
    if(mat.Cols() == 3)
    {
        return mat(0,0)*mat(1,1)*mat(2,2)+mat(0,1)*mat(1,2)*mat(2,0)+mat(1,0)*mat(2,1)*mat(0,2)
        -mat(0,2)*mat(1,1)*mat(2,0)-mat(1,0)*mat(0,1)*mat(2,2)-mat(2,1)*mat(1,2)*mat(0,0);

    }
    DebugStop();
    return T(0);
}

template<class TSHAPE>
void VerifyDivergenceCompatibility(const TPZVec<Fad<REAL>> &qsifad, const TPZFMatrix<Fad<REAL>> &gradX, TPZVec<int64_t> &nodeids, TPZFMatrix<REAL> &divphimaster, TPZFMatrix<REAL> &divphicomputed)
{
    const int nfaces = TSHAPE::NFacets;
    TPZManVector<int,27> orders(nfaces+1,pOrder);
    const int dim = TSHAPE::Dimension;
    TPZManVector<int,27> sideorient(nfaces,1);
    TPZShapeHDiv<TSHAPE> shapehdiv;
    TPZShapeData data;
    shapehdiv.Initialize(nodeids, orders, sideorient, data);
    TPZFMatrix<REAL> phi;
    TPZManVector<REAL,3> xiReal(dim);
    for(int i=0; i<dim; i++) xiReal[i]=qsifad[i].val();
    shapehdiv.Shape(xiReal, data, phi, divphimaster);
    TPZFMatrix<Fad<REAL>> dphiFAD(2,data.fDPhi.Cols()),phiFAD;
    shapehdiv.Shape(qsifad, data, phiFAD, dphiFAD);
    TPZFMatrix<Fad<REAL>> dirHDiv;
    auto detjac = DetJac(gradX);
    gradX.MultAdd(phiFAD, phiFAD, dirHDiv, 1./detjac, 0.);
    divphicomputed.Resize(dirHDiv.Cols(),1);
    for (int i = 0; i<divphicomputed.Rows(); i++) {
        divphicomputed(i,0) = 0.;
        for (int d=0; d<dim; d++) {
            divphicomputed(i,0) += dirHDiv(d,i).dx(d)*detjac.val();
        }
    }

}
}
