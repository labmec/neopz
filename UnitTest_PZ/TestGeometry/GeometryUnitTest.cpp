/**
 * @file GeometryUnitTest.cpp
 * @brief Define a Unit Test using Boost for all kind of geometries available in the library
 *
 */


#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "pzgmesh.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "pzquad.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"

#include "tpzquadraticcube.h"
#include "tpzquadratictetra.h"
#include "tpzquadraticprism.h"
#include "tpzquadraticpyramid.h"

#include "tpzarc3d.h"
#include "tpzellipse3d.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"
#include "TPZCylinderMap.h"
#include "tpzchangeel.h"
#include "TPZCurve.h"
#include "TPZSurface.h"

#include "TPZVTKGeoMesh.h"


#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.testgeom");
#endif

#include "fad.h"

#include <catch2/catch.hpp>
//#define NOISY //outputs x and grad comparisons
//#define NOISYVTK //prints all elements in .vtk format

std::string dirname = PZSOURCEDIR;
using namespace pzgeom;

/** 
 * @name Generate a geometric mesh with all topology of elements
 *
 * @{ 
 */

template<class T>
void AddElement(TPZGeoMesh &mesh, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
{
    int matid = mesh.NElements()+1;
    T::InsertExampleElement(mesh, matid, lowercorner, size);
    lowercorner[0] += size[0];
}

void FillGeometricMesh(TPZGeoMesh &mesh)
{
    TPZManVector<REAL,3> lowercorner(3,0.),size(3,1.); // Setting the first corner as the origin and the max element size is 1.0;

    AddElement<TPZGeoPoint>(mesh,lowercorner,size); // @omar:: It makes no sense to test gradx of a 0D element
                                                    // @phil : It helps to check whether vectors of proper dimension are created
    AddElement<TPZGeoLinear>(mesh,lowercorner,size);
    AddElement<TPZGeoTriangle>(mesh,lowercorner,size);
    AddElement<TPZGeoQuad>(mesh,lowercorner,size);
    AddElement<TPZGeoCube>(mesh,lowercorner,size);
    AddElement<TPZGeoTetrahedra>(mesh,lowercorner,size);
    AddElement<TPZGeoPrism>(mesh,lowercorner,size);
    AddElement<TPZGeoPyramid>(mesh,lowercorner,size);
    lowercorner[0] = 1.;
    lowercorner[1] = 2.;
    AddElement<TPZQuadraticLine>(mesh,lowercorner,size);
    AddElement<TPZQuadraticTrig>(mesh,lowercorner,size);
    AddElement<TPZQuadraticQuad>(mesh,lowercorner,size);
    AddElement<TPZQuadraticCube>(mesh,lowercorner,size);
    AddElement<TPZQuadraticTetra>(mesh,lowercorner,size);
    AddElement<TPZQuadraticPrism>(mesh,lowercorner,size);
    AddElement<TPZQuadraticPyramid>(mesh,lowercorner,size);
    lowercorner[0] = 1.;
    lowercorner[1] = 3.;
    AddElement<TPZGeoBlend<TPZGeoLinear> >(mesh,lowercorner,size);
    AddElement<TPZGeoBlend<TPZGeoTriangle> >(mesh,lowercorner,size);
    AddElement<TPZGeoBlend<TPZGeoQuad> >(mesh,lowercorner,size);
    AddElement<TPZGeoBlend<TPZGeoCube> >(mesh,lowercorner,size);
    AddElement<TPZGeoBlend<TPZGeoTetrahedra> >(mesh,lowercorner,size);
    AddElement<TPZGeoBlend<TPZGeoPrism> >(mesh,lowercorner,size);
    AddElement<TPZGeoBlend<TPZGeoPyramid> >(mesh,lowercorner,size);
    
    AddElement<TPZArc3D >(mesh,lowercorner,size);
    AddElement<TPZEllipse3D >(mesh,lowercorner,size);

    AddElement<TPZCylinderMap<TPZGeoTriangle>>(mesh,lowercorner,size);
    AddElement<TPZCylinderMap<TPZGeoQuad>>(mesh,lowercorner,size);
//    AddElement<TPZWavyLine>(mesh,lowercorner,size);
//    AddElement<TPZQuadTorus>(mesh,lowercorner,size);
//    AddElement<TPZTriangleTorus>(mesh,lowercorner,size);
//    AddElement<TPZQuadSphere<> >(mesh,lowercorner,size);
//    AddElement<TPZTriangleSphere<> >(mesh,lowercorner,size);
    mesh.BuildConnectivity();
}

void PlotRefinedMesh(TPZGeoMesh &gmesh,const std::string &filename)
{
    gRefDBase.InitializeAllUniformRefPatterns();
    int numref = 4;
    for (int iref=0; iref<numref; iref++) {
        int64_t nel = gmesh.NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh.Element(el);
            if (gel->HasSubElement()) {
                continue;
            }
            TPZStack<TPZGeoEl *> subels;
            gel->Divide(subels);
        }
    }
    std::ofstream out(filename);
    TPZVTKGeoMesh::PrintGMeshVTK(&gmesh, out);
}

/* @} */

/** 
 * @name Testing the conformity of gradx
 * @{ 
 */


/** @} */

TEMPLATE_TEST_CASE("TPZCylinder3D","[special_maps][geometry_tests]",
          pzgeom::TPZGeoTriangle,
          pzgeom::TPZGeoQuad) {

    constexpr REAL tol = 1e-10;

    auto CreateAndTestCylinder = [](auto x_pts,
                                    const TPZVec<REAL> &origin,
                                    const TPZVec<REAL> &axis,
                                    const REAL radius){
        using TGeo=TestType;
        TPZGeoMesh gmesh;
        TPZManVector<int64_t, TGeo::NNodes> nodeind(TGeo::NNodes,-1);
        TPZManVector<REAL,3> x(3,0.);
        for(auto in = 0; in < TestType::NNodes; in++){
            x_pts(in,x);
            nodeind[in] = gmesh.NodeVec().AllocateNewElement(); 
            gmesh.NodeVec()[nodeind[in]].Initialize(x,gmesh);
        }

        constexpr int matid{1};
        auto *gel = new TPZGeoElRefLess<TPZCylinderMap<TGeo>> (nodeind,matid,gmesh);

        gel->Geom().SetOrigin(origin, radius);
        gel->Geom().SetCylinderAxis(axis);
        gel->Geom().ComputeCornerCoordinates(gmesh);
        //let us test the corner nodes
        for(int i = 0; i < TGeo::NNodes; i++){
            TPZManVector<REAL,TGeo::Dimension> qsi(TGeo::Dimension,0);
            TGeo::ParametricDomainNodeCoord(i,qsi);

            TPZManVector<REAL,3> x_calc(3,0.), x_given(3,0.);
            x_pts(i,x_given);
            gel->X(qsi,x_calc);
            CAPTURE(x_calc);
            CAPTURE(x_given);
            REAL diff{0};
            for(auto ix = 0; ix < 3; ix++){
                diff += (x_given[ix] - x_calc[ix])*(x_given[ix] - x_calc[ix]);
            }
            REQUIRE((diff <= tol));
        }
    };
    
    const TPZVec<REAL> origin = {0,0,0};    
    SECTION("+z-oriented axis"){
        const TPZVec<REAL> axis = {0,0,1};
        constexpr REAL r{1};
        auto x_pts = [](int i, TPZVec<REAL> &x){
            switch(i){
            case 0:
                x = {1,0,0};
                break;
            case 1:
                x = {0,1,0};
                break;
            case 2:
                x = {0,1,1};
                break;
            case 3://only for quads
                x = {1,0,1};
                break;
            default:
                PZError<<__PRETTY_FUNCTION__
                       <<"\n invalid number of nodes!\n";
                DebugStop();
            }
        };
        CreateAndTestCylinder(x_pts,origin,axis,r);
    }
    SECTION("-z-oriented axis"){
        const TPZVec<REAL> axis = {0,0,-1};
        constexpr REAL r{1};
        auto x_pts = [](int i, TPZVec<REAL> &x){
            switch(i){
            case 0:
                x = {1,0,0};
                break;
            case 1:
                x = {0,1,0};
                break;
            case 2:
                x = {0,1,1};
                break;
            case 3://only for quads
                x = {1,0,1};
                break;
            default:
                PZError<<__PRETTY_FUNCTION__
                       <<"\n invalid number of nodes!\n";
                DebugStop();
            }
        };
        CreateAndTestCylinder(x_pts,origin,axis,r);
    }
    SECTION("x-oriented axis"){
        const TPZVec<REAL> axis = {1,0,0};
        constexpr REAL r{1};
        auto x_pts = [](int i, TPZVec<REAL> &x){
            switch(i){
            case 0:
                x = {0,1,0};
                break;
            case 1:
                x = {1,1,0};
                break;
            case 2:
                x = {1,0,1};
                break;
            case 3://only for quads
                x = {0,0,1};
                break;
            default:
                PZError<<__PRETTY_FUNCTION__
                       <<"\n invalid number of nodes!\n";
                DebugStop();
            }
        };
        CreateAndTestCylinder(x_pts,origin,axis,r);
    }
    SECTION("bigger radius"){
        const TPZVec<REAL> axis = {0,0,1};
        constexpr REAL r{3};
        auto x_pts = [](int i, TPZVec<REAL> &x){
            switch(i){
            case 0:
                x = {3,0,0};
                break;
            case 1:
                x = {0,3,0};
                break;
            case 2:
                x = {0,3,3};
                break;
            case 3://only for quads
                x = {3,0,3};
                break;
            default:
                PZError<<__PRETTY_FUNCTION__
                       <<"\n invalid number of nodes!\n";
                DebugStop();
            }
        };
        CreateAndTestCylinder(x_pts,origin,axis,r);
    }       
}

TEST_CASE("gradx_tests","[geometry_tests]") {
    TPZGeoMesh gmesh;
    FillGeometricMesh(gmesh);
    
    int npoints = 10;
    REAL tol;
    ZeroTolerance(tol);
    tol *= 20.;
    
//    std::ofstream file("gmesh.txt");
//    gmesh.Print(file);
    
    int nel = gmesh.NElements();
    for(int iel = 0; iel < nel; iel++){
        TPZGeoEl *gel = gmesh.Element(iel);
        int iel_dim = gel->Dimension();
		TPZManVector< REAL, 3 > qsi_r(iel_dim,0);
		TPZManVector<Fad<REAL>, 3 > qsi(iel_dim,0);

        for(int ip = 0; ip < npoints; ip++){
            TPZManVector<REAL,3> pt(iel_dim);
            gel->RandomPoint(pt);
            for(int i = 0; i < iel_dim; i++)
            {
                Fad<REAL> a(iel_dim,i,pt[i]);
                qsi[i] = a;
                qsi_r[i] = a.val();
            }
            
            // FAD
            TPZManVector<Fad<REAL> ,3> x(3,0);
            TPZFMatrix< Fad<REAL> > gradx;

            // REAL
            TPZManVector< REAL ,3> x_r(3,0);
            TPZFNMatrix< 9, REAL > gradx_r;
            gel->X(qsi, x);
            gel->GradX(qsi, gradx);
            
            gel->X(qsi_r, x_r);
            gel->GradX(qsi_r, gradx_r);
            
            int r = gradx_r.Rows();
            int c = gradx_r.Cols();
            
            for(int i = 0; i < r; i++ ){
#ifdef NOISY
                std::cout << " x = " << x_r[i] << std::endl;
                std::cout << " x fad = " << x[i] << std::endl;
#endif
                for(int j = 0; j < c; j++ ){
#ifdef NOISY
                    std::cout << " gradx = " << gradx_r(i,j) << std::endl;
                    std::cout << " gradx fad = " << x[i].dx(j) << std::endl;
#endif
                    REAL diff1 = gradx_r(i,j)-x[i].dx(j);
                    REAL diff2 = gradx_r(i,j)-gradx(i,j).val();
#ifdef REALfloat
                    bool gradx_from_x_fad_check = std::abs(diff1) < tol;
                    bool gradx_vs_gradx_fad_check = std::abs(diff2) < tol;
#else
                    bool gradx_from_x_fad_check = fabs(diff1) < tol;
                    bool gradx_vs_gradx_fad_check = fabs(diff2) < tol;
#endif
                    if(!gradx_from_x_fad_check || ! gradx_vs_gradx_fad_check)
                    {
                        std::cout << "gel type name " << gel->TypeName() << std::endl;
                    }
                    REQUIRE(gradx_from_x_fad_check);
                    REQUIRE(gradx_vs_gradx_fad_check);
                    
                }
            }
        }

        
    }
#ifdef NOISYVTK
    PlotRefinedMesh(gmesh,"AllElements.vtk");
#endif
    
    return;

}

TEST_CASE("changeel_tests","[geometry_tests]") {
    TPZGeoMesh gmesh;

    
    //radius of the circumference
    constexpr REAL radius{1.0};
    //center of the circumference
    TPZVec<REAL> xc = {0,0,0};
    //axis of the cylinder (must be unitary)
    TPZVec<REAL> axis = {0,0,1};
    //this string will be filled at each section to generate .vtk files
    std::string section_name = "";
    
    //simple lambda for creating new node and returning its index
    auto CreateNewNode = [&gmesh] (const TPZVec<REAL> &x){
        const int64_t nodeidx = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeidx].Initialize(x,gmesh);
        return nodeidx;
    };


    //lambda to check whether pts lie in arc/cylinder

    auto TestPts = [radius,xc,axis](TPZGeoEl *gel){
        constexpr REAL tol = std::numeric_limits<REAL>::epsilon();
        const int nsides = gel->NSides();
        const int dim = gel->Dimension();
        constexpr int p{4};
        auto intRule = gel->CreateSideIntegrationRule(nsides-1, p);
        const int npts = intRule->NPoints();
        
        TPZVec<REAL> xi(dim,0.), x(3,0.), dist(3,0.), dist_cross(3,0.);
        for(int ipt = 0; ipt < npts; ipt++){
            REAL w{0};
            intRule->Point(ipt, xi, w);
            gel->X(xi,x);
            TPZManVector<REAL,3> dist = {x[0]-xc[0],x[1]-xc[1],x[2]-xc[2]};
            //cross product
            dist_cross = {
                axis[1]*dist[2] - axis[2]*dist[1],
                axis[2]*dist[0] - axis[0]*dist[2],
                axis[0]*dist[1] - axis[1]*dist[0]
            };

            REAL norm = sqrt(dist_cross[0]*dist_cross[0] +
                             dist_cross[1]*dist_cross[1] +
                             dist_cross[2]*dist_cross[2]);
            REQUIRE(norm== Approx(radius).epsilon(tol));
        }
    };


    constexpr int line_mat{1};
    constexpr int trig_mat{1};
    
    //create linear element that will be replaced by an TPZArc3D
    TPZVec<int64_t> nodevec = {CreateNewNode({1,0,0}),CreateNewNode({0,1,0})};
    
    TPZGeoEl* lin_el =
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodevec,line_mat,gmesh);

    gmesh.BuildConnectivity();
    SECTION("Arc3D"){
        lin_el = TPZChangeEl::ChangeToArc3D(&gmesh, lin_el->Index(), xc, radius);
        REQUIRE(lin_el);
        TestPts(lin_el);
        section_name = "Arc3D";
    }

    SECTION("Blend"){
        lin_el = TPZChangeEl::ChangeToArc3D(&gmesh, lin_el->Index(), xc, radius);
        nodevec.Resize(3);
        nodevec[2] = CreateNewNode({0,0,0});
        TPZGeoEl* trig_el =
            new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodevec,trig_mat,gmesh);
        gmesh.BuildConnectivity();
        trig_el = TPZChangeEl::ChangeToGeoBlend(&gmesh, trig_el->Index());
        REQUIRE(trig_el);
        section_name = "Blend";
    }

    SECTION("Quad"){
        lin_el = TPZChangeEl::ChangeToArc3D(&gmesh, lin_el->Index(), xc, radius);
        nodevec.Resize(3);
        nodevec[2] = CreateNewNode({0,0,0});
        TPZGeoEl* trig_el =
            new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodevec,trig_mat,gmesh);

        gmesh.BuildConnectivity();
        trig_el = TPZChangeEl::ChangeToQuadratic(&gmesh, trig_el->Index());
        REQUIRE(trig_el);
        section_name = "Quad";
    }

    SECTION("Cylinder"){
        lin_el = TPZChangeEl::ChangeToArc3D(&gmesh, lin_el->Index(), xc, radius);
        nodevec.Resize(3);
        nodevec[2] = CreateNewNode({1,0,1});
        TPZGeoEl* trig_el =
            new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodevec,trig_mat,gmesh);

        gmesh.BuildConnectivity();
        trig_el = TPZChangeEl::ChangeToCylinder(&gmesh, trig_el->Index(),
                                                xc, axis, radius);
        REQUIRE(trig_el);
        TestPts(trig_el);
        section_name = "Cylinder";
    }
#ifdef NOISYVTK
    if(section_name.size() != 0){
        PlotRefinedMesh(gmesh, "change_el"+section_name+".vtk");
    }
#endif
}