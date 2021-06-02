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