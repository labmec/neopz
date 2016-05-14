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
#include "TPZCurve.h"
#include "TPZSurface.h"

#include "TPZVTKGeoMesh.h"


#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testgeom"));
#endif

// Using Unit Test of the Boost Library
#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz geometry_tests tests

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"
#include "boost/test/output_test_stream.hpp"

#endif

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

    AddElement<TPZGeoPoint>(mesh,lowercorner,size);
    AddElement<TPZGeoLinear>(mesh,lowercorner,size);
    AddElement<TPZGeoTriangle>(mesh,lowercorner,size);
    AddElement<TPZGeoQuad>(mesh,lowercorner,size);
    AddElement<TPZGeoCube>(mesh,lowercorner,size);
    AddElement<TPZGeoTetrahedra>(mesh,lowercorner,size);
    AddElement<TPZGeoPrism>(mesh,lowercorner,size);
    AddElement<TPZGeoPyramid>(mesh,lowercorner,size);
    lowercorner[0] = 1.;
    lowercorner[1] = 2.;
    AddElement<TPZGeoBlend<TPZGeoLinear> >(mesh,lowercorner,size);
    AddElement<TPZGeoBlend<TPZGeoTriangle> >(mesh,lowercorner,size);
    AddElement<TPZGeoBlend<TPZGeoQuad> >(mesh,lowercorner,size);
    AddElement<TPZGeoBlend<TPZGeoCube> >(mesh,lowercorner,size);
    AddElement<TPZGeoBlend<TPZGeoTetrahedra> >(mesh,lowercorner,size);
    AddElement<TPZGeoBlend<TPZGeoPrism> >(mesh,lowercorner,size);
    AddElement<TPZGeoBlend<TPZGeoPyramid> >(mesh,lowercorner,size);
    lowercorner[0] = 1.;
    lowercorner[1] = 1.;
    AddElement<TPZQuadraticLine>(mesh,lowercorner,size);
    AddElement<TPZQuadraticTrig>(mesh,lowercorner,size);
    AddElement<TPZQuadraticQuad>(mesh,lowercorner,size);
    AddElement<TPZQuadraticCube>(mesh,lowercorner,size);
    AddElement<TPZQuadraticTetra>(mesh,lowercorner,size);
    AddElement<TPZQuadraticPrism>(mesh,lowercorner,size);
    AddElement<TPZQuadraticPyramid>(mesh,lowercorner,size);
    mesh.BuildConnectivity();
}

void PlotRefinedMesh(TPZGeoMesh &gmesh,const std::string &filename)
{
    gRefDBase.InitializeAllUniformRefPatterns();
    int numref = 3;
    for (int iref=0; iref<numref; iref++) {
        long nel = gmesh.NElements();
        for (long el=0; el<nel; el++) {
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


#ifdef USING_BOOST


BOOST_AUTO_TEST_SUITE(geometry_tests)


BOOST_AUTO_TEST_CASE(gradx_tests) {
    
    TPZGeoMesh gmesh;
    FillGeometricMesh(gmesh);
    
    int npoints = 10;
    REAL tol = 1.0e-10;
    TPZManVector< REAL, 3 > qsi_r(3);
    TPZVec<Fad<REAL> > qsi(3);
    
    std::ofstream file("gmesh.txt");
    gmesh.Print(file);
    
    int nel = gmesh.NElements();
    for(int iel = 0; iel < nel; iel++){
        TPZGeoEl *gel = gmesh.Element(iel);
        int iel_dim = gel->Dimension();

        for(int ip = 0; ip < npoints; ip++){
            for(int i = 0; i < iel_dim; i++){
                REAL val = (REAL) rand() / (RAND_MAX);
                Fad<REAL> a(iel_dim,i,val);
                qsi[i] = a;
                qsi_r[i] = a.val();
            }
            
            TPZVec<Fad<REAL> > x(3);
            TPZFMatrix< REAL > gradxr;
            gel->X(qsi, x);
            gel->GradX(qsi_r, gradxr);
            int r = gradxr.Rows();
            int c = gradxr.Cols();
            for(int i = 0; i < r; i++ ){
                for(int j = 0; j < c; j++ ){
                    bool check = fabs(gradxr(i,j)-x[i].dx(j)) < tol;
                    BOOST_CHECK(check);
                    
                }
            }
        }

        
    }
                
    PlotRefinedMesh(gmesh,"AllElements.vtk");                
    
    return;

}





BOOST_AUTO_TEST_SUITE_END()


#endif

