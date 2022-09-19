/**
* @file
* @brief Contains Unit Tests for methods of the material classes.
*/

//#include "pzmatrix.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZHybridDarcyFlow.h"
#include "TPZShapeH1.h"
#include "pzshapequad.h"
#include "TPZShapeDisc.h"
#include "iostream"
#include "fstream"
#include <catch2/catch.hpp>

/**
 * @brief Create the stiffness matrix of a quadrilaterial in master element space
 * @param
 *
 */
static void computeStiffnessH1(TPZFMatrix<STATE> &ek)
{
    ek.Zero();
    TPZMaterialDataT<STATE> DataVec[3];
    TPZManVector<int64_t> ids = {0,1,2,3};
    TPZManVector<int> connectorders = {1,1,1,1,1};
    
    TPZShapeH1<pzshape::TPZShapeQuad>::Initialize(ids,connectorders,DataVec[0]);
    int64_t nshape = DataVec[0].fPhi.Rows();
//    int order = 0;
//    int dimension = 2;
//    pzshape::TPZShapeDisc::MShapeType shtype = pzshape::TPZShapeDisc::ETensorial;
//    pzshape::TPZShapeDisc::Initialize(dimension,order,shtype,DataVec[1]);
    TPZDarcyFlow material(1, 2);
    pzshape::TPZShapeQuad::IntruleType integrule(2);
    ek.Redim(nshape,nshape);
    TPZFMatrix<STATE> ef(nshape,1,0.);
    int npoints = integrule.NPoints();
    TPZManVector<REAL> pt(2);
    for(int p = 0; p< npoints; p++)
    {
        REAL weight;
        integrule.Point(p, pt, weight);
        TPZShapeH1<pzshape::TPZShapeQuad>::Shape(pt, DataVec[0]);
        DataVec[0].dphix = DataVec[0].fDPhi;
        material.Contribute(DataVec[0], weight, ek, ef);
    }
}

static void computeStiffnessHybrid(TPZFMatrix<STATE> &ek)
{
    ek.Zero();
    TPZVec<TPZMaterialDataT<STATE>> DataVec(3);
    TPZManVector<int64_t> ids = {0,1,2,3};
    TPZManVector<int> connectorders = {1,1,1,1,1};
    
    TPZShapeH1<pzshape::TPZShapeQuad>::Initialize(ids,connectorders,DataVec[0]);
    int64_t nshape = DataVec[0].fPhi.Rows();
    int order = 0;
    int dimension = 2;
    pzshape::TPZShapeDisc::MShapeType shtype = pzshape::TPZShapeDisc::ETensorial;
    pzshape::TPZShapeDisc::Initialize(order,dimension,shtype,DataVec[1]);
    pzshape::TPZShapeDisc::Initialize(order,dimension,shtype,DataVec[2]);
    TPZHybridDarcyFlow material(1, 2);
    pzshape::TPZShapeQuad::IntruleType integrule(2);
    ek.Redim(nshape+2,nshape+2);
    TPZFMatrix<STATE> ef(nshape+2,1,0.);
    int npoints = integrule.NPoints();
    TPZManVector<REAL> pt(dimension);
    for(int p = 0; p< npoints; p++)
    {
        REAL weight;
        integrule.Point(p, pt, weight);
        TPZShapeH1<pzshape::TPZShapeQuad>::Shape(pt, DataVec[0]);
        int degree = 0;
        pzshape::TPZShapeDisc::Shape(dimension,degree,pt,DataVec[1].fPhi,DataVec[1].fDPhi,shtype);
        pzshape::TPZShapeDisc::Shape(dimension,degree,pt,DataVec[2].fPhi,DataVec[2].fDPhi,shtype);
        DataVec[0].dphix = DataVec[0].fDPhi;
        material.Contribute(DataVec, weight, ek, ef);
    }
}



TEST_CASE("test_matriz_darcy","[material_tests]")
{
    TPZFNMatrix<16,STATE> stiff;
    computeStiffnessH1(stiff);
	REAL tol = 1.e-8;
	bool sym = stiff.VerifySymmetry(tol);
	REQUIRE(sym==1);		// Verify the symmetry of the stiffness matrix
    REAL RightStiff[4][4] = {
        {2./3.,-1./6,-1./3.,-1./6.},
        {-1./6.,2./3.,-1./6,-1./3.},
        {-1./3.,-1./6.,2./3.,-1./6},
        {-1./6,-1./3.,-1./6.,2./3.}
    };
	REAL dif;
	for (int i = 0 ; i < 4 ; i++)
	{
		for (int j = 0 ; j < 4 ; j++)
		{
			dif = fabs(stiff(i,j) - RightStiff[i][j]);
			REQUIRE(dif == Approx(0.0).margin(0.00001));
		}
	}
}

TEST_CASE("test_matriz_hybriddarcy","[material_tests]")
{
    TPZFNMatrix<36,STATE> stiff;
    computeStiffnessHybrid(stiff);
    REAL tol = 1.e-8;
    bool sym = stiff.VerifySymmetry(tol);
    REQUIRE(sym==1);        // Verify the symmetry of the stiffness matrix
    REAL RightStiff[6][6] = {
        {2./3.,-1./6,-1./3.,-1./6.,1,0},
        {-1./6.,2./3.,-1./6,-1./3.,1,0},
        {-1./3.,-1./6.,2./3.,-1./6,1,0},
        {-1./6,-1./3.,-1./6.,2./3.,1,0},
        {1,1,1,1,0,-4},
        {0,0,0,0,-4,0}
    };
    REAL dif;
    for (int i = 0 ; i < 6 ; i++)
    {
        for (int j = 0 ; j < 6 ; j++)
        {
            dif = fabs(stiff(i,j) - RightStiff[i][j]);
            REQUIRE(dif == Approx(0.0).margin(0.00001));
        }
    }
}
