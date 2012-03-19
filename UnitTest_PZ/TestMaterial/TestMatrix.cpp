/**
* @file
* @brief Contains Unit Tests for methods of the material classes.
*/

//#include "pzmatrix.h"
#include "pzelast3d.h"
#include "iostream"
#include "fstream"


#ifdef USING_BOOST

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN pz material tests

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#endif

/**
 * @brief Create the stiffness matrix of a cube from -1 and 1 in cartesian coordinates, with elastic material
 * @param
 * @note The stiffness is created with the following basis function: (x,0,0), (0,x,0), (0,0,x), (y,0,0), (0,y,0), (0,0,y), (z,0,0), (0,z,0), (0,0,z) 
 * 
 */
TPZFMatrix<REAL> CubeStiffness()
{
	int nummat = 1;
	REAL Ela = 1000, poisson = 0.2;
	TPZVec<REAL> force(3,0.);
	TPZElasticity3D elast3d(nummat, Ela, poisson, force);
	TPZMaterialData cubodata;
	TPZVec <REAL> pt(3,0.);
	TPZFMatrix<REAL> phi(3,1,0.);
	TPZFMatrix<REAL> dphi(3,3,0.);
	dphi(0,0) = 1;
	dphi(1,1) = 1;
	dphi(2,2) = 1;
	cubodata.phi = phi;
	cubodata.dphix = dphi;
	cubodata.x = pt;
	TPZFMatrix<REAL> ek(9,9,0), ef(9,1,0);
	REAL weight = 8.;
	elast3d.Contribute(cubodata, weight, ek, ef);
	return ek;
}

TPZFMatrix<REAL> ReadStiff(std::string &FileName)
{
	std::ifstream in(FileName.c_str());
	char buf[1024];
	in.getline(buf,1024);
	TPZFMatrix<REAL> RightStiff(9,9,0.);
	REAL temp;
	for (int i = 0 ; i < 9 ; i++) 
	{
		for (int j = 0 ; j < 9 ; j++) 
		{
			in >> temp; 
			RightStiff(i,j) = temp;
		}
	}
	return RightStiff;
}

#ifdef USING_BOOST

BOOST_AUTO_TEST_SUITE(material_tests)

BOOST_AUTO_TEST_CASE(test_tonto) 
{
	int a = 3;
	BOOST_CHECK_EQUAL(a,3);
}



BOOST_AUTO_TEST_CASE(test_matriz_rigidez_cubo)
{
	std::string name = "CubeStiffMatrix.txt";
	TPZFMatrix<REAL> RightStiff(ReadStiff(name)), stiff(CubeStiffness());
	REAL tol = 1.e-8;
	bool sym = stiff.VerifySymmetry(tol);
	std::cout << sym << std::endl;
	BOOST_CHECK_EQUAL(sym,1);		// Verify the symmetry of the stiffness matrix
	REAL dif;
	for (int i = 0 ; i < 9 ; i++) 
	{
		for (int j = 0 ; j < 9 ; j++) 
		{
			dif = fabs(stiff(i,j) - RightStiff(i,j));
			BOOST_CHECK_SMALL(dif, (REAL)0.01L);
		}
	}
	//ek.Print("ek: ");	
}

BOOST_AUTO_TEST_SUITE_END()

#endif