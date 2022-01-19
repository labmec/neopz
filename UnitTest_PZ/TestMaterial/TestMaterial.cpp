/**
* @file
* @brief Contains Unit Tests for methods of the material classes.
*/

//#include "pzmatrix.h"
#include "Elasticity/TPZElasticity3D.h"
#include "iostream"
#include "fstream"
#include <catch2/catch.hpp>

/**
 * @brief Create the stiffness matrix of a cube from -1 and 1 in cartesian coordinates, with elastic material
 * @param
 * @note The stiffness is created with the following basis function: (x,0,0), (0,x,0), (0,0,x), (y,0,0), (0,y,0), (0,0,y), (z,0,0), (0,z,0), (0,0,z) 
 * 
 */
TPZFMatrix<STATE> computeStressStrain()
{
	int nummat = 1;
	REAL Ela = 1000, poisson = 0.2;
	TPZVec<STATE> force(3,0.);
	TPZElasticity3D elast3d(nummat, Ela, poisson, force);
	TPZMaterialDataT<STATE> cubodata;
	TPZVec <REAL> pt(3,0.);
	TPZFMatrix<REAL> phi(3,1,0.);
	TPZFMatrix<REAL> dphi(3,3,0.);
	dphi(0,0) = 1;
	dphi(1,1) = 1;
	dphi(2,2) = 1;
	cubodata.fPhi = phi;
	cubodata.dphix = dphi;
	cubodata.x = pt;
	cubodata.axes.Redim(3, 3);
	cubodata.axes.Identity();
	TPZFMatrix<STATE> ek(9,9,0), ef(9,1,0);
	REAL weight = 8.;
	elast3d.Contribute(cubodata, weight, ek, ef);
	return ek;
}

TPZFMatrix<STATE> readStressStrain(std::string &FileName)
{
	std::ifstream in(FileName.c_str());
	char buf[1024];
	in.getline(buf,1024);
	TPZFMatrix<STATE> RightStiff(9,9,0.);
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

TEST_CASE("test_tonto","[material_tests]") 
{
	int a = 3;
	REQUIRE(a==3);
}



TEST_CASE("test_matriz_rigidez_cubo","[material_tests]")
{
	std::string name = "CubeStiffMatrix.txt";
	TPZFMatrix<STATE> RightStiff(readStressStrain(name)), stiff(computeStressStrain());
	REAL tol = 1.e-8;
	bool sym = stiff.VerifySymmetry(tol);
	std::cout << sym << std::endl;
	REQUIRE(sym==1);		// Verify the symmetry of the stiffness matrix
	REAL dif;
	for (int i = 0 ; i < 9 ; i++) 
	{
		for (int j = 0 ; j < 9 ; j++) 
		{
			dif = fabs(stiff(i,j) - RightStiff(i,j));
			REQUIRE(dif == Approx(0.0).margin(0.01));
		}
	}

}