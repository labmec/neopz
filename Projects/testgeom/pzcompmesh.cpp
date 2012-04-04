#include "pzgbmesh.h"
#include "TPZGeoElement.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzshapepoint.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"

#include "TPZRefPatternDataBase.h"
#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"

#include "pzartdiff.h"

#include "pzeuler.h"
#include "pzeulerconslaw.h"

using namespace std;

/** Make a uniform refinement of the gmesh, nDiv times */
void UniformRefine(TPZAutoPointer<TPZGeoMesh> gmesh, int nDiv);

int main() {
    
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	char expression[260];
    // First rectangular mesh
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;

	//TPZManVector<int> nx(6,2);   // subdivisions in X and in Y -> Then all the intervals are of  the 0.1 cm.
	TPZManVector<int> nx(2,1);   // subdivisions in X and in Y -> Then all the intervals are of  the 0.3 cm in X and 0.2 cm in Y.
	TPZManVector<REAL> x0(3,0.), x1(3,.6);  // Corners of the rectangular mesh
	x1[1] = 0.2;
	x1[2] = 0.;

	TPZGenGrid gen(nx,x0,x1);    // mesh generator 
	gen.SetElementType(0);       // type = 0 means rectangular elements
	gen.Read(gmesh);            // generating mesh in gmesh

	ofstream saida("malhateste.txt");
	strncpy(expression,"Malha inicial",strlen("Malha inicial")+1);
	gen.Print(expression,saida);
	
	// Second rectangular domain - subdividions and corners of the second rectangular mesh
    TPZAutoPointer<TPZGeoMesh> gmesh2 = new TPZGeoMesh;
	nx[0] = 10;   // subdivision x in 30 intervals
	nx[1] = 4;    // subdivision y in 8 intervals -> Then all the intervals are of the 0.1 cm. 
	x0[1] = .2;
	x1[0] = 3.;
	x1[1] = 1.;

	TPZGenGrid gen2(nx,x0,x1);   // second mesh generator
	gen2.SetElementType(0);      // type = 0 means rectangular elements, type = 1 means triangular elements
	// generating gmesh2 on data of the gen2 and merge gmesh into the gmesh2
	gen2.ReadAndMergeGeoMesh(gmesh2,gmesh);
	
	// setting bc condition -1 [no flux - is wall] from (0.,0.) until (2.,1.)
	x0[1] = 0.;
	x1[1] = .2;
	gen2.SetBC(gmesh2,x0,x1,-1);
	// setting bc condition -1 from (2.,3.) until (0.,3.)
	x0[0] = 3.;
	x0[1] = 1.;
	x1[0] = 0.;
	x1[1] = 1.;
	gen2.SetBC(gmesh2,x0,x1,-1);
	// setting bc condition -2 from (3.,.2) until (3.,1.)
	x1[0] = 3.;
	x1[1] = .2;
	gen2.SetBC(gmesh2,x1,x0,-2);
	// setting bc condition -3 from (0.,0.) until (0.,1.)
	x0[0] = 0.;
	x1[0] = x1[1] = 0.;
	gen2.SetBC(gmesh2, x0, x1, -3);
	
	// Uniform refinement of the geometrical mesh, two level
	int nDiv = 0;
	UniformRefine(gmesh2, nDiv);
	gmesh2->Print(saida);

	// Creating computational mesh associated with geometric mesh gmesh2
	// First we need to create material object
	int p = 2;   // interpolation order
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh2);
	cmesh->SetDimModel(2);
	REAL deltat = 0.01;
	REAL gamma = 1.4;
	TPZAutoPointer<TPZMaterial> mat = new TPZEulerConsLaw(1,deltat,gamma,2,SUPG_AD);
	cmesh->InsertMaterialObject(mat);
	// creating boundary condition object
//	cmesh->InsertMaterialObject((TPZMaterial *)(mat->CreateBC(mat,-1,TPZEulerEquation::EFreeSlip,val1,val2)));
//	cmesh->InsertMaterialObject((TPZMaterial *)(mat->CreateBC(mat,-2,TPZEulerEquation::EFreeSlip,val1,val2)));
//	cmesh->InsertMaterialObject((TPZMaterial *)(mat->CreateBC(mat,-3,TPZEulerEquation::EFreeSlip,val1,val2)));
	// Inserts Dirichlet boundary condition - wall then zero values for all variables
	int nstates = mat->NStateVariables();
	TPZFMatrix<REAL> val1(nstates,nstates,0.),val2(nstates,1,0.);
	cmesh->InsertMaterialObject((TPZMaterial *)(mat->CreateBC(mat,-1,0,val1,val2)));   // Dirichlet condition -  (wall)
	// Inserts Neumann boundary condition - input flux
	cmesh->InsertMaterialObject((TPZMaterial *)(mat->CreateBC(mat,-2,1,val1,val2)));   // Neumann condition - fluxo livre
	// Inserts Neumann boundary condition - output free flux
	cmesh->InsertMaterialObject((TPZMaterial *)(mat->CreateBC(mat,-3,1,val1,val2)));   // Neumann condition - fluxo livre
	
	// Making all the computational elements as discontinuous
//	TPZCreateApproximationSpace::SetAllCreateFunctionsDiscontinuous();
	cmesh->SetAllCreateFunctionsContinuous();
	cmesh->SetDefaultOrder(p);
	TPZCompElDisc::SetgOrder(p);

	// creating a computational elements and the degree of freedom nodes
	cmesh->AutoBuild();
	cmesh->Print(saida);
	
	/** PART 2 : Constructing three-dimensional mesh for testing */
	
	TPZGeoMesh *gmesh3;
	// Using ExtendGridDimension from Pre module
	TPZExtendGridDimension extendmesh(gmesh2 ,0.1);
	gmesh3 = extendmesh.ExtendedMesh();
	gmesh3->Print(saida);
	
	// constructing computational mesh
	TPZCompMesh *cmesh3 = new TPZCompMesh(gmesh3);
	cmesh3->SetDimModel(3);
	TPZAutoPointer<TPZMaterial> mat3 = new TPZEulerEquation(1,gamma);
	cmesh3->InsertMaterialObject(mat3);
	// Inserts Dirichlet boundary condition - wall then zero values for all variables
	int nstates3 = mat3->NStateVariables();
	TPZFMatrix<REAL> val_1(nstates3,nstates3,0.),val_2(nstates3,1,0.);
	cmesh3->InsertMaterialObject((TPZMaterial *)(mat3->CreateBC(mat3,-1,0,val_1,val_2)));   // Dirichlet condition -  (wall)
	// Inserts Neumann boundary condition - input flux
	cmesh3->InsertMaterialObject((TPZMaterial *)(mat3->CreateBC(mat3,-2,1,val_1,val_2)));   // Neumann condition - fluxo livre
	// Inserts Neumann boundary condition - output free flux
	cmesh3->InsertMaterialObject((TPZMaterial *)(mat3->CreateBC(mat3,-3,1,val_1,val_2)));   // Neumann condition - fluxo livre
	
	// Making all the computational elements as discontinuous
	//	TPZCreateApproximationSpace::SetAllCreateFunctionsDiscontinuous();
	cmesh3->SetDefaultOrder(p);
	//TPZCompElDisc::SetgOrder(p);
	
	// creating a computational elements and the degree of freedom nodes
	cmesh3->AutoBuild();
	cmesh3->Print(saida);
	
	saida.close();
	return 0;

}

void UniformRefine(TPZAutoPointer<TPZGeoMesh> gmesh, int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
}

