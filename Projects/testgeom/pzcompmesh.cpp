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

#include "pzeuler.h"

using namespace std;

// Program to reproduce the situation as issue 1 into CodeGoogle NeoPZ
#define REFPATTERNDIR "/Users/jorge/Labmec/GoogleCodes/neopz/Refine/RefPatterns"

void UniformRefine(TPZAutoPointer<TPZGeoMesh> gmesh, int nDiv);

int main() {
    
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	char expression[260];
    // First rectangular mesh
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;

	TPZManVector<int> nx(6,2);   // subdivisions in X and in Y -> Then all the intervals are of  the 0.1 cm.
	TPZManVector<REAL> x0(3,0.), x1(3,.6);  // Corners of the rectangular mesh
	x1[1] = 0.2;
	x1[2] = 0.;

	TPZGenGrid gen(nx,x0,x1);    // mesh generator 
	gen.SetElementType(0);       // type = 0 means rectangular elements
	gen.Read(*(gmesh.operator->()));            // generating mesh in gmesh

	ofstream saida("malhateste.txt");
	strncpy(expression,"Malha inicial",strlen("Malha inicial")+1);
	gen.Print(expression,saida);
	
	// Second rectangular domain - subdividions and corners of the second rectangular mesh
    TPZAutoPointer<TPZGeoMesh> gmesh2 = new TPZGeoMesh;
	nx[0] = 30;   // subdivision x in 30 intervals
	nx[1] = 8;    // subdivision y in 8 intervals -> Then all the intervals are of the 0.1 cm. 
	x0[1] = .2;
	x1[0] = 3.;
	x1[1] = 1.;

	TPZGenGrid gen2(nx,x0,x1);   // second mesh generator
	gen2.SetElementType(0);
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
	TPZAutoPointer<TPZMaterial> mat = new TPZEulerEquation(1,1.4);
	cmesh->InsertMaterialObject(mat);
	TPZFMatrix val1,val2;
	// creating boundary condition object
	cmesh->InsertMaterialObject((TPZMaterial *)(mat->CreateBC(mat,-1,TPZEulerEquation::EFreeSlip,val1,val2)));
	cmesh->InsertMaterialObject((TPZMaterial *)(mat->CreateBC(mat,-2,TPZEulerEquation::EFreeSlip,val1,val2)));
	cmesh->InsertMaterialObject((TPZMaterial *)(mat->CreateBC(mat,-3,TPZEulerEquation::EFreeSlip,val1,val2)));
	
	// Making all the computational elements as discontinuous
	TPZCompMesh::SetAllCreateFunctionsDiscontinuous();
//	TPZCompMesh::SetAllCreateFunctionsContinuous();
	cmesh->SetDefaultOrder(p);
	TPZCompElDisc::SetgOrder(p);

	// creating a computational elements and the degree of freedom nodes
	cmesh->AutoBuild();
	cmesh->Print(saida);
	
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

