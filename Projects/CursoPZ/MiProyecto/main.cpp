#include "pzlog.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzmat2dlin.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzfunction.h"

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include "pzshapelinear.h"

#include "TPZRefPatternTools.h"

#include <time.h>
#include <stdio.h>
#include <fstream>

using namespace std;
using namespace pzshape;

int anothertests = 0;
char saida[514];
int materialId = 4;

std::string Archivo = PZSOURCEDIR;

TPZGeoMesh *CreateGeoMesh(std::string &nome);
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<TPZManVector<REAL> > &points,REAL &distance,bool &isdefined);

// bi-dimensional problem for elasticity
int main_first() {

#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for quadrilaterals and triangles
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	gRefDBase.InitializeUniformRefPattern(ETriangle);

	Archivo += "/Projects/CursoPZ/MiProyecto/";
	Archivo += "MiPlaca.dump";
	TPZGeoMesh *gmesh;
	time_t sttime;
	time_t endtime;

	int nelem;
	REAL distance = 0.;
	bool isdefined = false;
	for(int ii=0;ii<2;ii++) {
		time (& sttime);
		// Creating geometric mesh
		gmesh = CreateGeoMesh(Archivo);
		if(!ii) {			
			// Refinando nas esquinas desejadas
			nelem=0;
			int nrefs = 5;
			TPZManVector<REAL> point(3,0.);
			TPZVec<TPZManVector<REAL> > points(3);
			points[0] = point;
			point[1] = -1.;
			points[1] = point;
			point[0] = 1.;
			points[2] = point;
			
			for(int i=0;i<nrefs;i++) {
				distance = 1./((i+1)*13);
				RefineGeoElements(2,gmesh,points,distance,isdefined);
			}
			// Constructing connectivities
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
		else
			// Refinamento uniforme para toda a malla
			UniformRefine(gmesh,2);
	
		// Creating computational mesh (approximation space and materials)
		int p;
		if(!ii) p = 5;
		else p = 2;
		TPZCompEl::SetgOrder(p);
		TPZCompMesh *cmesh = CreateMesh(gmesh);
		// Colocando a menor ordem para elementos subdivididos
		nelem = 0;
		while(!ii && nelem < cmesh->NElements()-1) {
			TPZCompEl *cel = cmesh->ElementVec()[nelem++];
			if(cel && cel->Reference()->Father()) {
				if(cel->Reference()->Father()->Father())
					((TPZInterpolatedElement*)cel)->PRefine(1);
				((TPZInterpolatedElement*)cel)->PRefine(3);
			}
		}
		cmesh->AutoBuild();
		cmesh->AdjustBoundaryElements();
		cmesh->CleanUpUnconnectedNodes();
		
		// Solving linear equations
		// Initial steps
		TPZAnalysis an (cmesh);
		TPZSkylineStructMatrix strskyl(cmesh);
		an.SetStructuralMatrix(strskyl);
		// Solver (is your choose) 
		TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
		direct->SetDirect(ECholesky);
		an.SetSolver(*direct);
		delete direct;
		direct = 0;

		/*
		 // Caso no simetrico
		 //	TPZFStructMatrix full(cmesh);
		 TPZBandStructMatrix full(cmesh);
		 an.SetStructuralMatrix(full);
		 an.Solution().Zero();
		 TPZStepSolver<REAL> step;
		 step.SetDirect(ELU);
		 an.SetSolver(step);
		 */
		an.Run();
		
		// Calculando o tempo que demorou para calcular em cada cenario 
		time (& endtime);
		int time_elapsed = endtime - sttime;
		std::cout << "\n\n\tHP-Adaptive Methods....step: " << ii+1 << " time elapsed " << time_elapsed << "\n\n\n";

		// Post processing
		TPZStack<std::string> scalarnames, vecnames;
		std::string filename;
		if(!ii) filename = "ElasticitySolutions.vtk";
		else filename += "ElasticitySolutions1.vtk";
		scalarnames.Push("POrder");
		scalarnames.Push("SigmaX");
		scalarnames.Push("SigmaY");
		scalarnames.Push("Pressure");
		scalarnames.Push("MaxStress");
		scalarnames.Push("TauXY");
		vecnames.Push("displacement");
		vecnames.Push("PrincipalStress1");
		vecnames.Push("PrincipalStress2");
		//vecnames.Push("POrder");
		an.DefineGraphMesh(2,scalarnames,vecnames,filename);
		
		an.PostProcess(0);
	}
	return 0;
}

int main() {
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for quadrilaterals and triangles
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	gRefDBase.InitializeUniformRefPattern(ETriangle);
	
	time_t sttime;
	time_t endtime;
	
    // First rectangular mesh:
	// The rectangular mesh has four corners: (0,-1,0), (1,-1,0), (1,0,0) and (0,0,0)
	// and was divides in two segments on X and two on Y, then hx = 0.5 and hy = 0.5
	// Has 4 elements, 9 connects
	cout << "Generating geometric mesh bi-dimensional ...\n";
	TPZManVector<REAL> point(3,0.), pointlast(3,0.);
	TPZGeoMesh* gmesh1 = new TPZGeoMesh;
	TPZManVector<REAL> x0(3,0.), x1(3,0.);  // Corners of the rectangular mesh. Coordinates of the first extreme are zeros.
	x0[1] = -1.; x1[0] = 1.;
	TPZManVector<int> nx(2,2);   // subdivisions in X and in Y. 
	TPZGenGrid gen1(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1  
	gen1.SetElementType(0);       // type = 0 means rectangular elements
	gen1.Read(gmesh1,materialId);             // generating grid in gmesh
	
	// Selecting base functions on vertices
	if(anothertests) {
		// Setting Chebyshev polynomials as orthogonal sequence generating shape functions
		TPZShapeLinear::fOrthogonal = &TPZShapeLinear::Legendre;
		sprintf(saida,"meshextrudedLeg.vtk");
		
	}
	else {
		sprintf(saida,"meshextrudedTChe.vtk");
	}
	
	int nelem;
	REAL radius = 0.;
	bool isdefined = false;
	
	for(int ii=0;ii<2;ii++) {
		// Constructing a geometric mesh
		TPZGeoMesh* gmesh = new TPZGeoMesh;
		x0[0] = -1.; x0[1] = 0.;
		x1[0] = 1.; x1[1] = 1.;
		nx[0] = 4; //nx[1] *= 2;
		TPZGenGrid gen(nx,x0,x1);
		gen.SetElementType(0);
		gen.ReadAndMergeGeoMesh(gmesh,gmesh1,materialId);
		// Inserting boundary elements with associated material
		// Bottom is fixed
		point[0] = 0.; point[1] = -1;
		pointlast[0] = 1.; pointlast[1] = -1.;
		gen.SetBC(gmesh,point,pointlast,1);
		// Top boundary has vertical force applied
		point[0] = -1; point[1] = 1.;
		pointlast[0] = 1.; pointlast[1] = 1.;
		gen.SetBC(gmesh,point,pointlast,2);
		// Vertical right boundary has horizontal force applied to left
		point[0] = 1; point[1] = -1.;
		pointlast[0] = 1.; pointlast[1] = 1.;
		gen.SetBC(gmesh,point,pointlast,3);

		// Initializing the process
		time (& sttime);
		if(!ii) {			
			// Refinando nas esquinas desejadas
			nelem=0;
			int nrefs = 1;
			point[0] = point[1] = point[2] = 0.;
			TPZVec<TPZManVector<REAL> > points(3);
			points[0] = point;
			point[1] = -1.;
			points[1] = point;
			point[0] = 1.;
			points[2] = point;
			
			for(int i=0;i<nrefs;i++) {
				RefineGeoElements(2,gmesh,points,radius,isdefined);
				radius *= .5;
			}
			// Constructing connectivities
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
		}
		else {
			// Refinamento uniforme para toda a malla
			UniformRefine(gmesh,1);
		}
	
		// Creating computational mesh (approximation space and materials)
		int p;
		if(!ii) p = 7;
		else p = 3;
		TPZCompEl::SetgOrder(p);
		TPZCompMesh *cmesh = CreateMesh(gmesh);
		// Disminuindo a ordem p dos elementos subdivididos
		// A cada nivel disminue em uma unidade o p, mas não será menor de 1.
		nelem = 0;
		TPZGeoEl *gelem;
		while(!ii && nelem < cmesh->NElements()-1) {
			REAL pCopy = p;
			TPZCompEl *cel = cmesh->ElementVec()[nelem++];
			if(cel) {
				gelem = cel->Reference();
				while(gelem) {
					gelem = gelem->Father();
					if(gelem) {
						if(pCopy != 1) pCopy--;
						((TPZInterpolatedElement*)cel)->PRefine(pCopy);
					}
				}
			}
		}
		cmesh->AutoBuild();
		cmesh->AdjustBoundaryElements();
		cmesh->CleanUpUnconnectedNodes();
		
		// Solving linear equations
		// Initial steps
		TPZAnalysis an (cmesh);
		TPZSkylineStructMatrix strskyl(cmesh);
		an.SetStructuralMatrix(strskyl);
		// Solver (is your choose) 
		TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
		direct->SetDirect(ECholesky);
		an.SetSolver(*direct);
		delete direct;
		direct = 0;
		
		an.Run();
		
		// Calculando o tempo que demorou para calcular em cada cenario 
		time (& endtime);
		int time_elapsed = endtime - sttime;
		std::cout << "\n\n\tHP-Adaptive Methods....step: " << ii+1 << " time elapsed " << time_elapsed << "\n\n\n";
		
		// Post processing
		TPZStack<std::string> scalarnames, vecnames;
		std::string filename;
		if(!ii) filename = "ElasticitySolutions.vtk";
		else filename += "ElasticitySolutions1.vtk";
		scalarnames.Push("POrder");
		scalarnames.Push("SigmaX");
		scalarnames.Push("SigmaY");
		scalarnames.Push("Pressure");
		scalarnames.Push("MaxStress");
		scalarnames.Push("TauXY");
		vecnames.Push("displacement");
		vecnames.Push("PrincipalStress1");
		vecnames.Push("PrincipalStress2");
		//vecnames.Push("POrder");
		an.DefineGraphMesh(2,scalarnames,vecnames,filename);
		
		an.PostProcess(1);
		
		delete cmesh;
		delete gmesh;
	}
	return 0;
}

void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<TPZManVector<REAL> > &points,REAL &distance,bool &isdefined) {
	TPZManVector<REAL> centerpsi(3), center(3);
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;

	int nelem = 0;
	int ngelem=gmesh->NElements();
	int i, npoints = points.NElements();
	// na esquina inferior esquerda Nó = (0,-1,0)
	while(nelem<ngelem) {
		gel = gmesh->ElementVec()[nelem++];
		if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
		gel->CenterPoint(8,centerpsi);
		gel->X(centerpsi,center);
		if(!isdefined) {
			TPZVec<REAL> FirstNode(3,0.);
			gel->CenterPoint(0,centerpsi);
			gel->X(centerpsi,FirstNode);
			distance = 1.1*TPZGeoEl::Distance(center,FirstNode);
			isdefined = true;
		}
		for(i=0;i<npoints;i++) {
			REAL semidiag = TPZGeoEl::Distance(center,points[i]);
			if(semidiag < distance) {
				gel->Divide(sub);
				break;
			}
		}
	}
}
//*******Shell to deforming************
TPZGeoMesh *CreateGeoMesh(std::string &archivo) {
	
	// Ejemplo uni-dimensional para la generacion de una malla para un reservatorio 
	TPZReadGIDGrid grid;
	TPZGeoMesh *meshgrid = grid.GeometricGIDMesh(archivo);
	if(!meshgrid->NElements())
		return 0;
	
	return meshgrid;
}

//*************************************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh) {
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();
	
    // Creating elasticity material
    TPZMaterial * mat = new TPZElasticityMaterial(4,2000.,0.3,0.,0.);
    cmesh->InsertMaterialObject(mat);
	
	// Creating four boundary condition
    TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);
	TPZMaterial *bcBottom, *bcRigth, *bcTop, *bcLeft;
	
	// Condicion livre - nada para hacer
    bcLeft = mat->CreateBC(mat,5,3,val1,val2);
    cmesh->InsertMaterialObject(bcLeft);
	// Condicion de Dirichlet fijando la posicion de la placa
	val1(1,1) = 1000000.;
    bcBottom = mat->CreateBC(mat,1,0,val1,val2);
	cmesh->InsertMaterialObject(bcBottom);
	// Condicion de aplicar una fuerza horizontal
	val1(1,1) = 0.;
	val2(1,0) = 10.;
    bcTop = mat->CreateBC(mat,2,1,val1,val2);
	cmesh->InsertMaterialObject(bcTop);
	// Aplicando fuerza zero
	val2(1,0) = 0.;
    bcRigth = mat->CreateBC(mat,3,1,val1,val2);
	cmesh->InsertMaterialObject(bcRigth);
	
	// Inserting boundary conditions into computational mesh
	
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

void UniformRefine(TPZGeoMesh* gmesh, int nDiv)
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
	// Re-constructing connectivities
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}
