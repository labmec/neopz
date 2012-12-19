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
#include <cmath>

using namespace std;
using namespace pzshape;

int anothertests = 0;
char saida[514];
int materialId = 4;

REAL ElasticityModulus = 1000000.;

TPZManVector<REAL> ErrorsTrue(100,0.);
TPZManVector<REAL> Errors(100,0.);

std::string Archivo = PZSOURCEDIR;

TPZGeoMesh *CreateGeoMesh(std::string &nome);
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,bool hasforcingfunction);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<TPZManVector<REAL> > &points,REAL &distance,bool &isdefined);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL> &points,REAL r,REAL &distance,bool &isdefined);

void RightTerm(const TPZVec<REAL> &x, TPZVec<REAL> &force);

void ExactSol(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);

void GetPointsOnCircunference(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL> > &Points);

void ComputeDisplacementError(REAL &error,REAL &errorL2,TPZCompMesh *cmesh);

/** Laplace equation on square - Volker John article 2000 */
int main() {
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for quadrilaterals and triangles
	gRefDBase.InitializeUniformRefPattern(EOned);
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	gRefDBase.InitializeUniformRefPattern(ETriangle);
	
	time_t sttime;
	time_t endtime;
	
	ofstream fileerrors("ErrorsHPProcess.txt");
	// To compute the errors
	TPZVec<REAL> ervec(100,0.0);
	TPZVec<REAL> ervecL2(100,0.0);

	int ii, nrefs = 10;
	for(ii=8;ii<nrefs;ii++) {
		// First rectangular mesh:
		// The rectangular mesh has four corners: (0,-1,0), (1,-1,0), (1,0,0) and (0,0,0)
		// and was divides in two segments on X and two on Y, then hx = 0.5 and hy = 0.5
		// Has 4 elements, 9 connects
		cout << "Generating geometric mesh bi-dimensional ...\n";
		TPZManVector<REAL> point(3,0.), pointlast(3,0.);
		TPZGeoMesh* gmesh = new TPZGeoMesh;
		TPZManVector<REAL> x0(3,0.), x1(3,1.);  // Corners of the rectangular mesh. Coordinates of the first extreme are zeros.
		x1[2] = 0.;
		TPZManVector<int> nx(8,8);   // subdivisions in X and in Y. 
		TPZGenGrid gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1  
		gen.SetElementType(0);       // type = 0 means rectangular elements
		gen.Read(gmesh,materialId);             // generating grid in gmesh
		
		// Inserting boundary elements with associated material
		// Bottom is fixed
		point[0] = 0.; point[1] = 0.;
		pointlast[0] = 1.; pointlast[1] = 0.;
		gen.SetBC(gmesh,point,pointlast,-1);
		// Top boundary has vertical force applied
		point[0] = 1.; point[1] = 0.;
		pointlast[0] = 1.; pointlast[1] = 1.;
		gen.SetBC(gmesh,point,pointlast,-1);
		// Vertical right boundary has horizontal force applied to left
		point[0] = 1.; point[1] = 1.;
		pointlast[0] = 0.; pointlast[1] = 1.;
		gen.SetBC(gmesh,point,pointlast,-1);
		// Vertical right boundary has horizontal force applied to left
		point[0] = 0.; point[1] = 1.;
		pointlast[0] = 0.; pointlast[1] = 0.;
		gen.SetBC(gmesh,point,pointlast,-1);
		
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
		bool isdefined = false;
		
		// Initializing the process
		time (& sttime);
//		if(!ii) {			
			// Refinando nas esquinas desejadas
			nelem=0;
			int npoints = 36;
			point[0] = point[1] = 0.5; point[2] = 0.0;
			REAL r = 0.25, radius = 0.25;
			TPZVec<TPZManVector<REAL> > Points(npoints);
			GetPointsOnCircunference(npoints,point,r,Points);
			
			for(int i=0;i<(ii+1);i++) {
				// To refine elements with center near to points than radius
//				RefineGeoElements(2,gmesh,Points,radius,isdefined);
				// Para refinar elementos con centro tan cerca de la circuferencia cuanto radius 
//				RefineGeoElements(2,gmesh,point,r,radius,isdefined);
				RefineGeoElements(2,gmesh,point,r,radius,isdefined);
				radius *= .6;
			}
			// Constructing connectivities
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
//		}
//		else {
			// Refinamento uniforme para toda a malla
//			UniformRefine(gmesh,4);
//		}
		
		// Creating computational mesh (approximation space and materials)
		int p = 7, pinit;
//		if(!ii) 
//		else p = 2;
		pinit = p;
		TPZCompEl::SetgOrder(p);
		TPZCompMesh *cmesh = CreateMesh(gmesh,true);
		// Disminuindo a ordem p dos elementos subdivididos
		// Primeiro sera calculado o mayor nivel de refinamento
		// A cada nivel disminue em uma unidade o p, mas não será menor de 1.
		int level, highlevel = 0;
		nelem = 0;
		while(nelem < cmesh->NElements()) {
			TPZCompEl *cel = cmesh->ElementVec()[nelem++];
			if(cel) {
				level = cel->Reference()->Level();
			}
			if(level > highlevel)
				highlevel = level;
		}
		nelem = 0;
		while(highlevel && nelem < cmesh->NElements()) {
			TPZCompEl *cel = cmesh->ElementVec()[nelem++];
			if(cel) {
				level = cel->Reference()->Level();
				if(level == highlevel)
					((TPZInterpolatedElement*)cel)->PRefine(1);
				else if(level == 0)
					((TPZInterpolatedElement*)cel)->PRefine(p);
				else {
					REAL porder = (p/highlevel);
					if(porder < 1)
						((TPZInterpolatedElement*)cel)->PRefine(1);
					else
						((TPZInterpolatedElement*)cel)->PRefine((int)(porder*(highlevel-level)));
				}
			}
		}
		cmesh->AutoBuild();
		cmesh->AdjustBoundaryElements();
		cmesh->CleanUpUnconnectedNodes();
		
		// Solving linear equations
		// Initial steps
		TPZAnalysis an (cmesh);
		an.SetExact(ExactSol);
		
		TPZSkylineStructMatrix strskyl(cmesh);
		an.SetStructuralMatrix(strskyl);
		// Solver (is your choose) 
		TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
		direct->SetDirect(ECholesky);
		an.SetSolver(*direct);
		delete direct;
		direct = 0;
		
		// Solving
		an.Run();

		// Computing error
		ComputeDisplacementError(ervec[ii],ervecL2[ii],cmesh);
		
		/*		TPZAdaptMesh adapt;
		adapt.SetCompMesh (cmesh);
		adapt.GetAdaptedMesh(valerror,valtruerror,ervec,ExactSol,truervec,effect,0);
		ErrorsTrue[ii-1] = valtruerror;
		Errors[ii-1] = valerror;
*/		
		// Calculando o tempo que demorou para calcular em cada cenario 
		time (& endtime);
		char pp[3];
		sprintf(pp,"%d",pinit);
		int time_elapsed = endtime - sttime;
		std::cout << "\n\n\tHP-Adaptive Methods....step: " << ii+1 << " time elapsed " << time_elapsed << "\n\n\n";
		
		// Post processing
		TPZStack<std::string> scalarnames, vecnames;
		std::string filename = "ElastSolutions";
		filename += "_p";
		filename += pp;
		filename += "_hL";
		sprintf(pp,"%d",ii);
		filename += pp;
		filename += ".vtk";
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
		
		delete cmesh;
		delete gmesh;
	}
	
	// Printing computed errors
	fileerrors << "Approximation Error: " << std::endl;
	for(ii=0;ii<nrefs;ii++)
		fileerrors << "Refinement " << ii << "  Error " << ervec[ii] << "  ErrorL2 " << ervecL2[ii] << std::endl;
	fileerrors << std::endl << std::endl;
	fileerrors.close();
	
	return 0;
}
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL> &point,REAL r,REAL &distance,bool &isdefined) {
	TPZManVector<REAL> centerpsi(3), center(3);
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	
	int nelem = 0;
	int ngelem=gmesh->NElements();
	// na esquina inferior esquerda Nó = (0,-1,0)
	while(nelem<ngelem) {
		gel = gmesh->ElementVec()[nelem++];
		if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
		gel->CenterPoint(gel->NSides()-1,centerpsi);
		gel->X(centerpsi,center);
		if(!isdefined) {
			TPZVec<REAL> FirstNode(3,0.);
			gel->CenterPoint(0,centerpsi);
			gel->X(centerpsi,FirstNode);
			REAL distancia = TPZGeoEl::Distance(center,FirstNode);
			if(distancia > distance) distance = distancia;
			isdefined = true;
		}
		REAL centerdist = TPZGeoEl::Distance(center,point);
		if(fabs(r-centerdist) < distance) {
			gel->Divide(sub);
		}
	}
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
		gel->CenterPoint(gel->NSides()-1,centerpsi);
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
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,bool hasforcingfunction) {
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();
	
    // Creating elasticity material
    TPZMaterial * mat = new TPZElasticityMaterial(4,ElasticityModulus,0.3,0.,0.);
	if(hasforcingfunction)
		mat->SetForcingFunction(new TPZDummyFunction<STATE>(RightTerm));
    cmesh->InsertMaterialObject(mat);
	
	// Creating four boundary condition
    TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);
	TPZMaterial *bcBottom, *bcRigth, *bcTop, *bcLeft;
	
	// Condicion livre - nada para hacer
    bcLeft = mat->CreateBC(mat,-5,0,val1,val2);
    cmesh->InsertMaterialObject(bcLeft);
	// Condicion de Dirichlet fijando la posicion de la placa
	if(!hasforcingfunction) val1(1,1) = 1000000.;
    bcBottom = mat->CreateBC(mat,-1,0,val1,val2);
	cmesh->InsertMaterialObject(bcBottom);
	// Condicion de aplicar una fuerza horizontal
	if(!hasforcingfunction) {
		val1(1,1) = 0.;
		val2(1,0) = 10.;
		bcTop = mat->CreateBC(mat,-2,1,val1,val2);
	}
	else 
		bcTop = mat->CreateBC(mat,-2,0,val1,val2);
	cmesh->InsertMaterialObject(bcTop);
	// Aplicando fuerza zero
	if(!hasforcingfunction) {
		val2(1,0) = 0.;
		bcRigth = mat->CreateBC(mat,-3,1,val1,val2);
	}
	else 
		bcRigth = mat->CreateBC(mat,-3,0,val1,val2);
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

void RightTerm(const TPZVec<REAL> &x, TPZVec<REAL> &force) {
//	REAL Epsilon = 1000000;
	REAL B = 16./M_PI;
	REAL F = 2*sqrt(ElasticityModulus);
	REAL G = -0.4375;

	REAL sum = x[0]*(x[0]-1) + x[1]*(x[1]-1);
	REAL prod = x[0]*(x[0]-1)*x[1]*(x[1]-1);

	REAL temp = F*(G-sum);
	REAL arctan = atan(temp);
	REAL den = (1+temp*temp)*(1+temp*temp);
	REAL num = 2*F*(sum*(2*F*F*prod*(8*G+1)-(1+F*F*G*G)+F*F*sum*(2*G-6*prod-sum))-2*prod*(F*F*G+5*F*F*G*G+5));

	force[0] = B*(sum*(M_PI+2*arctan)+(num/den));
	force[0] *= (-ElasticityModulus);
	
	force[1] = force[0];
}

void ExactSol(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
	REAL F = 2*sqrt(ElasticityModulus);
	REAL arc = F*((0.25*0.25) - (x[0] - 0.5)*(x[0] - 0.5) - (x[1] - 0.5)*(x[1] - 0.5));
	REAL prodx = x[0]*(x[0]-1.);
	REAL prody = x[1]*(x[1]-1.);
	REAL prod = prodx*prody;
	sol[0] = 8*prod*(1+(2./M_PI)*(atan(arc)));
	sol[1] = sol[0];
	REAL temp = prody*(2*x[0]-1.)*(M_PI + 2*atan(arc));
	REAL frac = 2*prod*F*(1.-2*x[0]);
	frac = frac/(1+arc*arc);
	dsol(0,0) = (8./M_PI)*(temp + frac);
	dsol(0,1) = dsol(0,0);
	temp = prodx*(2*x[1]-1.)*(M_PI + 2*atan(arc));
	frac = 2*prod*F*(1.-2*x[1]);
	frac = frac/(1+arc*arc);
	dsol(1,0) = (8./ M_PI)*(temp + frac);
	dsol(1,1) = dsol(1,0);
}

void GetPointsOnCircunference(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL> > &Points) {
	Points.Resize(npoints);
	TPZManVector<REAL> point(3,0.);
	REAL angle = (2*M_PI)/npoints;
	for(int i=0;i<npoints;i++) {
		point[0] = center[0]+radius*cos(i*angle);
		point[1] = center[1]+radius*sin(i*angle);
		Points[i] = point;
	}
}

void ComputeDisplacementError(REAL &error,REAL &errorL2,TPZCompMesh *cmesh) {
	int i, it, nelem = cmesh->NElements();

	TPZBlock<REAL> flux;
	
	int orderp = 4;
	TPZIntQuad ordem2dq(orderp,orderp);
	int npoints = ordem2dq.NPoints();
	TPZVec<REAL> point(3,0.), x(3,0.);
	REAL weight = 0.;

	TPZCompEl *cel;
	int var = 9;
	REAL DispMagnitude = 0.;
	TPZVec<REAL> SolCel(5,0.);
	TPZFMatrix<REAL> DSolCel(5,5);
	
	REAL errorLoc = 0., errorLocL1 = 0.;
	for(i=0;i<nelem;i++) {
		cel = cmesh->ElementVec()[i];
		if(!cel || cel->Dimension() != 2) continue;
//		for(i=0;i<error.NElements();i++)
//			errorLoc[i] = 0.0;
//		cel->EvaluateError(ExactSol,errorLoc,&flux);
//		for(i=0;i<error.NElements();i++)
//			error[i] += errorLoc[i];
		for(it=0;it<npoints;it++){
			ordem2dq.Point(it,point,weight);
			cel->Reference()->X(point,x);
			cel->Solution(point,var,SolCel);
			DispMagnitude = sqrt(SolCel[0]*SolCel[0]+SolCel[1]*SolCel[1]);
			ExactSol(x,SolCel,DSolCel);
			REAL temp = sqrt(SolCel[0]*SolCel[0]+SolCel[1]*SolCel[1]);
			errorLocL1 += weight * fabs(DispMagnitude - temp);
			errorLoc += weight * (DispMagnitude - temp)*(DispMagnitude - temp);
		}
		errorLocL1 *= cel->Reference()->Volume();
		errorLoc *= cel->Reference()->Volume();
	}
	error = errorLocL1;
	errorL2 = sqrt(errorLoc);
}

/////

/////

// bi-dimensional problem for elasticity on square domain
int main_GID() {
	
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
		TPZCompMesh *cmesh = CreateMesh(gmesh,false);
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

/** Laplace equation on L-domain */
int main_LDomain() {
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing uniform refinements for quadrilaterals and triangles
	//gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
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
			int nrefs = 3;
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
		TPZCompMesh *cmesh = CreateMesh(gmesh,false);
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
