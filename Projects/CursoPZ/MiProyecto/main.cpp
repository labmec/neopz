/**
 * @file
 * @brief Projeto elaborado para Minicurso na Pontificia Catolica del Peru - Lima - 2012/07
 */

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

#include "TPZParSkylineStructMatrix.h"

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
// Crea malla computacional sem forcingfunction quando hasforcingfunction = 0, ou toma diferentes forcingfuncition para diferentes
// valores de hasforcingfunction
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int hasforcingfunction);

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<TPZManVector<REAL> > &points,REAL &distance,bool &isdefined);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL> &points,REAL r,REAL &distance,bool &isdefined);

void RightTermCircle(const TPZVec<REAL> &x, TPZVec<REAL> &force);
void RightTermProduct(const TPZVec<REAL> &x, TPZVec<REAL> &force);

void ExactSolNull(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
	sol.Fill(0.0);
	dsol.Zero();
}
void ExactSolCircle(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);
void ExactSolProduct(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);

void GetPointsOnCircunference(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL> > &Points);

void ComputeDisplacementError(REAL &error,REAL &errorL2,TPZCompMesh *cmesh);

void formatTimeInSec(char *strtime,int timeinsec);

void GradientReconstructionByLeastSquares(TPZFMatrix<REAL> &gradients,TPZCompMesh *cmesh,int var,int n_var=0,bool continuous=false);

int problem = 1;

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
	int time_elapsed;
	char tempo[256];
	
	ofstream fileerrors("ErrorsHPProcess.txt");
	// To compute the errors
	TPZVec<REAL> ervec(100,0.0);
	TPZVec<REAL> ervecL2(100,0.0);
	
	// Printing computed errors
	fileerrors << "Approximation Error: " << std::endl;
	
	int i, ii, nrefs = 8;
	int jj, nthreads = 9;
	for(jj=2;jj<nthreads;jj++) {
		for(ii=2;ii<nrefs;ii++) {
			// Initializing the generation mesh process
			time (& sttime);
			// First rectangular mesh:
			// The rectangular mesh has four corners: (0,-1,0), (1,-1,0), (1,0,0) and (0,0,0)
			// and was divides in two segments on X and two on Y, then hx = 0.5 and hy = 0.5
			// Has 4 elements, 9 connects
			cout << "\nGenerating geometric mesh bi-dimensional ...\n" << "\tRefinement: " << ii << "   Threads: " << jj;
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
			
			// Refinando no local desejado
			nelem=0;
			int npoints = 360;
			point[0] = point[1] = 0.5; point[2] = 0.0;
			REAL r = 0.25, radius = 0.15;
			TPZVec<TPZManVector<REAL> > Points(npoints);
			GetPointsOnCircunference(npoints,point,r,Points);
			
			for(i=0;i<ii;i+=2) {
				// To refine elements with center near to points than radius
				//				RefineGeoElements(2,gmesh,Points,radius,isdefined);
				// Para refinar elementos con centro tan cerca de la circuferencia cuanto radius 
				RefineGeoElements(2,gmesh,point,r,radius,isdefined);
				RefineGeoElements(2,gmesh,point,r,radius,isdefined);
				radius *= 0.5;
			}
			if(i==ii) {
				RefineGeoElements(2,gmesh,point,r,radius,isdefined);
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
			int p = 5, pinit;
			pinit = p;
			TPZCompEl::SetgOrder(p);
			TPZCompMesh *cmesh = CreateMesh(gmesh,problem);
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
			
			// closed generation mesh process
			time (& endtime);
			time_elapsed = endtime - sttime;
			time_elapsed = endtime - sttime;
			formatTimeInSec(tempo, time_elapsed);
			std::cout << "  Time elapsed " << time_elapsed << " <-> " << tempo << "\n\n";

			// Solving linear equations
			// Initial steps
			std::cout << "Solving HP-Adaptive Methods....step: " << ii << "  Threads " << jj << "\n";
			TPZAnalysis an (cmesh);
			if(!problem)
				an.SetExact(ExactSolNull);
			else if(problem==1)
				an.SetExact(ExactSolCircle);
			else if(problem==2)
				an.SetExact(ExactSolProduct);
			
			TPZParSkylineStructMatrix strskyl(cmesh,jj);
			an.SetStructuralMatrix(strskyl);
			// Solver (is your choose) 
			TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
			direct->SetDirect(ECholesky);
			an.SetSolver(*direct);
			delete direct;
			direct = 0;
			
			// Initializing the solving process
			time (& sttime);
			// Solving
			an.Run();
			
			// Calculando o tempo que demorou para calcular em cada cenario 
			time (& endtime);
			time_elapsed = endtime - sttime;
			formatTimeInSec(tempo, time_elapsed);

			char pp[3];
			sprintf(pp,"%d",pinit);
			std::cout << "\t....step: " << ii << "  Threads " << jj << "  Time elapsed " << time_elapsed << " <-> " << tempo << "\n\n\n";
			
			// Computing error
			ComputeDisplacementError(ervec[ii],ervecL2[ii],cmesh);
			fileerrors << "Refinement: " << ii << "  Threads: " << jj << "  NEquations: " << cmesh->NEquations() << "  ErrorL1: " << ervec[ii] << "  ErrorL2: " 
			<< ervecL2[ii] << "  TimeElapsed: " << time_elapsed << " <-> " << tempo << std::endl;
			
			TPZVec<REAL> Solution(cmesh->NEquations(),0.);
			
			// Computing approximation of gradient
			/** 
			 * @brief Method to reconstruct a gradient after run Solve of the analysis
			 * @param cmesh Computational mesh with solution */
//			TPZFMatrix<REAL> gradients;
//			GradientReconstructionByLeastSquares(gradients,cmesh,0,0,0);
//			gradients.Print();
//			cout << std::endl << std::endl;
//			GradientReconstructionByLeastSquares(gradients,cmesh,0,0,1);
//			gradients.Print();

			// Post processing
			TPZStack<std::string> scalarnames, vecnames;
			std::string filename = "ElastSolutions";
			filename += "_p";
			filename += pp;
			filename += "_hL";
			sprintf(pp,"%d",ii);
			filename += pp;
			sprintf(pp,"_PAR%d",jj);
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
		break;
	}
	
	fileerrors << std::endl << std::endl;
	fileerrors.close();
}

void formatTimeInSec(char *strtime,int timeinsec) {
	if(!strtime) return;
	memset(strtime,0,strlen(strtime));
	//	strtime[0] = '\0';
	int anos=0, meses=0, dias=0, horas=0, minutos=0, segundos=0;
	while(1) {
		if(timeinsec < 60) {
			segundos = timeinsec;
			break;
		}
		else {
			timeinsec -= 60;
			minutos++;
			if(minutos > 59) {
				minutos -= 60;
				horas++;
				if(horas > 23) {
					horas -= 24;
					dias++;
					if(dias > 29) {
						dias -= 30;
						meses++;
						if(meses > 11) {
							meses -= 12;
							anos++;
						}
					}
				}
			}
		}
	}
	// Formating
	if(anos)
		sprintf(strtime,"%d a, %d m, %d d, %d:%d:%d",anos,meses,dias,horas,minutos,segundos);
	else {
		if(meses) 
			sprintf(strtime,"%d m, %d d, %d:%d:%d",meses,dias,horas,minutos,segundos);
		else {
			if(dias)
				sprintf(strtime,"%d d, %d:%d:%d",dias,horas,minutos,segundos);
			else
				sprintf(strtime,"%d:%d:%d",horas,minutos,segundos);
		}
	}
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
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int hasforcingfunction) {
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();
	
    // Creating elasticity material
    TPZMaterial * mat = new TPZElasticityMaterial(4,ElasticityModulus,0.3,0.,0.);
	switch(hasforcingfunction) {
		case 1:
			mat->SetForcingFunction(new TPZDummyFunction<STATE>(RightTermCircle));
			break;
		case 2:
			mat->SetForcingFunction(new TPZDummyFunction<STATE>(RightTermProduct));
			break;
		default:
			break;
	}
    cmesh->InsertMaterialObject(mat);
	// Make compatible dimension of the model and the computational mesh
	cmesh->SetDimModel(mat->Dimension());
	
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

void RightTermCircle(const TPZVec<REAL> &x, TPZVec<REAL> &force) {
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

void RightTermProduct(const TPZVec<REAL> &x, TPZVec<REAL> &force) {
	force[0] = 32.*(x[0]*(x[0]-1.)+x[1]*(x[1]-1.));
	force[1] = 0.;
}
void ExactSolProduct(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
	sol[0] = 16*x[0]*x[1]*(1-x[0])*(1-x[1]);
	sol[1] = 0.;
	dsol(0,0) = 16*(1.-2*x[0])*(1.-x[1])*x[1];
	dsol(0,1) = 0.;
	dsol(1,0) = 16*(1.-2*x[1])*(1.-x[0])*x[0];
	dsol(1,1) = 0.;
}

void ExactSolCircle(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
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
REAL PartialDerivateX(const TPZVec<REAL> &x) {
	REAL F = 2*sqrt(ElasticityModulus);
	REAL arc = F*((0.25*0.25) - (x[0] - 0.5)*(x[0] - 0.5) - (x[1] - 0.5)*(x[1] - 0.5));
	REAL prodx = x[0]*(x[0]-1.);
	REAL prody = x[1]*(x[1]-1.);
	REAL result = (8./M_PI)*prody*(2*x[0]-1);
	REAL temp = M_PI + 2*atan(arc);
	REAL frac = 2*F*prodx;
	frac = frac/(1+arc*arc);
	temp -= frac;
	return (result*temp);
}
REAL PartialDerivateY(const TPZVec<REAL> &x) {
	REAL F = 2*sqrt(ElasticityModulus);
	REAL arc = F*((0.25*0.25) - (x[0] - 0.5)*(x[0] - 0.5) - (x[1] - 0.5)*(x[1] - 0.5));
	REAL prodx = x[0]*(x[0]-1.);
	REAL prody = x[1]*(x[1]-1.);
	REAL result = (8./M_PI)*prodx*(2*x[1]-1);
	REAL temp = M_PI + 2*atan(arc);
	REAL frac = 2*F*prody;
	frac = frac/(1+arc*arc);
	temp -= frac;
	return (result*temp);
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
			if(!problem)
				ExactSolNull(x,SolCel,DSolCel);
			else if(problem==1)
				ExactSolCircle(x,SolCel,DSolCel);
			else if(problem==2)
				ExactSolProduct(x,SolCel,DSolCel);
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

/** Reconstrucción del gradiente utilizando la linearizacion (Taylor) de la solución para los centros de todos los elementos vecinos */
/** Formula: u(xbi,ybi,zbi) = u(xa,ya,za) + a*(xbi-xa) + b*(ybi-ya) + c*(zbi-za)  ->  donde Grad(u) ~= (a,b,c) */
/** (xa,ya,za) es el centro del elemento donde queremos aproximar o gradiente de u */
/** (xbi,ybi,zbi) son los centros de los elementos vecinos al elemento corriente por alguno de sus lados, e enumerados por i */
void GradientReconstructionByLeastSquares(TPZFMatrix<REAL> &gradients,TPZCompMesh *cmesh,int var,int n_var,bool continuous) {
	int i, nstates=0;
	TPZCompEl *cel;
	int dim = cmesh->Dimension();
	for(i=0;i<cmesh->NElements();i++) {
		cel = cmesh->ElementVec()[i];
		if(cel && cel->Dimension() == dim) {
			nstates = cel->Material()->NSolutionVariables(var);
			break;
		}
	}
	
	int nelem = cmesh->NElements();
    gradients.Redim(nelem,4*dim);
	
	int k, side;
	int counter = 0;
	
	TPZStack<TPZCompElSide> neighs;
	int nneighs;
	
	TPZManVector<REAL> normal(3,0.0);
	TPZManVector<REAL> centerpsi(3,0.0);
	TPZManVector<REAL> center(3,0.0), centerbeta(3,0.0);
	TPZManVector<REAL> solalfa(nstates,0.0), solbeta(nstates,0.0);
	
	TPZFMatrix<REAL> A(dim,dim);    // Linear System matrix
	TPZFMatrix<REAL> B(dim,1,0.);   // Linear System vector
	
	// Creando las matrices para aplicar el metodo de los minimos cuadrados
	TPZFMatrix<REAL> DeltaH(nneighs,dim,0.);
	TPZFMatrix<REAL> DeltaHTranspose(dim,nneighs,0.);
	TPZFMatrix<REAL> DifSol(nneighs,1,0.);
	REAL Grad;
	
	// Calculando el gradiente por elemento computacional
	for(i=0;i<nelem;i++) {
		cel = cmesh->ElementVec()[i];
		// Nada sera realizado para elementos con dimension diferente de la dimension del problema
		if(!cel || cel->Dimension()!=dim) continue;
		
		// Limpiando las matrizes
		A.Zero(); B.Zero();
		// Encontramos el centro del elemento corriente cel
		TPZGeoEl* gelalfa = cel->Reference();
		gelalfa->CenterPoint(gelalfa->NSides()-1,centerpsi);
		center.Fill(0.);
		gelalfa->X(centerpsi,center);
		cel->Solution(centerpsi,var,solalfa);
		
		// PREFERENCIAL PARA CASOS DE CALCULO CON FUNCIONES DISCONTINUAS - Pues utiliza los valores de la solución en los elementos vecinos
		if(!continuous) {
			neighs.Resize(0);
			// Procuramos todos los elementos vecinos a cel (sobre todos los lados) sin duplicados
//			for(side = cel->Reference()->NCornerNodes(); side < cel->NConnects(); side++) {
			for(side = 0; side < cel->NConnects(); side++) {
				TPZCompElSide celside(cel,side);
				celside.ConnectedElementList(neighs,1,0);
			}
			nneighs = neighs.NElements();
			// si no hay vecinos continuamos con el siguiente elemento
			if(!nneighs) continue;
			// si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente			
			// Para cada vecino calculamos los deltaH (desde su centro al centro del elemento corriente)
			// y el valor de la solucion en su centro solbeta
			DeltaH.Redim(nneighs,dim);
			DeltaHTranspose.Redim(dim,nneighs);
			DifSol.Redim(nneighs,1);
			// Montando la matriz de los deltas DeltaH y de las diferencias de las soluciones DifSol
			for(int ineighs=0;ineighs<nneighs;ineighs++) {
				TPZGeoEl* gelbeta = neighs[ineighs].Element()->Reference();
				if(!gelbeta)
					DebugStop();
				centerpsi.Fill(0.0);
				centerbeta.Fill(0.0);
				gelbeta->CenterPoint(gelbeta->NSides()-1,centerpsi);
				gelbeta->X(centerpsi,centerbeta);
				gelbeta->Reference()->Solution(centerpsi,var,solbeta);
				for(k=0;k<dim;k++)
					DeltaH(ineighs,k) = centerbeta[k] - center[k];
				DifSol(ineighs,0) = solbeta[n_var] - solalfa[n_var];
			}
		}
		else {
			int nsides = cel->NConnects()-1;
			// Para cada lado calculamos los deltaH (desde el centro del elemento al centro del lado de dimension menor a él
			// y el valor de la solucion en su centro solbeta
			DeltaH.Redim(nsides,dim);
			DeltaHTranspose.Redim(dim,nsides);
			DifSol.Redim(nsides,1);
			// Procuramos todos los puntos medios de cada lado del elemento y calculamos baseados en los valores de la solucion sobre ellos
			for(side = 0; side < nsides; side++) {
				centerpsi.Fill(0.0);
				centerbeta.Fill(0.0);
				cel->Reference()->CenterPoint(side,centerpsi);
				cel->Reference()->X(centerpsi,centerbeta);
				cel->Solution(centerpsi,var,solbeta);
				for(k=0;k<dim;k++)
					DeltaH(side,k) = centerbeta[k] - center[k];
				DifSol(side,0) = solbeta[n_var] - solalfa[n_var];
				
			}
		}
		// Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u) 
		DeltaH.Transpose(&DeltaHTranspose);
		B = DeltaHTranspose*DifSol;
		A = DeltaHTranspose*DeltaH;
		A.SolveDirect(B,ELU);
		
		// Normalizando el vector gradiente
		Grad = 0.0;
		for(k=0;k<dim;k++)
			Grad += (B(k,0)*B(k,0));
		// Almacenando los gradientes encontrados
		for(k=0;k<dim;k++) {
			if(!IsZero(B(k))) {
				gradients(counter,k) = B(k,0)/sqrt(Grad);
			}
			gradients(counter,dim+k) = center[k];
			if(!k) {
				REAL dudx = PartialDerivateX(center);
				REAL dudy = PartialDerivateY(center);
				REAL dist = sqrt(dudx*dudx + dudy*dudy);
				if(!IsZero(dist)) {
					gradients(counter,2*dim) = dudx/dist;
					gradients(counter,2*dim+1) = dudy/dist;
				}
				dist = sqrt((center[0]-0.5)*(center[0]-0.5)+(center[1]-0.5)*(center[1]-0.5));
				if(!IsZero(dist)) {
					gradients(counter,3*dim) = (0.5-center[0])/dist;
					gradients(counter,3*dim+1) = (0.5-center[1])/dist;
				}
			}
		}
		counter++;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////
//  //  FROM ADAPT_HP_JORGE
////////////////////////////////////////////////////////////////////////////////////////////////

#include "pzelast3d.h"

// Global variable
int gLMax;
int NUniformRefs = 2;
REAL alfa = M_PI/6.;
//bool anothertests = false;
int nstate = 2;

/** Printing level */
int gPrintLevel = 0;
//bool gDebug = false;

//TPZFNMatrix<16,REAL> Rot(4,4,0.),RotInv(4,4,0.);

void Exact(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);

void InitializeSolver(TPZAnalysis &an);
void InitialSolutionLinearConvection(TPZFMatrix<REAL> &InitialSol, TPZCompMesh *cmesh);
void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);

void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial=true, const int matidtodivided=1);

/**
 * @brief This project shows the creation of a rectangular mesh (two-dimensional) and the creation of a three-dimensional cube mesh using extrude method (ExtendMesh).
 */
int main_AdaptHP(int argc, char *argv[]) {
	
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	const int L = 4;
	gLMax = L-1;
	char saida[260];
	
	//-----------  INITIALIZING CONSTRUCTION OF THE MESHES
	
	int r, dim;
	
	// Initializing a ref patterns
	//gRefDBase.InitializeAllUniformRefPatterns();
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	
	// output files
    std::ofstream convergence("conv3d.txt");
    std::ofstream out("output.txt");
	
	/** Set polynomial order */
	int p, pmax = 2;
	for(p=1;p<pmax;p++) {
		// First rectangular mesh:
		// The rectangular mesh has four corners: (0,0,0), (1,0,0), (1,1,0) and (0,1,0)
		// and was divides in two segments on X and two on Y, then hx = 0.5 and hy = 0.5
		// Has 4 elements, 9 connects and 8 bc elements
		cout << "Generating geometric mesh bi-dimensional ...\n";
		TPZGeoMesh* gmesh = new TPZGeoMesh;
		TPZManVector<REAL> x0(3,0.), x1(3,1.);  // Corners of the rectangular mesh. Coordinates of the first extreme are zeros.
		x1[2] = 0.;
		TPZManVector<int> nx(3,3);   // subdivisions in X and in Y. 
		TPZGenGrid gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1  
		gen.SetElementType(0);       // type = 0 means rectangular elements
		gen.Read(gmesh);             // generating grid in gmesh
		
		// Applying hp adaptive techniques 2012/10/01
		if(anothertests) {
			// Setting Chebyshev polynomials as orthogonal sequence generating shape functions
			TPZShapeLinear::fOrthogonal = &TPZShapeLinear::Legendre;
			sprintf(saida,"meshextrudedLeg.vtk");
			
		}
		else {
			sprintf(saida,"meshextrudedTChe.vtk");
		}	
		
		// Refinement of the some element	
		TPZGeoEl *gel, *gel1, *gel2, *gel3;
		TPZVec<TPZGeoEl *> sub;
		TPZVec<TPZGeoEl *> subsub;
		gel = gmesh->ElementVec()[4];
		//	gel1 = gmesh->ElementVec()[1];
		//	gel2 = gmesh->ElementVec()[2];
		//	gel3 = gmesh->ElementVec()[3];
		gel->Divide(sub);
		//	sub[0]->Divide(subsub);
		//		sub[1]->Divide(subsub);
		sub[2]->Divide(subsub);
		//	sub[3]->Divide(subsub);
		/*	gel1->Divide(sub);
		 sub[0]->Divide(subsub);
		 sub[1]->Divide(subsub);
		 sub[2]->Divide(subsub);
		 sub[3]->Divide(subsub);
		 gel2->Divide(sub);
		 sub[0]->Divide(subsub);
		 sub[1]->Divide(subsub);
		 sub[2]->Divide(subsub);
		 sub[3]->Divide(subsub);
		 gel3->Divide(sub);
		 sub[0]->Divide(subsub);
		 sub[1]->Divide(subsub);
		 sub[2]->Divide(subsub);
		 sub[3]->Divide(subsub);
		 */	
		// Constructing connectivities
		gmesh->ResetConnectivities();
		gmesh->BuildConnectivity();
		gmesh->Print();
		// Printing COMPLETE initial geometric mesh 
		PrintGeoMeshVTKWithDimensionAsData(gmesh,saida);
		
		TPZCompEl::SetgOrder(p);
		// Creating computational mesh
		TPZCompMesh *comp = new TPZCompMesh(gmesh);
		
		// Creating and inserting materials into computational mesh
		TPZMaterial * mat = new TPZElasticityMaterial(1,1.e5,0.2,.5,0);   // two-dimensional
		comp->InsertMaterialObject(mat);
		dim = mat->Dimension();
		nstate = mat->NStateVariables();
		
		// Boundary conditions
		// Dirichlet
		TPZFMatrix<REAL> val1(3,3,0.),val2(3,1,5.);
		val1(0,0) = 1.;
		TPZMaterial *bnd = mat->CreateBC(mat,-1,0,val1,val2);
		comp->InsertMaterialObject(bnd);
		// Neumann
		val2(0,0)=30.; val2(1,0) = 10.;
		bnd = mat->CreateBC(mat,-2,1,val1,val2);
		comp->InsertMaterialObject(bnd);
		
		// Constructing and adjusting computational mesh
		comp->AutoBuild();
		comp->AdjustBoundaryElements();   // Adjust boundary elements and higher level of refinement, clean elements but not connects into them
		comp->CleanUpUnconnectedNodes();  // Clean connects not connected at least one element enabled.
		comp->Print();
		
		
		//--- END construction of the meshes
		/** Variable names for post processing */
		TPZStack<std::string> scalnames, vecnames;
		if(mat->NSolutionVariables(mat->VariableIndex("POrder")) == 1)
			scalnames.Push("POrder");
		else
			vecnames.Push("POrder");
		if(mat->NSolutionVariables(mat->VariableIndex("Error")) == 1)
			scalnames.Push("Error");
		else
			vecnames.Push("Error");
		if(mat->NSolutionVariables(mat->VariableIndex("state")) == 1)
			scalnames.Push("state");
		else
			vecnames.Push("state");
		
		if(nstate == 1) {
			scalnames.Push("TrueError");
			scalnames.Push("EffectivityIndex");
		}else if(nstate == 2) {
			scalnames.Push("sig_x");
			scalnames.Push("sig_y");
			scalnames.Push("tau_xy");
		}
		if(nstate == 3) {
			scalnames.Push("StressX");
			scalnames.Push("StressY");
			scalnames.Push("StressZ");
			vecnames.Push("PrincipalStress");
			vecnames.Push("PrincipalStrain");
		}
		// END Determining the name of the variables
		
		// INITIAL POINT FOR SOLVING AND APPLYING REFINEMENT
		for(r=0;r<NUniformRefs;r++) {
			// Printing computational mesh to information
			if(comp->NElements() < 200)
				comp->Print(std::cout);
			else {
				std::cout << "Computacional mesh : NElements = " << comp->NElements() << "\t NConnects = " << comp->NConnects() << std::endl;
			}
			
			// Introduzing exact solution depending on the case
			TPZAnalysis an (comp);
			an.SetExact(Exact);		
			{   // To print solution
				std::stringstream sout;
				int angle = (int) (alfa*180./M_PI + 0.5);
				if(anothertests) sout << "Leg_";
				sout << "hptestAngo" << angle << "." << r << ".vtk";
				an.DefineGraphMesh(dim,scalnames,vecnames,sout.str());
			}
			std::string MeshFileName;
			{   // To print computational mesh
				std::stringstream sout;
				int angle = (int) (alfa*180./M_PI + 0.5);
				if(anothertests) sout << "Leg_";
				sout << "meshAngle" << angle << "." << r << ".vtk";
				MeshFileName = sout.str();
			}
			comp->SetName("Malha computacional adaptada");
			
			// Solve using symmetric matrix then using Cholesky (direct method)
			TPZSkylineStructMatrix strskyl(comp);
			an.SetStructuralMatrix(strskyl);
			
			TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
			direct->SetDirect(ECholesky);
			an.SetSolver(*direct);
			delete direct;
			direct = 0;
			
			int neq = comp->NEquations();
			//	an.NEquations();
			//		an.Solution().Print();
			an.Run();
			
			// Computing approximation of gradient
			/** 
			 * @brief Method to reconstruct a gradient after run Solve of the analysis
			 * @param cmesh Computational mesh with solution */
			TPZFMatrix<REAL> gradients;
			GradientReconstructionByLeastSquares(gradients,comp,0,0,true);
			gradients.Print();
			
			// Post processing
			an.PostProcess(1,dim);
			{
				std::ofstream out(MeshFileName.c_str());
				comp->LoadReferences();
				TPZVTKGeoMesh::PrintCMeshVTK(comp->Reference(), out, false);
			}
			
			return 0;
		}
	}
}

int main_AdaptHP_3D(int argc, char *argv[]) {
	
#ifdef LOG4CXX
	if (argc > 1) {
		std::string logpath ( argv[1] );
		cout << "initializing LOG usign the following configuration file " << logpath << endl;
		InitializePZLOG ( logpath );	
	} else {
		cout << "initializing LOG\n";
		InitializePZLOG();
	}
#endif
	const int L = 4;
	gLMax = L-1;
	char saida[260];
	
	//-----------  INITIALIZING CONSTRUCTION OF THE MESHES
	
	int r, dim;
	
	// Initializing a ref patterns
	gRefDBase.InitializeAllUniformRefPatterns();
	//gRefDBase.InitializeRefPatterns();
	//gRefDBase.InitializeUniformRefPattern(EOned);
	//gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	//gRefDBase.InitializeUniformRefPattern(ETriangle);
	// Inserting a special file with refinement pattern 
	std::string filename = REFPATTERNDIR;
	filename += "/3D_Hexa_Rib_Side_16_16_18_18.rpt";
	//filename += "/3D_Hexa_Face_20.rpt";
	
	TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
	if(!gRefDBase.FindRefPattern(refpat))
	{
		gRefDBase.InsertRefPattern(refpat);
	}
	refpat->InsertPermuted();
	
	// output files
    std::ofstream convergence("conv3d.txt");
    std::ofstream out("output.txt");
	
    // First rectangular mesh:
	// The rectangular mesh has four corners: (0,0,0), (1,0,0), (1,1,0) and (0,1,0)
	// and was divides in two segments on X and two on Y, then hx = 0.5 and hy = 0.5
	// Has 4 elements, 9 connects and 8 bc elements
	cout << "Generating geometric mesh bi-dimensional ...\n";
    TPZGeoMesh* gmesh = new TPZGeoMesh;
	TPZManVector<REAL> x0(3,0.), x1(3,1.);  // Corners of the rectangular mesh. Coordinates of the first extreme are zeros.
	TPZManVector<int> nx(2,2);   // subdivisions in X and in Y. 
	TPZGenGrid gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1  
	gen.SetElementType(0);       // type = 0 means rectangular elements
	gen.Read(gmesh);             // generating grid in gmesh
	
	// Extending geometric mesh (two-dimensional) to three-dimensional geometric mesh
	// The elements are hexaedras(cubes) over the quadrilateral two-dimensional elements
	cout << "Generating geometric mesh three-dimensional (extruding) ...\n";
	TPZExtendGridDimension gmeshextend(gmesh,0.5);
	TPZGeoMesh *gmesh3D = gmeshextend.ExtendedMesh(2,2,2);
	
	// Applying hp adaptive techniques 2012/10/01
	if(anothertests) {
		// Setting Chebyshev polynomials as orthogonal sequence generating shape functions
		TPZShapeLinear::fOrthogonal = &TPZShapeLinear::Legendre;
		sprintf(saida,"meshextrudedLeg.vtk");
		
	}
	else {
		sprintf(saida,"meshextrudedTChe.vtk");
	}	
	// Uniform Refinement - Some times for three dimensional elements
	//	UniformRefinement(2,gmesh3D,3);
	// Refinement of the some element	
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	TPZVec<TPZGeoEl *> subsub;
	int nele;
	//	for(int ii=0;ii<3;ii++) {
	//		int ngelem = gmesh->NElements()-1;
	nele = 0;
	//		for(;nele<ngelem;nele++) {
	gel = gmesh3D->ElementVec()[nele];
	
	//			if(gel->Dimension() != 3) continue;
	//	gel->SetRefPattern(refpat);
	
	gel->Divide(sub);
	//			int jj = 0;
	//			for(jj=0;jj<4;jj++) {
	//				gel = sub[jj];
	//				gel->Divide(subsub);
	//			}
	//		}
	//			TPZVec<REAL> coord(3,0.);
	//			TPZVec<REAL> point(3,-1.);
	//			gel->X(point,coord);
	//			if(!IsZero(coord[0]) || !IsZero(coord[1]) || !IsZero(coord[2])) continue;
	//			if(gel->Dimension() != 3) continue;
	//			gel->SetRefPattern(refpat);
	gel->Divide(sub);
	//			for(int jj=0;jj<4;jj++) {
	//				gel = sub[jj];
	//				if(gel->Dimension() != 3) continue;
	//				gel->SetRefPattern(refpat);
	//				gel->Divide(subsub);
	//				for(int kk=0;kk<gel->NSubElements();kk++)
	//					subsub[kk]->SetRefPattern(refpat);
	//			}
	//			gel = subsub[0];
	//		gel->Divide(sub);
	//		gel = sub[0];
	//		gel->Divide(subsub);
	//		gel = subsub[0];
	//		gel->Divide(sub);
	//			break;
	//		}
	//	}
	// Constructing connectivities
	gmesh3D->ResetConnectivities();
	gmesh3D->BuildConnectivity();
	gmesh3D->Print();
	// Printing COMPLETE initial geometric mesh 
	PrintGeoMeshVTKWithDimensionAsData(gmesh3D,saida);
	
    // Creating computational mesh
    TPZCompMesh *comp = new TPZCompMesh(gmesh3D);
	/** Set polynomial order */
	int p = 2;
    TPZCompEl::SetgOrder(p);
	
	TPZVec<REAL> forces(3,0.);
	// Creating and inserting materials into computational mesh
    //TPZMaterial * mat = new TPZElasticityMaterial(1,1.e5,0.2,0,0);   // two-dimensional
	TPZMaterial *mat = new TPZElasticity3D(1,1.e5,0.2,forces);          // three-dimensional
    comp->InsertMaterialObject(mat);
	dim = mat->Dimension();
	nstate = mat->NStateVariables();
	
    // Boundary conditions
    // Dirichlet
    TPZFMatrix<REAL> val1(3,3,0.),val2(3,1,5.);
	val1(0,0) = val1(1,1) = 1.;
    TPZMaterial *bnd = mat->CreateBC(mat,-1,0,val1,val2);
    comp->InsertMaterialObject(bnd);
	// Neumann
    val2(0,0)=30.; val2(1,0) = 10.;
    bnd = mat->CreateBC(mat,-2,1,val1,val2);
    comp->InsertMaterialObject(bnd);
    
    // Constructing and adjusting computational mesh
    comp->AutoBuild();
    comp->AdjustBoundaryElements();   // Adjust boundary elements and higher level of refinement, clean elements but not connects into them
    comp->CleanUpUnconnectedNodes();  // Clean connects not connected at least one element enabled.
	
	//--- END construction of the meshes
	
	/** Variable names for post processing */
    TPZStack<std::string> scalnames, vecnames;
	if(mat->NSolutionVariables(mat->VariableIndex("POrder")) == 1)
		scalnames.Push("POrder");
	else
		vecnames.Push("POrder");
	if(mat->NSolutionVariables(mat->VariableIndex("Error")) == 1)
		scalnames.Push("Error");
	else
		vecnames.Push("Error");
	if(mat->NSolutionVariables(mat->VariableIndex("state")) == 1)
		scalnames.Push("state");
	else
		vecnames.Push("state");
    
    if(nstate == 1) {
        scalnames.Push("TrueError");
        scalnames.Push("EffectivityIndex");
    }else if(nstate == 2) {
        scalnames.Push("sig_x");
        scalnames.Push("sig_y");
        scalnames.Push("tau_xy");
    }
    if(nstate == 3) {
        scalnames.Push("StressX");
        scalnames.Push("StressY");
        scalnames.Push("StressZ");
		vecnames.Push("PrincipalStress");
		vecnames.Push("PrincipalStrain");
    }
	
	// END Determining the name of the variables
	
	// TO MAKE MERGE ANOTHER DOMAIN
	if(anothertests) {
		// Second rectangular domain - subdivisions and corners of the second rectangular mesh
		TPZAutoPointer<TPZGeoMesh> gmesh2 = new TPZGeoMesh;
		x0[1] = 0.2;                 // left and right extremes of the new geo mesh. Coordinates: (0.,0.2,0.0) (3.,1.,0.) 
		x1[0] = 3.; x1[1] = 1.;
		nx[0] = 15; nx[1] = 8;       // subdivision in X and Y. hx = 0.2 and hy = 0.1
		TPZGenGrid gen2(nx,x0,x1);   // second mesh generator
		gen2.SetElementType(0);      // type = 0 means rectangular elements, type = 1 means triangular elements
		
		// Generating gmesh2 with last data and after this the gmesh is merged into the gmesh2. But gmesh is unmodified
		// if exist boundary elements into the mesh merged it will be deleted
		gen2.ReadAndMergeGeoMesh(gmesh2,gmesh);
		
		// setting bc condition -1 [no flux - is wall] from (0.0, 0.0) until (3.0, 0.2)
		x0[1] = 0.0; x1[1] = 0.2;
		gen2.SetBC(gmesh2,x0,x1,-1);
		// setting bc condition -1 from (3.0, 1.0) until (0.0, 1.0)
		x0[0] = 3.; x0[1] = 1.;
		x1[0] = 0.;	x1[1] = 1.;
		gen2.SetBC(gmesh2,x0,x1,-1);
		// setting bc condition -2 [free flux] from (3.0, 1.0) until (3.0, 0.2)
		x1[0] = 3.;	x1[1] = .2;
		gen2.SetBC(gmesh2,x1,x0,-2);
		// setting bc condition -2 [free flux] from (0.0, 1.0) until (0.0, 0.0)
		x0[0] = 0.;
		x1[0] = x1[1] = 0.;
		gen2.SetBC(gmesh2, x0, x1, -2);
		
#ifdef DEBUG
		sprintf(saida,"original.vtk");
		PrintGeoMeshVTKWithDimensionAsData(gmesh,saida);
		sprintf(saida,"meshes.vtk");
		PrintGeoMeshVTKWithDimensionAsData(gmesh2.operator->(),saida);
#endif
		// Extending geometric mesh (two-dimensional) to three-dimensional geometric mesh
		// The elements are hexaedras(cubes) over the quadrilateral two-dimensional elements
		TPZExtendGridDimension gmeshextend(gmesh2,0.3);
		TPZGeoMesh *gmesh3D = gmeshextend.ExtendedMesh(2,-5,-6);
#ifdef DEBUG
		sprintf(saida,"meshextrudedend.vtk");
		PrintGeoMeshVTKWithDimensionAsData(gmesh3D,saida);
#endif
		dim = 3;
	}
	
	// INITIAL POINT FOR SOLVING AND APPLYING REFINEMENT
	for(r=0;r<NUniformRefs;r++) {
		// Printing computational mesh to information
		if(comp->NElements() < 200)
			comp->Print(std::cout);
		else {
			std::cout << "Computacional mesh : NElements = " << comp->NElements() << "\t NConnects = " << comp->NConnects() << std::endl;
		}
		
		// Introduzing exact solution depending on the case
		TPZAnalysis an (comp);
		an.SetExact(Exact);		
		{   // To print solution
			std::stringstream sout;
			int angle = (int) (alfa*180./M_PI + 0.5);
			if(anothertests) sout << "Leg_";
			sout << "hptestAngo" << angle << "." << r << ".vtk";
			an.DefineGraphMesh(dim,scalnames,vecnames,sout.str());
		}
		std::string MeshFileName;
		{   // To print computational mesh
			std::stringstream sout;
			int angle = (int) (alfa*180./M_PI + 0.5);
			if(anothertests) sout << "Leg_";
			sout << "meshAngle" << angle << "." << r << ".vtk";
			MeshFileName = sout.str();
		}
		comp->SetName("Malha computacional adaptada");
		
		// Solve using symmetric matrix then using Cholesky (direct method)
		TPZSkylineStructMatrix strskyl(comp);
		an.SetStructuralMatrix(strskyl);
		
		TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
		direct->SetDirect(ECholesky);
		an.SetSolver(*direct);
		delete direct;
		direct = 0;
		
		an.Run();
		
		// Post processing
		an.PostProcess(1,dim);
		{
			std::ofstream out(MeshFileName.c_str());
			comp->LoadReferences();
			TPZVTKGeoMesh::PrintCMeshVTK(comp->Reference(), out, false);
		}
		
		/*
		 REAL valerror =0.;
		 REAL valtruerror=0.;
		 TPZVec<REAL> ervec,truervec,effect;
		 
		 TPZAdaptMesh adapt;
		 adapt.SetCompMesh (comp);
		 
		 std::cout << "\n\n\n\nEntering Auto Adaptive Methods... step " << r << "\n\n\n\n";
		 
		 time_t sttime;
		 time (& sttime);
		 TPZCompMesh *adptmesh;
		 
		 adptmesh = adapt.GetAdaptedMesh(valerror,valtruerror,ervec,Exact,truervec,effect,0);
		 
		 time_t endtime;
		 time (& endtime);
		 
		 int time_elapsed = endtime - sttime;
		 std::cout << "\n\n\n\nExiting Auto Adaptive Methods....step " << r
		 << "time elapsed " << time_elapsed << "\n\n\n\n";
		 
		 int prt;
		 std::cout << "neq = " << comp->NEquations() << " error estimate = " << valerror
		 << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;
		 
		 #ifdef LOG4CXX
		 if (loggerconv->isDebugEnabled())
		 {
		 std::stringstream sout;
		 sout << "neq = " << comp->NEquations() << " error estimate = " << valerror
		 << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;
		 LOGPZ_DEBUG(loggerconv, sout.str())
		 }
		 #endif
		 
		 convergence  << comp->NEquations() << "\t"
		 << valerror << "\t"
		 << valtruerror << "\t"
		 << ( valtruerror / valerror ) <<  "\t"
		 << sttime <<std::endl;
		 for (prt=0;prt<ervec.NElements();prt++){
		 std::cout <<"error " << ervec[prt] << "  truerror = " << truervec[prt] << "  Effect " << effect[prt] << std::endl;
		 // convergence << '\t' << ervec[prt] << '\t' << truervec[prt] << "  Effect " << effect[prt] <<  std::endl;
		 //  adptmesh->Print(cout);
		 }
		 
		 std::cout.flush();
		 comp->Reference()->ResetReference();
		 comp->LoadReferences();
		 adapt.DeleteElements(comp);
		 delete comp;
		 comp = adptmesh;
		 
		 comp->CleanUpUnconnectedNodes();
		 */
	}
	/* Uniform refinement. Two times
	 UniformRefinement(2,gmesh3D,3);
	 
	 #ifdef DEBUG
	 sprintf(saida,"meshrefined.vtk");
	 PrintGeoMeshVTKWithDimensionAsData(gmesh3D,saida);
	 #endif
	 
	 // Creating a computational mesh (interpolation space)
	 TPZCompMesh *cmesh = CreateMeshMultires(gmesh3D);
	 #ifdef DEBUG
	 sprintf(saida,"aftercmesh.vtk");
	 PrintGeoMeshVTKWithDimensionAsData(gmesh3D,saida);
	 #endif
	 
	 // To work with the temporal variable.
	 REAL timeStep;
	 timeStep = ComputeTimeStep(1.,L,cmesh->Reference());
	 
	 #ifdef DEBUG
	 {
	 ofstream malhas("malhas.vtk");
	 cmesh->Print(malhas);
	 }
	 #endif
	 
	 TPZExplFinVolAnal an(cmesh, cout);
	 
	 InitializeSolver(an);
	 const double PhysicalTime = 0.1;
	 int niter = PhysicalTime/timeStep+1;
	 cout << "L = " << L << endl;
	 cout << "\nnequations = " << cmesh->NEquations();
	 cout << "\nNiter = " << niter << "\n";
	 
	 TPZFMatrix<REAL> InitialSol;
	 InitialSolutionLinearConvection(InitialSol,cmesh);
	 
	 an.SetInitialSolution(InitialSol);
	 
	 an.Set(timeStep,niter,1e-10);
	 an.SetSaveFrequency(niter/6,0);
	 
	 double Epsl = 1.e12;
	 an.MultiResolution( Epsl );
	 #ifdef DEBUG
	 sprintf(saida,"meshInitialSol.vtk");
	 TPZVTKGeoMesh::PrintGMeshVTK(gmesh3D,saida,0);
	 #endif
	 
	 return EXIT_SUCCESS;
	 */
	return 0;
}

/** Exact solutions to calculate the rate of convergence */

static REAL onethird = 0.33333333333333333;
static REAL PI = 3.141592654;

void Exact(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
    TPZManVector<REAL,3> x2(x);
	//    TransformInvX(x2,RotInv);
  	REAL r = sqrt(x2[0]*x2[0]+x2[1]*x2[1]);
  	REAL theta = atan2(x2[1],x2[0]);
  	REAL rexp = pow(r,onethird);
  	sol[0] = rexp*sin(onethird*(theta+PI/2));
    TPZFNMatrix<3,REAL> grad(4,1,0.),grad2(4,1,0.);
  	grad(0,0) = onethird*sin(onethird*(PI/2.-2.*theta))/(rexp*rexp);
  	grad(1,0) = onethird*cos(onethird*(PI/2.-2.*theta))/(rexp*rexp);
	//    Rot.Multiply(grad, grad2);
    dsol(0,0) = grad2(0,0);
    dsol(1,0) = grad2(1,0);
}

void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename) {
	int i, size = gmesh->NElements();
	TPZChunkVector<int> DataElement;
	DataElement.Resize(size);
	// Making dimension of the elements as data element
	for(i=0;i<size;i++) {
		if(gmesh->ElementVec()[i])
			DataElement[i] = (gmesh->ElementVec()[i])->Dimension();
		else
			DataElement[i] = -999;
	}
	// Printing geometric mesh to visualization in Paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filename, DataElement);
}

void InitializeSolver(TPZAnalysis &an) {
	TPZStepSolver<REAL> step;
	TPZBandStructMatrix matrix(an.Mesh());
	an.SetStructuralMatrix(matrix);
	step.SetDirect(ELU);
	an.SetSolver(step);
}

void InitialSolutionLinearConvection(TPZFMatrix<REAL> &InitialSol, TPZCompMesh * cmesh){
	InitialSol.Redim(cmesh->NEquations(),1);
	InitialSol.Zero();
	for(int iel = 0; iel < cmesh->NElements(); iel++){
		TPZCompEl * cel = cmesh->ElementVec()[iel];
		if(!cel) continue;
		TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cel);
		if(!disc) continue;
		if(disc->NConnects() == 0) continue;
		int bl = disc->Connect(0).SequenceNumber();
		int blpos = cmesh->Block().Position(bl);
		int blocksize = cmesh->Block().Size(bl);
		
		TPZGeoEl * gel = cel->Reference();
		TPZVec<REAL> xi(3), xVec(3);
		gel->CenterPoint(gel->NSides()-1,xi);
		gel->X(xi,xVec);
		double x = xVec[0];
		double y = xVec[1];
		double u = 0.125;
		
		double xCircle = 0.25;
		double yCircle = 0.5;
		double R = 0.1;
		if( (x-xCircle)*(x-xCircle)+(y-yCircle)*(y-yCircle) <= R*R ) u = 1.;
		
		InitialSol(blpos+blocksize-20+0,0) = u;
		InitialSol(blpos+blocksize-20+1,0) = 0.;
		InitialSol(blpos+blocksize-20+2,0) = 0.;
		InitialSol(blpos+blocksize-20+3,0) = 0.;
		InitialSol(blpos+blocksize-20+4,0) = 0.;
		
	}//for iel
	
	TPZVec<REAL> celerity(3,0.);
	celerity[0] = 1.;
#ifdef LinearConvection
	TPZEulerEquation::SetLinearConvection(cmesh, celerity);
#endif
	
}//method

void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial, const int matidtodivided) {
	TPZManVector<TPZGeoEl*> filhos;
    for(int D=0; D<nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {    
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
			if(!gel || gel->HasSubElement())
				continue;
			if(dim > 0 && gel->Dimension() != dim) continue;
			if(!allmaterial){
				if(gel->MaterialId() == matidtodivided){
					gel->Divide(filhos);
				}
			}
			else{
				gel->Divide(filhos);
			}
        }
    }
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}
