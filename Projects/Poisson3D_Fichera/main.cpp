/**
 * @file
 */

#include "pzshapelinear.h"

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZVTKGeoMesh.h"

//#include "MultiResMesh.h"
#include "pzbstrmatrix.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzpoisson3d.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"

#include "TPZRefPatternTools.h"

#include "TPZParSkylineStructMatrix.h"

#include <stdio.h>

//#include "pzexplfinvolanal.h"
#include "adapt.h"

#include "pzlog.h"

#include "pzgeoelbc.h"

#ifdef LOG4CXX
static LoggerPtr loggerconv(Logger::getLogger("pz.adaptivity.conv"));
static LoggerPtr loggerpoint(Logger::getLogger("pz.adaptivity.points"));
#endif

using namespace std;
using namespace pzshape;

/**
 * @addtogroup Tutorials
 * @{
 */

// Global variable
int gLMax;
int NUniformRefs = 2;
REAL alfa = M_PI/6.;
bool anothertests = false;
int nstate = 2;

// output files  -> Because it has many energy faults
std::ofstream out("output.txt");

/** Printing level */
int gPrintLevel = 0;
bool gDebug = false;

//TPZFNMatrix<16,REAL> Rot(4,4,0.),RotInv(4,4,0.);
void Exact(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);

void ExactSolin(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);

void BCSolin(const TPZVec<REAL> &x, TPZVec<REAL> &sol);

void Ff(const TPZVec<REAL> &x, TPZVec<REAL> &f);

//void InitializeSolver(TPZAnalysis &an);
void InitialSolutionLinearConvection(TPZFMatrix<REAL> &InitialSol, TPZCompMesh *cmesh);
void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);

void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial=true, const int matidtodivided=1);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL> &points,REAL r,REAL &distance,bool &isdefined);

TPZGeoMesh *ConstructingFicheraCorner(REAL L, REAL H,bool print = false);
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction);

void ComputeSolutionError(REAL &error,REAL &errorL2,TPZCompMesh *cmesh,int dim);
void formatTimeInSec(char *strtime,int timeinsec);

int main(int argc, char *argv[]) {
	
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing a ref patterns
	gRefDBase.InitializeAllUniformRefPatterns();
    //	gRefDBase.InitializeRefPatterns();
	// To compute processing times
	time_t sttime;
	time_t endtime;
	int time_elapsed;
	char tempo[256];
    
	// To compute the errors
	ofstream fileerrors("ErrorsHPProcess.txt");
    char saida[256];
    
	TPZVec<REAL> ervec(100,0.0);
	// Printing computed errors
	fileerrors << "Approximation Error: " << std::endl;
	
	//-----------  INITIALIZING CONSTRUCTION OF THE MESHES
	REAL InitialL = 1.0, InitialH = 1.;
	int i, nref, NRefs = 7;
	int nthread, NThreads = 4;
	int dim = 3;
    
    for(nref=2;nref<NRefs;nref++) {
		//        for(nthread=1;nthread<NThreads;nthread++) {
		if(nref > 4) nthread = NThreads;
		else nthread = 1;
		
		// Initializing the generation mesh process
        time(&sttime);
        // Constructing geometric mesh as Fichera corner using hexahedra
        TPZGeoMesh *gmesh3D = ConstructingFicheraCorner(InitialL, InitialH);
        // h_refinement
        TPZManVector<REAL> point(3,0.);
        REAL r = 0.0, radius = 0.9;
        bool isdefined = false;
        for(i=0;i<nref+1;i++) {
            // Para refinar elementos con centro tan cerca de la circuferencia cuanto radius 
            RefineGeoElements(3,gmesh3D,point,r,radius,isdefined);
            radius *= 0.6;
        }
        sprintf(saida,"meshextrudedmerged_%d.vtk",nref);
        PrintGeoMeshVTKWithDimensionAsData(gmesh3D,saida);
        
        // Creating computational mesh
        /** Set polynomial order */
        int p = 2, pinit;
        pinit = p;
        //            TPZCompEl::SetgOrder(p);
        TPZCompMesh *cmesh = CreateMesh(gmesh3D,dim,1);
		cmesh->SetName("Computational mesh for Fichera problem");
		dim = cmesh->Dimension();
        
		/* Primeiro sera calculado o mayor nivel de refinamento. Remenber, the first level is zero level.
		// A cada nivel disminue em uma unidade o p, mas não será menor de 1.
		int level, highlevel = 0;
		int nelem = 0;
		while(nelem < cmesh->NElements()) {
			TPZCompEl *cel = cmesh->ElementVec()[nelem++];
			if(cel) {
				level = cel->Reference()->Level();
			}
			if(level > highlevel)
				highlevel = level;
		}
		// Identifying maxime interpolation order
		if(highlevel>p-1) pinit = p;
		else pinit = highlevel+1;
		// Put order 1 for more refined element and (highlevel - level)+1 for others, but order not is greater than initial p
		nelem = 0;
		while(highlevel && nelem < cmesh->NElements()) {
			TPZCompEl *cel = cmesh->ElementVec()[nelem++];
            if(!cel) continue;
			level = cel->Reference()->Level();
			p = (highlevel-level)+1;
			if(p > pinit) p = pinit;
			((TPZInterpolatedElement*)cel)->PRefine(p);
		}*/
		cmesh->SetAllCreateFunctionsContinuous();
		cmesh->AutoBuild();
		cmesh->AdjustBoundaryElements();
        cmesh->ExpandSolution();
        
		// closed generation mesh process
		time (& endtime);
		time_elapsed = endtime - sttime;
		formatTimeInSec(tempo, time_elapsed);
		out << "  Time elapsed " << time_elapsed << " <-> " << tempo << "\n\n";
		
		//--- END construction of the meshes
        
		// Solving linear equations
		// Initial steps
		out << "Solving HP-Adaptive Methods....step: " << nref << "  Threads " << nthread << "\n";
        
		/** Variable names for post processing */
        TPZStack<std::string> scalarnames, vecnames;
		scalarnames.Push("POrder");
		scalarnames.Push("Solution");
		scalarnames.Push("KDuDx");
		scalarnames.Push("KDuDy");
		scalarnames.Push("KDuDz");
		scalarnames.Push("NormKDu");
		scalarnames.Push("Pressure");
		
		vecnames.Push("Derivate");
		vecnames.Push("Flux");
		vecnames.Push("MinusKGradU");
		// END Determining the name of the variables
        
        // Introduzing exact solution
        TPZAnalysis an (cmesh);
        
        // Solve using symmetric matrix then using Cholesky (direct method)
        TPZParSkylineStructMatrix strskyl(cmesh,nthread);
        an.SetStructuralMatrix(strskyl);
        
        TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
        direct->SetDirect(ECholesky);
        an.SetSolver(*direct);
        delete direct;
        direct = 0;
        
		// Initializing the solving process
		time (& sttime);
        an.Run();
        time(&endtime);
        time_elapsed = endtime - sttime;
        formatTimeInSec(tempo,time_elapsed);
        
        char pp[3];
        sprintf(pp,"%d",pinit);
        out << "\t....step: " << nref << "  Threads " << nthread << "  Time elapsed " << time_elapsed << " <-> " << tempo << "\n\n\n";
        
        // Post processing
        std::string filename = "Poisson3DSol";
        filename += "_p";
        filename += pp;
        filename += "_hL";
        sprintf(pp,"%d",nref);
        filename += pp;
		if(nthread > 1) {
            sprintf(pp,"_PAR%d",nthread);
            filename += pp;
		}
        filename += ".vtk";
        an.DefineGraphMesh(dim,scalarnames,vecnames,filename);
        
        an.PostProcess(0,dim);
        
        // Computing error
        an.SetExact(ExactSolin);
        fileerrors << "Refinement: " << nref << "  Dimension: " << dim << "  NEquations: " << cmesh->NEquations();
        an.PostProcess(ervec,out);
        for(int rr=0;rr<ervec.NElements();rr++)
            fileerrors << "  Error_" << rr+1 << ": " << ervec[rr]; 
        fileerrors << "  TimeElapsed: " << time_elapsed << " <-> " << tempo << std::endl;
        
        delete cmesh;
        delete gmesh3D;
        break;
        //        }
	}
    out.close();
    fileerrors.close();
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
//////////   FICHERA CORNER - Problem as Anders Solin Presentation   ///////////////////
////////////////////////////////////////////////////////////////////////////////////////

void ComputeSolutionError(REAL &error,REAL &errorL2,TPZCompMesh *cmesh,int dim) {
	int i, it, nelem = cmesh->NElements();
	
	TPZBlock<REAL> flux;
	
	TPZVec<REAL> point(3,0.), x(3,0.);
	REAL weight = 0.;
	
	TPZInterpolatedElement *cel;
	int var = 0;                       // For solution of "state" or "Displacement"
	TPZVec<REAL> SolCel(5,0.0);
	TPZVec<REAL> SolCelExact(5,0.0);
	TPZFMatrix<REAL> DSolCel(5,5,0.0);
    
	error = 0.; errorL2 = 0.;
    REAL errorLoc, errorLocL1;
	for(i=0;i<nelem;i++) {
        errorLoc = errorLocL1 = 0.;
		cel = (TPZInterpolatedElement *)cmesh->ElementVec()[i];
		if(!cel || cel->Dimension() != dim) continue;
		int npoints = cel->GetIntegrationRule().NPoints();
		for(it=0;it<npoints;it++){
			cel->GetIntegrationRule().Point(it,point,weight);
			cel->Reference()->X(point,x);
			cel->Solution(point,var,SolCel);
			ExactSolin(x,SolCelExact,DSolCel);
			errorLocL1 += weight * fabs(SolCel[0] - SolCelExact[0]);
			errorLoc += weight * (SolCel[0] - SolCelExact[0])*(SolCel[0] - SolCelExact[0]);
		}
		errorLocL1 *= fabs(cel->Reference()->Volume());
		errorLoc *= fabs(cel->Reference()->Volume());
        error += fabs(errorLocL1);
        errorL2 += errorLoc;
	}
	errorLoc = errorL2;
    errorL2 = sqrt(errorLoc);
}

void ExactSolin(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
	REAL quad_r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	sol[0] = sqrt( sqrt (quad_r) );
	if(!IsZero(sol[0])) {
		REAL den = sol[0]*sol[0]*sol[0];
		dsol(0,0) = .5*x[0]/den;
		dsol(1,0) = .5*x[1]/den;
		dsol(1,0) = .5*x[2]/den;
	}
	else {
		dsol(0,0) = dsol(1,0) = dsol(2,0) = 0.;
	}
}

void BCSolin(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol) {
	REAL quad_r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	bcsol[0] = sqrt( sqrt (quad_r) );	
}

void Ff(const TPZVec<REAL> &x, TPZVec<REAL> &f) {
	REAL quad_r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	REAL raiz = sqrt( sqrt(quad_r));
	f[0] = -3./(4.*(raiz*raiz*raiz));
}

TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction) {
	
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();
	
	// Creating Poisson material
	TPZMaterial *mat = new TPZMatPoisson3d(1,dim);
	TPZVec<REAL> convd(3,0.);
	((TPZMatPoisson3d *)mat)->SetParameters(1.,0.,convd);
	if(hasforcingfunction) {
		mat->SetForcingFunction(new TPZDummyFunction<STATE>(Ff));
	}
	cmesh->InsertMaterialObject(mat);
	// Make compatible dimension of the model and the computational mesh
	cmesh->SetDimModel(mat->Dimension());
	cmesh->SetAllCreateFunctionsContinuous();
    
	// Boundary conditions
	// Dirichlet
	TPZAutoPointer<TPZFunction<STATE> > FunctionBC = new TPZDummyFunction<STATE>(BCSolin);
	TPZFMatrix<REAL> val1(dim,dim,0.),val2(dim,1,0.);
	TPZMaterial *bnd = mat->CreateBC(mat,-1,0,val1,val2);
	bnd->SetForcingFunction(FunctionBC);
	cmesh->InsertMaterialObject(bnd);
	
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->ExpandSolution();
	return cmesh;
}

TPZGeoMesh *ConstructingFicheraCorner(REAL InitialL, REAL InitialH, bool print) {
	TPZManVector<REAL> x0(3,-InitialL), x1(3,-InitialL);  // Corners of the rectangular mesh.
	TPZManVector<int> nx(2,2);   // subdivisions in X and in Y. 
	
	// Mesh as a hexahedra with large edge (2 units) and was cutting a small hexahedra with edge with 1 unit:
	// 1) The bottom rectangular mesh has four corners: (1,-1,-1), (1,1,-1), (-1,1,-1) and (-1,-1,-1)
	// and was divides in four segments on X and four on Y, then hx = 0.5 and hy = 0.5
	if(print) out << "Generating (1) bottom geometric mesh bi-dimensional ...\n";
	TPZGeoMesh* gmesh = new TPZGeoMesh;
	x0[0] = InitialL; x0[1] = x0[2] = -InitialL;
	x1[0] = x1[2] = - InitialL; x1[1] = InitialL;
	TPZGenGrid gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1  
	gen.SetElementType(0);       // type = 0 means rectangular elements
	gen.Read(gmesh);             // generating grid in gmesh
	// Extending geometric mesh (two-dimensional) to three-dimensional geometric mesh
	// The elements are hexaedras(cubes) over the quadrilateral two-dimensional elements
	if(print) out << "... (Extruding) first geometric mesh three-dimensional...\n";
	TPZExtendGridDimension gmeshextend(gmesh,InitialH);
	TPZGeoMesh *gmesh3D = gmeshextend.ExtendedMesh(1,-1,2);
	gmesh3D->SetName("First Mesh Extruded");
	// Dirichlet boundary condition for inferior cubes
	TPZGeoElBC gelbcD(gmesh3D->ElementVec()[0],21,-1);   // Dirichlet condition for hexahedra//
	TPZGeoElBC gelbcD1(gmesh3D->ElementVec()[0],24,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD2(gmesh3D->ElementVec()[1],21,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD3(gmesh3D->ElementVec()[1],22,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD4(gmesh3D->ElementVec()[2],23,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD5(gmesh3D->ElementVec()[2],24,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD6(gmesh3D->ElementVec()[2],25,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD7(gmesh3D->ElementVec()[3],22,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD8(gmesh3D->ElementVec()[3],23,-1);   // Dirichlet condition for hexahedra
	
	// 2) The left upper rectangular mesh has four corners: (1,-1,0), (1,0,0), (-1,0,0) and (-1,-1,0)
	// and was divides in one segment on X and two on Y, then hx = hy = 1.
	TPZGeoMesh *gmesh2 = new TPZGeoMesh;
	if(print) out << "Generating (2) left upper geometric mesh bi-dimensional ...\n";
	x0[0] = InitialL; x0[1] = -InitialL; x0[2] = 0.;
	x1[0] = -InitialL; x1[1] = x1[2] = 0.;
	nx[0] = 2; nx[1] = 1;
	gen.SetData(nx,x0,x1,0);
	gen.Read(gmesh2);             // generating grid in gmesh
	// Extending geometric mesh (two-dimensional) to three-dimensional geometric mesh
	// The elements are hexaedras(cubes) over the quadrilateral two-dimensional elements
	if(print) out << "... (Extruding) first geometric mesh three-dimensional...\n";
	TPZExtendGridDimension gmeshextend2(gmesh2,InitialH);
	TPZGeoMesh *gmesh3D2 = gmeshextend2.ExtendedMesh(1,2,-1);
	gmesh3D2->SetName("Second Mesh Extruded");
	// Dirichlet boundary condition for inferior cubes
	TPZGeoElBC gelbcD10(gmesh3D2->ElementVec()[0],21,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD11(gmesh3D2->ElementVec()[0],23,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD12(gmesh3D2->ElementVec()[0],24,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD13(gmesh3D2->ElementVec()[1],21,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD14(gmesh3D2->ElementVec()[1],22,-1);   // Dirichlet condition for hexahedra
	
	// 3) The right upper rectangular mesh has four corners: (0,0,0), (0,1,0), (-1,1,0) and (-1,0,0)
	// and was divides in one segment on X and two on Y, then hx = hy = 1.
	TPZGeoMesh *gmesh3 = new TPZGeoMesh;
	if(print) out << "Generating (3) left upper geometric mesh bi-dimensional ...\n";
	x0[0] = x0[1] = x0[2] = 0.;
	x1[0] = -InitialL; x1[1] = InitialL; x1[2] = 0.;
	nx[0] = nx[1] = 1;
	gen.SetData(nx,x0,x1,0);       // type = 0 means rectangular elements
	gen.Read(gmesh3);             // generating grid in gmesh
	// Extending geometric mesh (two-dimensional) to three-dimensional geometric mesh
	// The elements are hexaedras(cubes) over the quadrilateral two-dimensional elements
	if(print) out << "... (Extruding) first geometric mesh three-dimensional...\n";
	TPZExtendGridDimension gmeshextend3(gmesh3,InitialH);
	TPZGeoMesh *gmesh3D3 = gmeshextend3.ExtendedMesh(1,2,-1);
	gmesh3D3->SetName("Third Mesh Extruded");
	// Dirichlet boundary condition for inferior cubes
	TPZGeoElBC gelbcD20(gmesh3D3->ElementVec()[0],22,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD21(gmesh3D3->ElementVec()[0],23,-1);   // Dirichlet condition for hexahedra
	TPZGeoElBC gelbcD22(gmesh3D3->ElementVec()[0],24,-1);   // Dirichlet condition for hexahedra
	
	// Merge two meshes
	gen.MergeGeoMesh(gmesh3D,gmesh3D2,1);
	gen.MergeGeoMesh(gmesh3D,gmesh3D3,1);
	
	// Cleaning unnecessary meshes
	delete gmesh3D2;
	delete gmesh3D3;
	
	return gmesh3D;
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

#include "TPZGenSpecialGrid.h"

int main2() {
	// To visualization of the geometric mesh
	std::ofstream fgeom("GeoMeshByTolerance.vtk");
	std::ofstream fgeom2("GeoMeshByNRefinements.vtk");
	
	/** --- To test a polygonalized sphere using a tolerance defined */
	TPZVec<REAL> center(3,-3.);
	REAL radius = 0.5;
	REAL tol = 0.002;
	TPZGeoMesh *ggrid = TPZGenSpecialGrid::GeneratePolygonalSphereFromOctahedron(center,radius,tol);
	TPZCompMesh *cgrid = new TPZCompMesh(ggrid);
	TPZMaterial * mat = new TPZElasticityMaterial(1,1.e5,0.2,0,0);
	cgrid->InsertMaterialObject(mat);
	cgrid->AutoBuild();
	std::cout << "N Elements = " << cgrid->NElements() << std::endl << "N G Elements = " << ggrid->NElements() << std::endl;
	TPZVTKGeoMesh::PrintGMeshVTK(ggrid,fgeom);
	
	/** --- To test a polygonalized sphere using number of refinements */
	int nrefs = 5;
	TPZGeoMesh *ggrid2 = TPZGenSpecialGrid::GeneratePolygonalSphereFromOctahedron(center,radius,nrefs);
	TPZCompMesh *cgrid2 = new TPZCompMesh(ggrid2);
	cgrid2->InsertMaterialObject(mat);
	cgrid2->AutoBuild();
	std::cout << "N Elements = " << cgrid2->NElements() << std::endl << "N G Elements = " << ggrid2->NElements() << std::endl;
	TPZVTKGeoMesh::PrintGMeshVTK(ggrid2,fgeom2);
	
	////  ----  END SPHERE  -----  */
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
#ifdef LOG4CXX
    if (loggerpoint->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Point " << x2;
        LOGPZ_DEBUG(loggerpoint, sout.str())
    }
#endif
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

/** @} */
//void GradientReconstructionByLeastSquares(TPZFMatrix<REAL> &gradients,TPZCompMesh *cmesh,int var,int n_var=0,bool continuous=false);

/**
 * @brief This project shows the creation of a rectangular mesh (two-dimensional) and the creation of a three-dimensional cube mesh using extrude method (ExtendMesh).
 */
int main2D(int argc, char *argv[]) {
	
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
		TPZGeoEl *gel;   // *gel1, *gel2, *gel3;
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
		comp->ExpandSolution();  // Clean connects not connected at least one element enabled.
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
			
			//int neq = comp->NEquations();
			//	an.NEquations();
			//		an.Solution().Print();
			an.Run();
			
			// Computing approximation of gradient
			/** 
			 * @brief Method to reconstruct a gradient after run Solve of the analysis
			 * @param cmesh Computational mesh with solution */
			//	TPZFMatrix<REAL> gradients;
			//	GradientReconstructionByLeastSquares(gradients,comp,0,0,true);
			//	gradients.Print();
			
			// Post processing
			an.PostProcess(1,dim);
			{
				std::ofstream out(MeshFileName.c_str());
				comp->LoadReferences();
				TPZVTKGeoMesh::PrintCMeshVTK(comp->Reference(), out, false);
			}
			
		}
	}
	return 0;
}

/////

/////

//	gRefDBase.InitializeRefPatterns();
// Inserting a special file with refinement pattern 
/* std::string filename = REFPATTERNDIR;
 filename += "/3D_Hexa_Rib_Side_16_16_18_18.rpt";
 //filename += "/3D_Hexa_Face_20.rpt";
 
 TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
 if(!gRefDBase.FindRefPattern(refpat))
 {
 gRefDBase.InsertRefPattern(refpat);
 }
 refpat->InsertPermuted();
 */

// INITIAL POINT FOR SOLVING AND APPLYING REFINEMENT
/*
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
 }
 // Uniform refinement. Two times
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
