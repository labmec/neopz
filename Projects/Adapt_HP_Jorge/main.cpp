#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZVTKGeoMesh.h"

#include "MultiResMesh.h"
#include "pzlog.h"
#include "pzbstrmatrix.h"

#include <stdio.h>

#include "pzexplfinvolanal.h"
#include "adapt.h"

using namespace std;

// Global variable
int gLMax;

void InitializeSolver(TPZAnalysis &an);
void InitialSolutionLinearConvection(TPZFMatrix<REAL> &InitialSol, TPZCompMesh *cmesh);
void PrintGeoMeshVTKWithDimensionAsData(TPZGeoMesh *gmesh,char *filename);

/** 
 * @brief Make uniform refinement of the geometric mesh. This method finishs making the new connectivities for the elements (ResetConnectivities() and BuildConnectivity()) 
 * @param nDiv Number of divisions wished
 * @param gmesh Geometric mesh contained the elements
 * @param dim Dimension of the geometric elements will be divided
 * @param allmaterial If true divides all geometric element, if false check the material match with allmaterial to divide the element
 * @param matidtodivided To divide only elements with this material id
 */
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial=true, const int matidtodivided=0);

int main(int argc, char *argv[]) {
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
#ifdef DEBUG
	char saida[260];
#endif

	cout << "Generating geometric mesh bi-dimensional ...\n";
    // First rectangular mesh:
	// The rectangular mesh has four corners: (0,0,0), (1,0,0), (1,1,0) and (0,1,0)
	// and was divides in two segments on X and two on Y, then hx = 0.5 and hy = 0.5
	// Has 4 elements, 9 connects and 8 bc elements
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;

	TPZManVector<REAL> x0(3,0.), x1(3,0.);  // Corners of the rectangular mesh. Coordinates of the first extreme are zeros.
	x1[0] = 0.6; x1[1] = 0.2;    // coordinates of the upper right extreme of the rectangular mesh: (0.,0.,0.) (0.6,0.2,0.0)

	TPZManVector<int> nx(3,2);   // subdivisions in X and in Y. 
	TPZGenGrid gen(nx,x0,x1);    // mesh generator. On X we has three segments and on Y two segments. Then: hx = 0.2 and hy = 0.1  
	gen.SetElementType(0);       // type = 0 means rectangular elements
	gen.Read(gmesh);             // generating grid in gmesh

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
	PrintGeoMeshVTKWithDimensionAsData(gmesh.operator->(),saida);
	sprintf(saida,"meshes.vtk");
	PrintGeoMeshVTKWithDimensionAsData(gmesh2.operator->(),saida);
#endif

	// Extending geometric mesh (two-dimensional) to three-dimensional geometric mesh
	// The elements are hexaedras(cubes) over the quadrilateral two-dimensional elements
	TPZExtendGridDimension gmeshextend(gmesh2,0.3);
	
	TPZGeoMesh *gmesh3D = gmeshextend.ExtendedMesh(2,-5,-6);
#ifdef DEBUG
	sprintf(saida,"meshextruded.vtk");
	PrintGeoMeshVTKWithDimensionAsData(gmesh3D,saida);
#endif
	
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
/*
void PrintGeoMeshVTKWithData(TPZGeoMesh *gmesh,char *filename,int var) {
	int nrows = data.Rows();
	int ncols = data.Cols();
	int i, j, nvars;
	std::map<int, TPZAutoPointer <TPZMaterial> >::iterator m;
	for(m=gmesh->Reference()->MaterialVec().begin();m != gmesh->Reference()->MaterialVec().end();m++) {
		if(m->second->Id() > 0) {
			nvars = m->second->NStateVariables();
			break;
		}
	}
	if(m == gmesh->Reference()->MaterialVec().end() || var > nvars) {
		DebugStop();
	}
	TPZChunkVector<double> DataElement;
	int k = 0;
	char newname[1024];
	memset(newname,0,sizeof(newname));
	sprintf(newname,"%s",filename);
	char *p;
	DataElement.Resize(nrows/nvars);
	for(i=0;i<ncols;i++) {
		for(j=0;j<nrows;j+=nvars)
			DataElement[k++] = data.GetVal(j+var,i);

		p = strrchr(newname,'_');
		if(p) {
			sprintf(p,"%i.vtk",i);
		}
		else if(i) {
			p = strrchr(newname,'.');
			sprintf(newname,"_%i.vtk",i);
		}
		// Printing geometric mesh to visualization in Paraview
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh, newname, DataElement);
	}
}	*/

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
