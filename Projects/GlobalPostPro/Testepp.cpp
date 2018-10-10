/**
 * @file
 * @brief Implements a flux calculus for a Poisson problem, using a variational formulation in the post-processing
 */

#include <pzvec.h>
#include <pzgmesh.h>
#include <pzcompel.h>
#include <pzgeoel.h>
#include <pzquad.h>
#include <pzmat2dlin.h>
#include <TPZGeoElement.h>
#include <pzskylstrmatrix.h>
#include <pzcmesh.h>
#include "TPZStream.h"
#include "mat2dpospro.h"
#include "pzbndcond.h"
#include <TPZMatLaplacian.h>
//#include "tpzdifureac.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "pzanalysis.h"
#include "pzstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzbndmat.h"
#include "pzstepsolver.h"
#include "TPZVTKGeoMesh.h"
#include "pzfstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "tpzhierarquicalgrid.h"
#include "pzfunction.h"

void ParametricfunctionX(const TPZVec<REAL> &par, TPZVec<REAL> &X);
void ParametricfunctionY(const TPZVec<REAL> &par, TPZVec<REAL> &X);

#include <iostream>
#include <fstream>
using namespace std;

// cria a malha geométrica, dados os comprimentos em x e em y, o argumento bool é para escolher a malha com triangulos ou quadrilateros
TPZGeoMesh * GetMesh(REAL Lx,REAL Ly, bool triang_elements);
TPZGeoMesh * QuadDomain(int n, REAL t, REAL dt);

//nDiv indica o número de divisões feitas em cada direção
void UniformRefine(TPZGeoMesh* gmesh, int nDiv);

//cria a malha computational
TPZCompMesh * ComputationalMesh(TPZGeoMesh * gmesh, int pOrder);
void SolveSystem(TPZAnalysis &an, TPZCompMesh *fCmesh, bool symmetric_matrix);

//pos processamento da solucao
void OutSolution(TPZAnalysis &an, std::string plotfile);
void SourceTerm(const TPZVec<REAL> &X, TPZVec<STATE> &Result);
void U(const TPZVec<REAL> &X, REAL time, TPZVec<STATE> &Result, TPZFMatrix<STATE> &GradU);

void ExactSigma(const TPZVec<REAL> &X, TPZVec<STATE> &Result, TPZFMatrix<STATE> &Derivative);
//os ids das bc se colocam com iinteiros negativos
int matId = 1;


int dirichlet = 0;
int neumann = 1;
REAL const Pi = 4.*atan(1.);

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.Margui"));
#endif

int main() {
    
#ifdef LOG4CXX
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    FileName = dirname + "/Projects/GlobalPostPro/";
    FileName += "MarguiLog.cfg";
    InitializePZLOG(FileName);
#endif  
    
    REAL t=0.0;
    REAL dt= 0.1;
    int n = 10;
    TPZGeoMesh * gmesh = QuadDomain(n, t, dt);
    
    int porder = 3;
    TPZCompMesh * cmesh = ComputationalMesh(gmesh,porder);
    
    TPZAnalysis *an = new TPZAnalysis(cmesh);
    bool sym = false;
    SolveSystem(*an, cmesh, sym);
    
    string Margui("Result.vtk");
    OutSolution(*an, Margui);
    return 0;

}

TPZGeoMesh* QuadDomain(int n, REAL t, REAL dt)
{
     // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh1 = new TPZGeoMesh;
    GeoMesh1->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh1->NodeVec()[0]=Node;
    
    TPZVec<int64_t> Topology(1,0);
    int elid=0;
    int matid=1;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh1);
    GeoMesh1->BuildConnectivity();
    GeoMesh1->SetDimension(0);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew1.txt");
        GeoMesh1->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew1.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh1,Dummyfile, true);
    }
    
    
    TPZHierarquicalGrid CreateGridFrom(GeoMesh1);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc = new TPZDummyFunction<REAL>(ParametricfunctionX, 5);
    CreateGridFrom.SetParametricFunction(ParFunc);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh2 = CreateGridFrom.ComputeExtrusion(t, dt, n);
    
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew2.txt");
        GeoMesh2->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew2.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh2,Dummyfile, true);
    }
    
    
    
    TPZHierarquicalGrid CreateGridFrom2(GeoMesh2);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc2 = new TPZDummyFunction<REAL>(ParametricfunctionY, 5);
    CreateGridFrom2.SetParametricFunction(ParFunc2);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh3 = CreateGridFrom2.ComputeExtrusion(t, dt, n);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew3.txt");
        GeoMesh3->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew3.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh3,Dummyfile, true);
    }

    return GeoMesh3;
}

void ParametricfunctionX(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = par[0];
    X[1] = 0.0;
    X[2] = 0.0;
}

void ParametricfunctionY(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0;
    X[1] = par[0];
    X[2] = 0.0;
}

TPZGeoMesh *GetMesh (REAL Lx,REAL Ly, bool triang_elements){
  
    return 0;
// 	int Qnodes = 4;
// 	
// 	TPZGeoMesh * gmesh = new TPZGeoMesh;
// 	gmesh->SetMaxNodeId(Qnodes-1);
// 	gmesh->NodeVec().Resize(Qnodes);
// 	TPZVec<TPZGeoNode> Node(Qnodes);
// 	  
// 	TPZVec <int64_t> TopolQuad(4);
//        TPZVec <int64_t> TopolTriang(3);
// 	TPZVec <int64_t> TopolLine(2);
// 	
// 	//indice dos nos
// 	int64_t id = 0;
// 	REAL valx, valy;
// 	REAL x00 = -1., y00 = -1.; //Valores mínimos en cada coordenada del domínio
//     
// 	for(int xi = 0; xi < Qnodes/2; xi++)
// 	{
// 		//valx = -1. + xi*Lx; //quadrilatero [-1,1] x [-1,1]
// 		//valy = -1.; 
// 		valx = x00 + xi*Lx;
// 		valy = y00;		
// 		Node[id].SetNodeId(id);
// 		Node[id].SetCoord(0,valx);//coord X
// 		Node[id].SetCoord(1,valy);//coord Y
// 		gmesh->NodeVec()[id] = Node[id];
// 		id++;
// 	}
// 	
// 	for(int xi = 0; xi < Qnodes/2; xi++)
// 	{
// 		//valx = Lx - xi*Lx - 1;
// 		//valy = -1. + Ly;
// 		valx = Lx - xi*Lx +x00;
// 		valy = y00+Ly;
// 		Node[id].SetNodeId(id);
// 		Node[id].SetCoord(0,valx);//coord X
// 		Node[id].SetCoord(1,valy);//coord Y
// 		gmesh->NodeVec()[id] = Node[id];
// 		id++;
// 	}
// 	
// 	//indice dos elementos
// 	id = 0;
//     
//     if(triang_elements==true)
//     {
//         TopolTriang[0] = 0;
//         TopolTriang[1] = 1;
//         TopolTriang[2] = 2;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
//         id++;
//         
//         TopolTriang[0] = 2;
//         TopolTriang[1] = 3;
//         TopolTriang[2] = 0;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
//         id++;
//         
//         TopolLine[0] = 1;
//         TopolLine[1] = 0;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
//         id++;
//         
//         TopolLine[0] = 1;
//         TopolLine[1] = 2;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
//         id++;
//         
//         TopolLine[0] = 2;
//         TopolLine[1] = 3;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
//         id++;
//         
//         TopolLine[0] = 3;
//         TopolLine[1] = 0;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
//     }
//     else{
//         
//         TopolQuad[0] = 0;
//         TopolQuad[1] = 1;
//         TopolQuad[2] = 2;
//         TopolQuad[3] = 3;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
//         id++;
//         
// 	//elementos de contorno
//         TopolLine[0] = 0;
//         TopolLine[1] = 1;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
//         id++;
//         
//         TopolLine[0] = 1;
//         TopolLine[1] = 2;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
//         id++;
//         
//         TopolLine[0] = 2;
//         TopolLine[1] = 3;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
//         id++;
//         
//         TopolLine[0] = 3;
//         TopolLine[1] = 0;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
//     }
//     
// 	gmesh->BuildConnectivity();
// 	return gmesh;
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

TPZCompMesh *ComputationalMesh(TPZGeoMesh * gmesh, int pOrder)
{
    /// criar materiais
    int dim = gmesh->Dimension();

    // 	Material being used
    Mat2Dpospro * material = new Mat2Dpospro(matId,dim);
    TPZMaterial * mat(material);

    TPZCompEl::SetgOrder(pOrder);
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
// Setting up the explicit functions
     TPZAutoPointer<TPZFunction<STATE> > spatialf;
     TPZAutoPointer<TPZFunction<STATE> > spatialU;     
     TPZDummyFunction<STATE> * dum1 = new TPZDummyFunction<STATE>(SourceTerm, 5);
     TPZDummyFunction<STATE> * dum2 = new TPZDummyFunction<STATE>(U, 5);
//      dum1->SetPolynomialOrder(20);
//      dum2->SetPolynomialOrder(20);     
     spatialf = dum1;
     spatialU = dum2;
     material->SetForcingFunction(spatialf);
     material->SetTimeDependentForcingFunction(spatialU);
     REAL alpha = 1.0;
     REAL delta = 1.0;
     material->SetParameters(alpha,delta);
     
// Setting up the exact function     
    TPZAutoPointer<TPZFunction<STATE> > SigmaExact;
    SigmaExact = new TPZDummyFunction<STATE>(ExactSigma, 5);
    material->SetForcingFunctionExact(SigmaExact);    

    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->InsertMaterialObject(mat);    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

void SolveSystem(TPZAnalysis &an, TPZCompMesh *fCmesh, bool symmetric_matrix)
{
    if(symmetric_matrix){
        TPZSkylineStructMatrix skmat(fCmesh);
        an.SetStructuralMatrix(skmat);
        TPZStepSolver<STATE> direct;
        direct.SetDirect(ELDLt);
        an.SetSolver(direct);
    }
    else{
        TPZSkylineNSymStructMatrix sknmat(fCmesh);
        an.SetStructuralMatrix(sknmat);
        TPZStepSolver<STATE> direct;
        direct.SetDirect(ELU);
        an.SetSolver(direct);
    }
    an.Run();

}

void OutSolution(TPZAnalysis &an, std::string plotfile){
    
	TPZManVector<std::string,10> scalnames(2), vecnames(3);
	scalnames[0] = "U";
	scalnames[1] = "F";
	vecnames[0]= "Sigma";
	vecnames[1]= "SigmaExact";
	vecnames[2]= "GradU";	
	
	const int dim = an.Mesh()->Dimension();
	int div = 3;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

void SourceTerm(const TPZVec<REAL> &X, TPZVec<STATE> &Result){ //função do lado direito
    
    double x = X[0];
    double y = X[1];
    Result[0] = 2*Pi*Pi*sin(Pi*x)*sin(Pi*y);
    
}
void U(const TPZVec<REAL> &X, REAL time, TPZVec<STATE> &Result, TPZFMatrix<STATE> &GradU){
    
    double x = X[0];
    double y = X[1];
    
    if( GradU.Rows() != 2 )
    {
      cout << " Needs 2D Grad " << endl;
      DebugStop();
    }
    
    Result[0]= sin(Pi*x)*sin(Pi*y);
    GradU(0,0) = Pi*cos(Pi*x)*sin(Pi*y);// X direction
    GradU(1,0) = Pi*sin(Pi*x)*cos(Pi*y);// Y direction
}

void ExactSigma(const TPZVec<REAL> &X, TPZVec<STATE> &Result, TPZFMatrix<STATE> &Derivative) {

    double x = X[0];
    double y = X[1];

    if( Result.size() != 2 )
    {
      cout << " Needs 2D Grad " << endl;
      DebugStop();
    }
    
    Result[0]= Pi*cos(Pi*x)*sin(Pi*y);
    Result[1]= Pi*sin(Pi*x)*cos(Pi*y);

}
