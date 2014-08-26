#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzpoisson3d.h"
//#include "pzhybridpoisson.h"
#include "pzpoisson3dreferred.h"
#include "pzelasmat.h"
#include "pzelasthybrid.h"
#include "pzmat1dlin.h"
#include "TPZMatLaplacianLagrange.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzfunction.h"
#include "pzgraphmesh.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "pzgengrid.h"

#include "TPZMHMeshControl.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphysics"));
#endif

using namespace std;


const int matInterno = 1;
const int matCoarse = 2;
const int matInterface = 3;
const int matInterfaceCorase = 4;

const int dirichlet = 0;
const int neumann = 1;
const int mixed = 2;

int const bc1=-1;
int const bc2=-2;
int const bc3=-3;
int const bc4=-4;
int const bc5=-5;


TPZGeoMesh *MalhaGeom(REAL Lx, REAL Ly);
void RefinamentoUniforme(TPZAutoPointer<TPZGeoMesh> gmesh, int nref,TPZVec<int> dims);
void InsertMaterialObjects(TPZCompMesh &cmesh);
void GetElIndexCoarseMesh(TPZAutoPointer<TPZGeoMesh>  gmesh, std::set<long> &coarseindex);

///problema de Steklov
TPZGeoMesh *GMeshSteklov(bool triang_elements);
void RefinamentoSingular(TPZAutoPointer<TPZGeoMesh> gmesh,int nref);
void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
static void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    TPZFMatrix<STATE> du(2,1);
    SolExataSteklov(loc,result,du);
}


int main(int argc, char *argv[])
{	
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    
    TPZAutoPointer<TPZGeoMesh> gmesh = GMeshSteklov(false);
	ofstream arg0("gmesh0.txt");
	gmesh->Print(arg0);
    
    
    //-------- construindo malha coarse ----------
    
    //1 refinamento uniforme
    TPZVec<int> dims(2,0);
    dims[0]=1; dims[1]=2;
    int nref = 1;
    RefinamentoUniforme(gmesh, nref, dims);
    
    nref = 4;
    RefinamentoSingular(gmesh, nref);
    
    
    ofstream arg1("gmesh1.txt");
	gmesh->Print(arg1);
    //index dos elementos da malha coarse
    std::set<long> coarseindex;
    GetElIndexCoarseMesh(gmesh, coarseindex);
    
    TPZAutoPointer<TPZGeoMesh> gmesh2 = new TPZGeoMesh(gmesh);
    
    dims.Resize(1, 0);
    dims[0]=2;
    nref = 2;
    RefinamentoUniforme(gmesh2, nref, dims);
	
	return EXIT_SUCCESS;
}

TPZGeoMesh *MalhaGeom(REAL Lx, REAL Ly)
{
    int Qnodes = 4;
	long dim = 2;
    
	TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <long> TopolQuad(4);
	TPZVec <long> TopolLine(2);
	
	//indice dos nos
	long id = 0;
	REAL valx;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = Lx - xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,Ly);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
	//indice dos elementos
	id = 0;
    
    //elementos internos
    TopolQuad[0] = 0;
	TopolQuad[1] = 1;
	TopolQuad[2] = 2;
	TopolQuad[3] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matInterno,*gmesh);
	id++;
    
    //elementos de contorno
	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
	id++;
	
	TopolLine[0] = 1;
	TopolLine[1] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
	id++;
	
	TopolLine[0] = 2;
	TopolLine[1] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
	id++;
	
	TopolLine[0] = 3;
	TopolLine[1] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
	id++;
    
    //construir a malha
	gmesh->BuildConnectivity();
	
	return gmesh;
}

TPZGeoMesh *GMeshSteklov(bool triang_elements)
{
    TPZManVector<int,2> nx(2,2);
    REAL extent = 1;
    nx[1] =1;
    TPZManVector<REAL,3> x0(3,0.),x1(3,extent);
    x0[0] = -extent;
    TPZGenGrid gengrid(nx,x0,x1);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    if(triang_elements)
    {
        gengrid.SetElementType(1);
    }
    gengrid.Read(gmesh);
    
    //elementos de contorno
    TPZManVector<REAL,3> firstpoint(3,0.),secondpoint(3,0.);
    firstpoint[0] = extent;
    secondpoint[0] = extent;
    secondpoint[1] = extent;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc2);
    gengrid.SetBC(gmesh,6,bc3);
    gengrid.SetBC(gmesh,7,bc4);
    firstpoint[0] = -extent;
    secondpoint[0] = 0.;
    secondpoint[1] = 0.;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc5);
    firstpoint = secondpoint;
    secondpoint[0] = extent;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc1);
    
    return gmesh;
}

void RefinamentoUniforme(TPZAutoPointer<TPZGeoMesh> gmesh, int nref,TPZVec<int> dims)
{
    
    int ir, iel, k;
    int nel=0, dim=0;
    int ndims = dims.size();
	for(ir = 0; ir < nref; ir++ )
    {
		TPZVec<TPZGeoEl *> filhos;
        nel = gmesh->NElements();
        
		for (iel = 0; iel < nel; iel++ )
        {
			TPZGeoEl * gel = gmesh->ElementVec()[iel];
            if(!gel) DebugStop();
            
            dim = gel->Dimension();
            
            for(k = 0; k<ndims; k++)
            {
                if(dim == dims[k])
                {
                    gel->Divide (filhos);
                    break;
                }
            }
		}
	}
    
}


void InsertMaterialObjects(TPZCompMesh &cmesh)
{
	/// criar materiais
	int dim = cmesh.Dimension();
    TPZMatLaplacianLagrange *material1 = new TPZMatLaplacianLagrange(matInterno,dim);
    
	TPZMaterial * mat1(material1);
    
    TPZMat1dLin *materialCoarse = new TPZMat1dLin(matCoarse);
    TPZFNMatrix<1,STATE> xk(1,1,0.),xb(1,1,0.),xc(1,1,0.),xf(1,1,0.);
    materialCoarse->SetMaterial(xk, xc, xb, xf);
    
    cmesh.InsertMaterialObject(materialCoarse);
    
    //    REAL diff = -1.;
    //	REAL conv = 0.;
    //	TPZVec<REAL> convdir(3,0.);
    //	REAL flux = 8.;
    //
    //	material1->SetParameters(diff, conv, convdir);
    //	material1->SetInternalFlux( flux);
    
	cmesh.InsertMaterialObject(mat1);
	
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,1.), val2(2,1,0.);
	
    //BC -1
    TPZMaterial * BCondD1 = material1->CreateBC(mat1, bc1,neumann, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet1 = new TPZDummyFunction<REAL>(Dirichlet);
    BCondD1->SetForcingFunction(bcmatDirichlet1);
    cmesh.InsertMaterialObject(BCondD1);
    
    //BC -2
	TPZMaterial * BCondD2 = material1->CreateBC(mat1, bc2,neumann, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet2 = new TPZDummyFunction<REAL>(Dirichlet);
    BCondD2->SetForcingFunction(bcmatDirichlet2);
    cmesh.InsertMaterialObject(BCondD2);
    
    //BC -3
	TPZMaterial * BCondD3 = material1->CreateBC(mat1, bc3,neumann, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet3 = new TPZDummyFunction<REAL>(Dirichlet);
    BCondD3->SetForcingFunction(bcmatDirichlet3);
    cmesh.InsertMaterialObject(BCondD3);
    
    //BC -4
	TPZMaterial * BCondD4 = material1->CreateBC(mat1, bc4,neumann, val1, val2);
    TPZAutoPointer<TPZFunction<REAL> > bcmatDirichlet4 = new TPZDummyFunction<REAL>(Dirichlet);
    BCondD4->SetForcingFunction(bcmatDirichlet4);
    cmesh.InsertMaterialObject(BCondD4);
    
    //BC -5: dirichlet nulo
    TPZMaterial * BCondD5 = material1->CreateBC(mat1, bc5,neumann, val1, val2);
    cmesh.InsertMaterialObject(BCondD5);
}

void GetElIndexCoarseMesh(TPZAutoPointer<TPZGeoMesh>  gmesh, std::set<long> &coarseindex)
{
    int nel = gmesh->NElements();
    int iel;
    int hassubel=0;
    int dim = gmesh->Dimension();
    int eldim;
    for(iel = 0; iel<nel; iel++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        if(!gel) DebugStop();
        
        hassubel = gel->HasSubElement();
        eldim = gel->Dimension();
        if(!hassubel && eldim ==dim)
        {
            coarseindex.insert(gel->Index());
        }
    }
    
}


void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    const REAL n = 0;
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    const REAL r = sqrt(x*x+y*y);
    const REAL t = atan2(y,x);
    const REAL sol = pow((REAL)2,0.25 + n/2.)*pow(r,0.5 + n)*cos((0.5 + n)*t);
    u[0] = sol;
    
    du(0,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(x*cos((0.5 + n)*atan2(y,x)) + y*sin((0.5 + n)*atan2(y,x)));
    du(1,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(y*cos((0.5 + n)*atan2(y,x)) - x*sin((0.5 + n)*atan2(y,x)));
    
}

void RefinamentoSingular(TPZAutoPointer<TPZGeoMesh> gmesh,int nref)
{
    long nnodes = gmesh->NNodes();
    long in;
    for (in=0; in<nnodes; in++) {
        TPZGeoNode *gno = &gmesh->NodeVec()[in];
        if (abs(gno->Coord(0))< 1.e-6 && abs(gno->Coord(1)) < 1.e-6) {
            break;
        }
    }
    if (in == nnodes) {
        DebugStop();
    }
    TPZGeoElSide gelside;
    long nelem = gmesh->NElements();
    for (long el = 0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        int ncorner = gel->NCornerNodes();
        for (int ic=0; ic<ncorner; ic++) {
            long nodeindex = gel->NodeIndex(ic);
            if (nodeindex == in) {
                gelside = TPZGeoElSide(gel, ic);
                break;
            }
        }
        if (gelside.Element()) {
            break;
        }
    }
    if (!gelside.Element()) {
        DebugStop();
    }
    for (int iref = 0; iref <nref; iref++) {
        TPZStack<TPZGeoElSide> gelstack;
        gelstack.Push(gelside);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            gelstack.Push(neighbour);
            neighbour = neighbour.Neighbour();
        }
        long nstack = gelstack.size();
        for (long ist=0; ist < nstack; ist++) {
            if (!gelstack[ist].Element()->HasSubElement()) {
                TPZVec<TPZGeoEl *> subel;
                gelstack[ist].Element()->Divide(subel);
            }
        }
    }
}

