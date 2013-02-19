//
//  File.cpp
//  PZ
//
//  Created by Agnaldo Farias on 7/31/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//


#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzreferredcompel.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzinterpolationspace.h"

#include "pzpoisson3d.h"
#include "pzelasmat.h"
#include "pzmat1dlin.h"
#include "pzelastpressure.h"
#include "pznlfluidstructure2d.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"

#include "pzanalysis.h"
#include "pznonlinanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"
#include "pzreducedspace.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzbuildmultiphysicsmesh.h"

#include "TPZVTKGeoMesh.h"

#include <iostream>
#include <math.h>

#include <fstream>
#include <sstream>

#include "toolstransienttime.h"

#include "tpzchangeel.h"

using namespace std;

int const matId1 =1; //elastic
int const matId2 =2; //pressure

int const bcmixedx=-1;
int const bcneumannzero=-2;
//int const bcneumannzero=-3;
int const bcmixedy=-4; 
int const bcfluxIn=-5; //bc pressure
int const bcfluxOut=-6; //bc pressure

int const dirichlet =0;
int const neumann = 1;
int const mixed =2;
int const freey = 3;
int const freex = 4;

int const  neum_elast =10;
int const  dir_elast =11;
int const  mix_elast = 20;
int const neum_pressure = 21;
int const dir_pressure = 22;


REAL const Pi = 4.*atan(1.);

TPZGeoMesh *GMesh(int nh,REAL Lx, REAL Ly, REAL Lmax, int nfrac);
TPZGeoMesh *GMesh2(int nh, REAL L);

TPZGeoMesh * PlaneMesh(REAL lf, REAL ldom, REAL hdom, REAL lmax);

TPZCompMesh *CMeshElastic(TPZGeoMesh *gmesh, int pOrder);
TPZCompMeshReferred *CMeshReduced(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);
void InsertMultiphysicsMaterials(TPZCompMesh *cmesh);

TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZNLFluidStructure2d * &mymaterial);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void SaidaPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
void SaidaPressao(TPZCompMesh * cmesh);

void MySolve(TPZAnalysis &an, TPZCompMesh *Cmesh);
void PosProcessamento1(TPZAnalysis &an, std::string plotfile);

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.reducedspace.data"));
#endif

//Dados sem adimensionalizar
REAL LxD = 2.;
REAL LyD = 1.;
REAL LD = 1.;
REAL HD = 1.;
REAL ED = 0.3E5;
REAL nu = 0.25;
REAL viscD = 0.001;
REAL signD = 1.E5;
REAL QinD =  1.E-3;
REAL tD = 10.;
REAL deltaTD = 1.;
//

//#define malhaAgnaldo
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    int p = 2;
	//primeira malha
	
	// geometric mesh (initial)
//#ifdef malhaAgnaldo
	TPZGeoMesh * gmesh = GMesh(4,LxD,LyD,LxD/2,1);
//#else
 //   TPZGeoMesh * gmesh = PlaneMesh(LxD/2., LxD, LyD, LxD/8.);
//#endif
    
    //computational mesh elastic
/** >>> Resolvendo um problema modelo de elastica linear para utilizar a solucao 
    como espaco de aproximacao do problema nao linear (acoplado) */
	TPZCompMesh * cmesh_elast = CMeshElastic(gmesh, p);
    TPZAnalysis an1(cmesh_elast);
    MySolve(an1, cmesh_elast);
                
//    string plotfile("saidaSolution_mesh1.vtk");
//    PosProcessamento1(an1, plotfile);
    
    
/** >>> Passando a malha computacional jah processada do problema modelo elastico para a malha do problema elastico que serah acoplado.
    Esta malha que foi passada servirah como espaco de aproximacao da malha do problema elastico que serah acoplado */
    TPZCompMeshReferred * cmesh_referred = CMeshReduced(gmesh, cmesh_elast, p);
    cmesh_referred->ComputeNodElCon();
    TPZFStructMatrix fstr(cmesh_referred);
    TPZFMatrix<STATE> rhs;//(1);
    TPZAutoPointer<TPZMatrix<STATE> > strmat = fstr.CreateAssemble(rhs,NULL);
	
/** >>> Criando a malha computacional para o problema de fluxo de fluido */
    TPZCompMesh * cmesh_pressure = CMeshPressure(gmesh, p);
    
/** >>> Problema multifisico : acoplamento das 2 malhas computacionais anteriores (elastica e fluxo) */
    //multiphysic mesh
    TPZVec<TPZCompMesh *> meshvec(2);
	meshvec[0] = cmesh_referred;
	meshvec[1] = cmesh_pressure;
    
    gmesh->ResetReference();
	TPZNLFluidStructure2d * mymaterial = NULL;
/** >>> Serah utilizado o material criado (pznlfluidstructure2D do Agnaldo) */
    TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh, meshvec, mymaterial);
    mphysics->SetDefaultOrder(p);
    
    REAL deltaT = deltaTD; //second
    mymaterial->SetTimeStep(deltaT);
    REAL maxTime = tD;
    
    
/** >>> Metodo de resolucao de problema transient */
    //TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    TPZAnalysis *an = new TPZAnalysis(mphysics);
    TPZFMatrix<REAL> InitialSolution = an->Solution();
    ToolsTransient::SolveSistTransient(deltaT, maxTime, InitialSolution, an, mymaterial, meshvec, mphysics);
    
    return 0;
}

ofstream outfile1("SaidaPressao1.nb");
void SaidaPressao(TPZCompMesh * cmesh){
    outfile1<<" Saida2 = {";
    for(int i = 0;  i< cmesh->ElementVec().NElements(); i++)
    {
        TPZCompEl * cel =cmesh->ElementVec()[i];
        TPZInterpolatedElement * sp = dynamic_cast <TPZInterpolatedElement*>(cel);
        if(!sp) continue;
        TPZVec<REAL> qsi(1,0.),out(3,0.);
        TPZMaterialData data;
        sp->InitMaterialData(data);
        for(int j = 0; j < 1; j++){
            qsi[0] = -1.+2.*i/10.;
            sp->ComputeShape(qsi, data);
            sp->ComputeSolution(qsi, data);
            TPZVec<REAL> SolP = data.sol[0]; 
            cel->Reference()->X(qsi,out);
            outfile1<<"{" << out[0]<< ", " << data.sol[0]<<"}";
        }
        if(i!= cmesh->ElementVec().NElements()-1) outfile1<<", ";
        if(i== cmesh->ElementVec().NElements()-1) outfile1<<"};"<<std::endl;
    }
    outfile1<<"ListPlot[Saida2,Joined->True]"<<endl;
}

TPZGeoMesh *GMesh(int nh,REAL Lx, REAL Ly, REAL Lmax,int nfrac){
    
    int Qnodes = 4 + 2*nfrac;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolQuad(4);
    TPZVec <int> TopolTriang(3);
	TPZVec <int> TopolLine(2);
    TPZVec <int> TopolPoint(1);
	
	//indice dos nos
	int id = 0;
	REAL valx;
	for(int xi = 0; xi <= nfrac; xi++)
	{
		valx = xi*Lmax/nfrac;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0,valx);//coord X
		Node[id].SetCoord(1,0.);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
    for(int xi = 0; xi <2; xi++){
        
        REAL valy = xi*Ly;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0,Lx);//coord X
        Node[id].SetCoord(1,valy);//coord Y
        gmesh->NodeVec()[id] = Node[id];
        id++;
    }
    
	for(int xi = 0; xi <= nfrac; xi++)
	{
		valx = Lmax - xi*Lmax/nfrac;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0,valx);//coord X
		Node[id].SetCoord(1,Ly);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
	//indice dos elementos
	id = 0;
    
    for(int i=0; i<nfrac; i++){
        TopolQuad[0] = i;
        TopolQuad[1] = i+1;
        TopolQuad[2] = (Qnodes-1)-(i+1);
        TopolQuad[3] = (Qnodes-1)-i;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId1,*gmesh);
        id++;
    }

    TopolQuad[0] = nfrac;
    TopolQuad[1] = nfrac+1;
    TopolQuad[2] = nfrac+2;
    TopolQuad[3] = nfrac+3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId1,*gmesh);
    id++;

    for(int i=0; i<nfrac; i++){
        TopolLine[0] = i;
        TopolLine[1] = i+1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId2,*gmesh);
        id++;
    }
        
    TopolLine[0] = nfrac;
    TopolLine[1] = nfrac+1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcmixedx,*gmesh);
    id++;
    
    
    for(int i=nfrac+1; i<Qnodes-1; i++){
        TopolLine[0] = i;
        TopolLine[1] = i+1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcneumannzero,*gmesh);
        id++;
    }
    
    TopolLine[0] = Qnodes-1;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bcmixedy,*gmesh);
    id++;
    
    TopolPoint[0] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolPoint,bcfluxIn,*gmesh);
    id++;
    
    TopolPoint[0] = nfrac;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolPoint,bcfluxOut,*gmesh);

    
	gmesh->BuildConnectivity();
    
    //refinamento uniforme
    for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			//if (gel->Dimension() == 2) gel->Divide (filhos);
            gel->Divide (filhos);
		}//for i
	}//ref
    
     //gRefDBase.InitializeRefPatterns();
    
//    int nrefdir = 2;
//    std::set<int> matNewm;
//    matNewm.insert(matId2);
//    for(int ii = 0; ii < nrefdir; ii++)
//    {
//        int nels = gmesh->NElements();
//        for(int iel = 0; iel < nels; iel++)
//        {
//            TPZRefPatternTools::RefineDirectional(gmesh->ElementVec()[iel], matNewm);
//        }
//    }

	return gmesh;
}

TPZGeoMesh *GMesh2(int nh, REAL L){
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(1);
	gmesh->NodeVec().Resize(2);
	TPZVec<TPZGeoNode> Node(2);
    TPZVec <int> TopolLine(2);
    TPZVec <int> TopolPoint(1);
	
	//indice dos nos
	int id = 0;
	REAL valx;
	for(int xi = 0; xi < 2; xi++)
	{
		valx = xi*L;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0,valx );//coord X
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}

    id=0;
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId2,*gmesh);
    id++;
    
    TopolPoint[0] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bcfluxIn,*gmesh);
    id++;

    TopolPoint[0] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bcfluxOut,*gmesh);
    
    gmesh->BuildConnectivity();
    
    //refinamento uniforme
    for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			//if (gel->Dimension() == 2) gel->Divide (filhos);
            gel->Divide (filhos);
		}//for i
	}//ref

    
    return gmesh;
}

#define usingRefUnif
TPZGeoMesh * PlaneMesh(REAL lf, REAL ldom, REAL hdom, REAL lmax)
{
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    int ndivfrac = int(lf/lmax + 0.5);
    int ndivoutfrac = int((ldom - lf)/lmax + 0.5);
    int ndivh = int(hdom/lmax + 0.5);
    
    int ncols = ndivfrac + ndivoutfrac + 1;
    int nrows = ndivh + 1;
    int nnodes = nrows*ncols;
    
    gmesh->NodeVec().Resize(nnodes);
    
    REAL deltadivfrac = lf/ndivfrac;
    REAL deltadivoutfrac = (ldom-lf)/ndivoutfrac;
    REAL deltandivh = hdom/ndivh;
    
    int nid = 0;
    for(int r = 0; r < nrows; r++)
    {
        for(int c = 0; c < ncols; c++)
        {
            REAL x, y;
            if(c <= ndivfrac)
            {
                x = c*deltadivfrac;
            }
            else
            {
                x = ndivfrac*deltadivfrac + (c-ndivfrac)*deltadivoutfrac;
            }
            y = r*deltandivh;
            
            TPZVec<REAL> coord(3,0.);
            coord[0] = x;
            coord[1] = y;
            gmesh->NodeVec()[r*ncols + c].SetCoord(coord);
            gmesh->NodeVec()[r*ncols + c].SetNodeId(nid);
            nid++;
        }
    }
    
    TPZGeoEl * gel = NULL;
    TPZVec<int> topol(4);
    int indx = 0;
    for(int r = 0; r < nrows-1; r++)
    {
        for(int c = 0; c < ncols-1; c++)
        {
            topol[0] = r*(ncols) + c;
            topol[1] = r*(ncols) + c + 1;
            topol[2] = r*(ncols) + c + 1 + ncols;
            topol[3] = r*(ncols) + c + ncols;
            
            gel = gmesh->CreateGeoElement(EQuadrilateral, topol, matId1, indx);
            gel->SetId(indx);
            indx++;
        }
    }
    
    gmesh->BuildConnectivity();
    
    int nelem = gmesh->NElements();
    for(int el = 0; el < nelem; el++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[el];
        
        //south BC
        TPZGeoElSide sideS(gel,4);
        TPZGeoElSide neighS(sideS.Neighbour());
        if(sideS == neighS)
        {
            if(el < ndivfrac)
            {
                gel->CreateBCGeoEl(4, matId2);
            }
            else
            {
                gel->CreateBCGeoEl(4, bcmixedx);
            }
        }
        
        //east BC
        TPZGeoElSide sideE(gel,5);
        TPZGeoElSide neighE(sideE.Neighbour());
        if(sideE == neighE)
        {
            gel->CreateBCGeoEl(5, bcneumannzero);
        }
        
        //north BC
        TPZGeoElSide sideN(gel,6);
        TPZGeoElSide neighN(sideN.Neighbour());
        if(sideN == neighN)
        {
            gel->CreateBCGeoEl(6, bcneumannzero);
        }
        
        //west BC
        TPZGeoElSide sideW(gel,7);
        TPZGeoElSide neighW(sideW.Neighbour());
        if(sideW == neighW)
        {
            gel->CreateBCGeoEl(7, bcmixedy);
        }
    }
    
    topol.Resize(1);
    topol[0] = 0;
    gel = gmesh->CreateGeoElement(EPoint, topol, bcfluxIn, indx);
    indx++;
    
    topol[0] = ndivfrac;
    gel = gmesh->CreateGeoElement(EPoint, topol, bcfluxOut, indx);
    
    gmesh->BuildConnectivity();
    
#ifdef usingQPoints
    TPZGeoElSide pt(gel,0);
    TPZGeoElSide ptneigh(pt.Neighbour());
    while(pt != ptneigh)
    {
        if(ptneigh.Element()->HasSubElement() == false)
        {
            int neighSide = ptneigh.Side();
            TPZGeoEl * ptneighEl = TPZChangeEl::ChangeToQuarterPoint(gmesh, ptneigh.Element()->Id(), neighSide);
            ptneigh = ptneighEl->Neighbour(neighSide);
        }
        else
        {
            ptneigh = ptneigh.Neighbour();
        }
    }
#endif
    
#ifdef usingRefUnif
    int nrefUnif = 2;
    for(int ref = 0; ref < nrefUnif; ref++)
    {
        nelem = gmesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            if(gmesh->ElementVec()[el]->Dimension() < 1) continue;
            if(gmesh->ElementVec()[el]->HasSubElement()) continue;
            if(gmesh->ElementVec()[el]->MaterialId() == matId2)
            {
                TPZVec<TPZGeoEl*> sons;
                gmesh->ElementVec()[el]->Divide(sons);
                continue;
            }
            for(int s = 0; s < gmesh->ElementVec()[el]->NSides(); s++)
            {
                TPZGeoElSide gelside(gmesh->ElementVec()[el],s);
                TPZGeoElSide neighside(gelside.Neighbour());
                bool refinedAlready = false;
                while(neighside != gelside)
                {
                    if(neighside.Element()->MaterialId() == matId2)
                    {
                        TPZVec<TPZGeoEl*> sons;
                        gmesh->ElementVec()[el]->Divide(sons);
                        refinedAlready = true;
                        break;
                    }
                    neighside = neighside.Neighbour();
                }
                if(refinedAlready == true)
                {
                    break;
                }
            }
        }
    }
#endif
    
#ifdef usingRefdir
    std::set<int> matDir;
    //matDir.insert(__2DfractureMat_inside);
    matDir.insert(bcfluxOut);
    int nrefDir = 1;
    for(int ref = 0; ref < nrefDir; ref++)
    {
        nelem = gmesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            if(!gmesh->ElementVec()[el]) continue;
            if(gmesh->ElementVec()[el]->Dimension() < 1) continue;
            if(gmesh->ElementVec()[el]->HasSubElement()) continue;
            TPZRefPatternTools::RefineDirectional(gmesh->ElementVec()[el], matDir);
        }
    }
#endif
    
    
    return gmesh;
}

TPZCompMesh *CMeshElastic(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
	
    TPZVec<REAL> force(dim,0.);
    REAL E = ED;
	REAL poisson = nu;
    //int planestress = 1;
    int planestrain = 0;
    
    TPZElasticityMaterial *material;
	material = new TPZElasticityMaterial(matId1, E, poisson, force[0], force[1], planestrain);

    TPZMaterial * mat(material);
    material->NStateVariables();
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    REAL big = material->gBigNumber;
    REAL sign = signD;
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    val2(0,0)=0;
    val2(1,0)=sign;
    TPZMaterial * BCond1 = material->CreateBC(mat, matId2,neumann, val1, val2);
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    val1(1,1) = big;
    TPZMaterial * BCond2 = material->CreateBC(mat, bcmixedx,mixed, val1, val2);
    
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    TPZMaterial * BCond3 = material->CreateBC(mat, bcneumannzero,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, bcneumannzero,dirichlet, val1, val2);
     
    val1.Redim(2,2);
    val2.Redim(2,1);
    val1(0,0) = big;
    TPZMaterial * BCond5 = material->CreateBC(mat, bcmixedy,mixed, val1, val2);
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->InsertMaterialObject(BCond1);
	cmesh->InsertMaterialObject(BCond2);
	cmesh->InsertMaterialObject(BCond3);
	cmesh->InsertMaterialObject(BCond4);
	cmesh->InsertMaterialObject(BCond5);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
        
	return cmesh;
}

TPZCompMeshReferred *CMeshReduced(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int pOrder){
    
    /// criar materiais
	int dim = 2;
    
    TPZVec<REAL> force(dim,0.);
    REAL E = ED;
	REAL poisson = nu;
    //int planestress = 1;
    int planestrain = 0;
    
    TPZElasticityMaterial *material;
	material = new TPZElasticityMaterial(matId1, E, poisson, force[0], force[1], planestrain); 
	material->NStateVariables();
    
    
    TPZCompMeshReferred *cmeshreferred = new TPZCompMeshReferred(gmesh);
    
    cmeshreferred->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmeshreferred->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    REAL big = material->gBigNumber;
    REAL sign = signD;
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    val2(0,0)=0;
    val2(1,0)=sign;
    TPZMaterial * BCond1 = material->CreateBC(mat, matId2,neumann, val1, val2);
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    val1(1,1) = big;
    TPZMaterial * BCond2 = material->CreateBC(mat, bcmixedx,mixed, val1, val2);
    
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    TPZMaterial * BCond3 = material->CreateBC(mat, bcneumannzero,dirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, bcneumannzero,dirichlet, val1, val2);
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    val1(0,0) = big;
    TPZMaterial * BCond5 = material->CreateBC(mat, bcmixedy,mixed, val1, val2);

    
    int numsol = cmesh->Solution().Cols();
    cmeshreferred->AllocateNewConnect(numsol, 1, 1);
    
	TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmeshreferred);
    
    cmeshreferred->InsertMaterialObject(BCond1);
    cmeshreferred->InsertMaterialObject(BCond2);
    cmeshreferred->InsertMaterialObject(BCond3);
    cmeshreferred->InsertMaterialObject(BCond4);
    cmeshreferred->InsertMaterialObject(BCond5);
    
	cmeshreferred->SetDefaultOrder(pOrder);
    cmeshreferred->SetDimModel(dim);
	
    gmesh->ResetReference();
	//Ajuste da estrutura de dados computacional
	cmeshreferred->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    cmeshreferred->LoadReferred(cmesh);
   
    return cmeshreferred;
}

TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder){
    
    /// criar materiais
	int dim = 1;
	
    TPZMat1dLin *material;
	material = new TPZMat1dLin(matId2);
    
    TPZFMatrix<REAL> xk(1,1,1.);
    TPZFMatrix<REAL> xc(1,1,0.);
    TPZFMatrix<REAL> xb(1,1,0.);
    TPZFMatrix<REAL> xf(1,1,-2.);
    material->SetMaterial(xk,xc,xb,xf);
    
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    REAL vazao = -1.;
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    val2(0,0)=vazao;
    TPZMaterial * BCond1 = material->CreateBC(mat, bcfluxIn, neumann, val1, val2);
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    //TPZMaterial * BCond2 = material->CreateBC(mat, matId2, neumann, val1, val2);
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    REAL pressao = 0.75;
    val2(0,0)=pressao;
    TPZMaterial * BCond3 = material->CreateBC(mat, bcfluxOut,dirichlet, val1, val2);
    
    cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
    cmesh->InsertMaterialObject(BCond1);
   // cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;    
}

void MySolve(TPZAnalysis &an, TPZCompMesh *Cmesh)
{			
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(Cmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	//Saida de Dados: solucao e  grafico no VT
	ofstream file("Solutout");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

void PosProcessamento1(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(4), vecnames(1);
	vecnames[0] = "displacement";
    scalnames[0] = "SigmaX";
	scalnames[1] = "SigmaY";
   scalnames[2] = "sig_x"; 
    scalnames[3] = "sig_y";
    
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZNLFluidStructure2d * &mymaterial){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int MatId=1;
    int dim =2;
    mymaterial = new TPZNLFluidStructure2d(MatId,dim);
    
    //data elasticity
    REAL fx = 0.;
    REAL fy = 0.;
    REAL E = ED;
	REAL poisson = nu;
    //int planestress = 1;
    int planestrain = 0;
    mymaterial->SetElasticParameters(E, poisson, fx, fy);
    mymaterial->SetfPlaneProblem(planestrain);

    //data pressure
    //TPZFMatrix<REAL> xk(1,1,1.);
    //TPZFMatrix<REAL> xf(1,1,0.);
    REAL hw = HD;
    REAL xvisc = viscD;
    REAL xql = 0.;
    mymaterial->SetParameters(hw, xvisc, xql);
    
    
    TPZMaterial *mat(mymaterial);
    mphysics->InsertMaterialObject(mat);
    
    mymaterial->NStateVariables();
    ///Inserir condicao de contorno
    REAL big = mymaterial->gBigNumber;
    REAL sign = signD;
    
    TPZFMatrix<REAL> val1(3,2,0.), val2(3,1,0.);
    val2(0,0)=0.;
    val2(1,0)=sign;
    TPZMaterial * BCond1 = mymaterial->CreateBC(mat, matId2,neumann, val1, val2);
    
    val2.Redim(3,1);
    val1.Redim(3,2);
    val1(1,1) = big;
    TPZMaterial * BCond2 = mymaterial->CreateBC(mat, bcmixedx,mix_elast, val1, val2);
    
    val2.Redim(3,1);
    val1.Redim(3,2);
    TPZMaterial * BCond3 = mymaterial->CreateBC(mat, bcneumannzero,dir_elast, val1, val2);
    TPZMaterial * BCond4 = mymaterial->CreateBC(mat, bcneumannzero,dir_elast, val1, val2);
    
    val2.Redim(3,1);
    val1.Redim(3,2);
    val1(0,0) = big;
    TPZMaterial * BCond5 = mymaterial->CreateBC(mat, bcmixedy, mix_elast, val1, val2);
    
    REAL vazao = QinD;
    val2.Redim(3,1);
    val1.Redim(3,2);
    val2(2, 0)=vazao;
    TPZMaterial * BCond6 = mymaterial->CreateBC(mat, bcfluxIn, neum_pressure, val1, val2);
    
    
    val2.Redim(3,1);
    val1.Redim(3,2);
    val2(2, 0)=0.0;
    TPZMaterial * BCond7 = mymaterial->CreateBC(mat, bcfluxOut, neum_pressure, val1, val2);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    mphysics->InsertMaterialObject(BCond5);
    mphysics->InsertMaterialObject(BCond6);
    mphysics->InsertMaterialObject(BCond7);

    
    mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    //    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional Multiphysic\n ";
    //        mphysics->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str());
    //	}
    //#endif
    //    
    return mphysics;
}


void InsertMultiphysicsMaterials(TPZCompMesh *cmesh)
{
    TPZMat1dLin *material;
	material = new TPZMat1dLin(matId2);
    
    TPZFMatrix<REAL> xk(1,1,1.);
    TPZFMatrix<REAL> xc(1,1,0.);
    TPZFMatrix<REAL> xb(1,1,0.);
    TPZFMatrix<REAL> xf(1,1,0.);
    material->SetMaterial(xk,xc,xb,xf);
    
    TPZMaterial * mat(material);
    material->NStateVariables();
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    REAL vazao = 0.2;
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    val2(0,0)=vazao;
    TPZMaterial * BCond1 = material->CreateBC(mat, bcneumannzero,neumann, val1, val2);

    cmesh->InsertMaterialObject(BCond1);
    
    /// criar materiais
	int dim = 2;
    
    TPZVec<REAL> force(dim,0.);
    REAL E = 1.43996*1e10;
	REAL poisson = 0.2;
    int planestress = 0;
    
	TPZElasticityMaterial *matelas;
	matelas = new TPZElasticityMaterial(matId1, E, poisson, force[0], force[1], planestress); 
    cmesh->InsertMaterialObject(matelas);
    
    ///Inserir condicao de contorno
    REAL big = material->gBigNumber;
    REAL sign = 1.;
    
    TPZFMatrix<REAL> val11(2,2,0.), val21(2,1,0.);
    val21(0,0)=0;
    val21(1,0)=sign;
    TPZMaterial * BCond4 = material->CreateBC(mat, matId2, neumann, val11, val21);
    
    TPZFMatrix<REAL> val12(2,2,0.), val22(2,1,0.);
    val12(1,1) = big;
    TPZMaterial * BCond2 = material->CreateBC(mat, bcmixedx,mixed, val12, val22);
    
    TPZFMatrix<REAL> val13(2,2,0.), val23(2,1,0.);
    val13(0,0) = big;
    TPZMaterial * BCond3 = material->CreateBC(mat, bcneumannzero,dirichlet, val13, val23);
    
    
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);

}

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(2), vecnames(0);
	scalnames[0]  = "Pressure";
    scalnames[1]  = "MinusKGradP";
    
	
	const int dim = 2;
	int div =0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
    
}
