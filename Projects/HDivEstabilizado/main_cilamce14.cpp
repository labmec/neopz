
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPBSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"

#include "pzanalysis.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"

#include "pzanalysis.h"

#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

#include "pzhdivfull.h"
#include "pzaxestools.h"
#include "TPZCopySolve.h"
#include "pzstrmatrix.h"

#include <iostream>
#include <math.h>
using namespace std;

int MatId =1;

int bcdirichlet = 0;
int bcneumann = 1;

int BC0=-1;
int BC1=-2;
int BC2=-3;
int BC3=-4;
int BC4=-5;
int BC5=-6;


//meio homogeneo
TPZGeoMesh *GMesh2(REAL Lx, REAL Ly,bool triang_elements);

//meio heterogeneo
TPZGeoMesh *GMesh3(bool triang_elements);

TPZCompMesh *CMeshFlux2(int pOrder, TPZGeoMesh *gmesh);
TPZCompMesh *CMeshPressure2(int pOrder,TPZGeoMesh *gmesh);
TPZCompMesh *CMeshMixed2(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh);

void RefinamentoUnif(TPZGeoMesh* gmesh, int nDiv);
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads, bool direct);
void PosProcessMultph(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void PosProcessFluxo(TPZAnalysis &an, std::string plotfile);

//solucao exata meio homogeneo
void SolExataMista(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void SolExataPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);

//lado direito da equacao - meio homogeneo
void ForcingMista(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//Para condicao de contorno - meio homogeneo
void NeumannAcima(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void NeumannAbaixo(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void DirichletEsquerda(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//solucoes exata meio heterogeneo
void SolFluxoHeter(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void SolPressaoHeter(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
void ForcingHeter(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

void NeumannAcimaXMaiorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void NeumannAcimaXMenorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void NeumannAbaixoXMaiorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void NeumannAbaixoXMenorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void DirichletXIgualMenosUm(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void DirichletXIgualUm(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void PermeabilityTensor(const TPZVec<REAL> &pt, TPZVec<STATE> &kabs, TPZFMatrix<STATE> &tensorK);

REAL const pi = 4.*atan(1.);
bool fTriang = false;
bool IsStab = true;
bool IsContinuou = true;
bool Useh2 = true;
REAL Delta1 = 0.5;
REAL Delta2 = 0.5;
bool IsFullHdiv = false;
bool IsHomogeneo = false;

//erros
void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out);
void ErrorL22(TPZCompMesh *l2mesh, std::ostream &out);

void ComputeFluxError(TPZCompMesh *cmesh, std::ostream &out);
void ComputePressureError(TPZCompMesh *cmesh, std::ostream &out);

TPZFMatrix<STATE> * ComputeInverse(TPZCompMesh * mphysics);

void NEquationsCondensed(TPZCompMesh *cmesh, int64_t &neqglob,int64_t &neqcond, bool ismisto);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.hdiv"));
#endif

#include "pztransfer.h"
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    FileName = dirname + "/Projects/HDivEstabilizado/";
    FileName += "Hdivlog.cfg";
    InitializePZLOG(FileName);
#endif
    
    REAL Lx=0., Ly=0.;
    if(IsHomogeneo==true){
        Lx = 1.;
        Ly = 0.5;
    }
    
    std::string outputfile("out");
    
    //ofstream saidaerro( "../erros-hdiv-estab.txt",ios::app);
    ofstream saidaerro( "erros-hdiv-estab.txt");
    
    saidaerro<<"\n";
    saidaerro <<"INFORMACOES DO TESTE: (0)->False e (1)->True \n";
    saidaerro <<"Malha Triangular: "<<fTriang<<"\n";
    saidaerro <<"Problema eh homogeneo: "<< IsHomogeneo<<"\n";
    saidaerro <<"Metodo estabilizado: " <<IsStab<< ".    Usa parametro h2: "<< Useh2<<"\n";
    saidaerro <<"Espaco do fluxo eh full Hdiv: "<<IsFullHdiv<<"\n";
    saidaerro <<"Espaco da pressao eh continuo: "<<IsContinuou<<"\n\n";

    for(int p = 1; p<5; p++)
    {
        int pq = p;
        int pp;
        if(fTriang==true){
            pp = pq-1;
        }else{
            pp = p;
        }
        
        int ndiv;
        saidaerro<<"\n CALCULO DO ERRO, COM ORDEM POLINOMIAL pq = " << pq << " e pp = "<< pp <<endl;
        for (ndiv = 1; ndiv < 6; ndiv++)
        {
            
            //std::cout << "p order " << p << " number of divisions " << ndiv << std::endl;
            
            saidaerro<<"\n<<<<<< Numero de divisoes uniforme ndiv = " << ndiv <<" >>>>>>>>>>> "<<endl;
            
            TPZGeoMesh *gmesh = NULL;
            if(IsHomogeneo==true){
                gmesh = GMesh2(Lx, Ly,fTriang);
            }else{
                gmesh = GMesh3(fTriang);
            }
           // ofstream arg("gmesh1.txt");
            //gmesh->Print(arg);
            
            RefinamentoUnif(gmesh, ndiv);
            
            TPZCompMesh *cmesh1 = CMeshFlux2(pq,gmesh);
            TPZCompMesh *cmesh2 = CMeshPressure2(pp,gmesh);
//            
//            ofstream arg1("cmeshflux.txt");
//            cmesh1->Print(arg1);
//            
//            ofstream arg2("cmeshpressure.txt");
//            cmesh2->Print(arg2);

//            ofstream arg4("gmesh2.txt");
//            gmesh->Print(arg4);
            
            
            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(2);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            
            TPZCompMesh * mphysics = CMeshMixed2(meshvec,gmesh);
            
            mphysics->ExpandSolution();
            mphysics->CleanUpUnconnectedNodes();
            
            
//            ofstream arg4("gmeshMulti.txt");
//            mphysics->Print(arg4);
            
//            saidaerro << "Number of equations of flux " << cmesh1->NEquations() << std::endl;
//            saidaerro << "Number of equations of pressure " << cmesh2->NEquations() << std::endl;
//            saidaerro << "Number of equations TOTAL " << mphysics->NEquations() << "\n\n";
            
//            int64_t neq_flux, neq_pres, neqcond_flux, neqcond_pres;
//            NEquationsCondensed(cmesh1, neq_flux, neqcond_flux);
//            NEquationsCondensed(cmesh2, neq_pres, neqcond_pres);
//            saidaerro << "Number of equations total flux: " <<neq_flux<< "\n";
//            saidaerro << "Number of equations condensadas flux: " <<neqcond_flux<< "\n";
//            saidaerro << "Number of equations total pressao: " <<neq_pres<< "\n";
//            saidaerro << "Number of equations condensadas pressao: " << neqcond_pres<< "\n\n";
            
            int64_t neq_misto, neqcond_misto;
            NEquationsCondensed(mphysics, neq_misto, neqcond_misto,true);
            saidaerro << "Numero total de equacoes: " <<neq_misto<< "\n";
            saidaerro << "Numero de equacoes condensaveis: " <<neqcond_misto<< "\n";
            saidaerro << "Numero de equacoes final: " << neq_misto - neqcond_misto<< "\n\n";

            int numthreads = 0;
            std::cout << "Number of threads " << numthreads << std::endl;

            
            bool opti = true;
            TPZAnalysis  * an = new TPZAnalysis(mphysics,opti);
            ResolverSistema(*an, mphysics,numthreads, true);
            
            
//            ofstream arg5("cmeshmultiphysics.txt");
//            mphysics->Print(arg5);

            //Calculo do erro
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            TPZVec<REAL> erros;
    
            //saidaerro<<"\nErro da simulacao multifisica do fluxo (q)" <<endl;
            ErrorHDiv2(cmesh1, saidaerro);
            
            //saidaerro<<"\nErro da simulacao multifisica da pressao (p)" <<endl;
            ErrorL22(cmesh2, saidaerro);
            
            //Plot da solucao aproximada
            //string plotfile("Solution_mphysics.vtk");
            
//            char buf[256] ;
//            sprintf(buf,"ProblemaJuanGC_porder%d_h%d.vtk",p,ndiv);
//            PosProcessMultph(meshvec,  mphysics, *an, buf);

        }
   
    }
    
	return EXIT_SUCCESS;
}

TPZGeoMesh *GMesh2(REAL Lx, REAL Ly,bool triang_elements){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
	TPZVec <int64_t> TopolLine(2);
    TPZVec <int64_t> TopolPoint(1);
	
	//indice dos nos
	int64_t id = 0;
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
    
    if(triang_elements==true)
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,MatId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,MatId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC0,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC1,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC3,*gmesh);
    }
    else{
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 2;
        TopolQuad[3] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,MatId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC0,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC1,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC3,*gmesh);
    }
    
	gmesh->BuildConnectivity();
    
    #ifdef LOG4CXX
    	if(logger->isDebugEnabled())
    	{
            std::stringstream sout;
            sout<<"\n\n Malha Geometrica Inicial\n ";
            gmesh->Print(sout);
          LOGPZ_DEBUG(logger,sout.str())
    	}
    #endif
    
	return gmesh;
}


TPZGeoMesh *GMesh3(bool triang_elements){
    
    int Qnodes = 6;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
	TPZVec <int64_t> TopolLine(2);
//    TPZVec <int64_t> TopolPoint(1);
	
	//indice dos nos
	int64_t id = 0;
	REAL valx;
    REAL inix = -1.;
    REAL valy = -1.;
	for(int ix = 0; ix < Qnodes/2; ix++)
	{
		valx = inix + ix;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx);//coord X
		Node[id].SetCoord(1 ,valy);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
    inix = 1.;
    valy = 1.;
	for(int ix = 0; ix < Qnodes/2; ix++)
	{
		valx = inix - ix;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx);//coord X
		Node[id].SetCoord(1 ,valy);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
    
    if(triang_elements==true)
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 5;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,MatId,*gmesh);
        id++;
        
        TopolTriang[0] = 1;
        TopolTriang[1] = 4;
        TopolTriang[2] = 5;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,MatId,*gmesh);
        id++;
        
        TopolTriang[0] = 1;
        TopolTriang[1] = 2;
        TopolTriang[2] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,MatId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 3;
        TopolTriang[2] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,MatId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC0,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC1,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC3,*gmesh);
        id++;
        
        TopolLine[0] = 4;
        TopolLine[1] = 5;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC4,*gmesh);
        id++;
        
        TopolLine[0] = 5;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC5,*gmesh);
        id++;
    }
    else{
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 4;
        TopolQuad[3] = 5;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,MatId,*gmesh);
        id++;
        
        TopolQuad[0] = 1;
        TopolQuad[1] = 2;
        TopolQuad[2] = 3;
        TopolQuad[3] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,MatId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC0,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC1,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC3,*gmesh);
        id++;
        
        TopolLine[0] = 4;
        TopolLine[1] = 5;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC4,*gmesh);
        id++;
        
        TopolLine[0] = 5;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,BC5,*gmesh);
        id++;
    }
    
	gmesh->BuildConnectivity();
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"\n\n Malha Geometrica Inicial\n ";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
	return gmesh;
}


void RefinamentoUnif(TPZGeoMesh* gmesh, int nDiv)
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
    //	gmesh->BuildConnectivity();
}


TPZCompMesh *CMeshFlux2(int pOrder,TPZGeoMesh *gmesh)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(MatId,dim);
	TPZMaterial * mat(material);
	material->NStateVariables();
    
    //    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
	
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
	
    if(IsFullHdiv){
        cmesh->SetAllCreateFunctionsHDivFull();
    }
    else{
		cmesh->SetAllCreateFunctionsHDiv();
    }

    cmesh->InsertMaterialObject(mat);
    cmesh->SetDefaultOrder(pOrder);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    if(IsHomogeneo==true)
    {
//        TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcdirichlet, val1, val2);
//        TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcdirichlet, val1, val2);
//        TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcdirichlet, val1, val2);
//        TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcdirichlet, val1, val2);
//        
//        cmesh->InsertMaterialObject(BCond0);
//        cmesh->InsertMaterialObject(BCond1);
//        cmesh->InsertMaterialObject(BCond2);
//        cmesh->InsertMaterialObject(BCond3);
        
        ///Criar condicoes de contorno
        TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
        TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcneumann, val1, val2);
        TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcneumann, val1, val2);
        TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcneumann, val1, val2);
        TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcdirichlet, val1, val2);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
        bcmatNeumannAcima = new TPZDummyFunction<STATE>(NeumannAcima, 5);
        BCond2->SetForcingFunction(bcmatNeumannAcima);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAbaixo;
        bcmatNeumannAbaixo = new TPZDummyFunction<STATE>(NeumannAbaixo, 5);
        BCond0->SetForcingFunction(bcmatNeumannAbaixo);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatDirichletEsquerda;
        bcmatDirichletEsquerda = new TPZDummyFunction<STATE>(DirichletEsquerda, 5);
        BCond3->SetForcingFunction(bcmatDirichletEsquerda);
        
        cmesh->InsertMaterialObject(BCond0);
        cmesh->InsertMaterialObject(BCond1);
        cmesh->InsertMaterialObject(BCond2);
        cmesh->InsertMaterialObject(BCond3);
    }else{
        TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcdirichlet, val1, val2);
        TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcdirichlet, val1, val2);
        TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcdirichlet, val1, val2);
        TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcdirichlet, val1, val2);
        TPZMaterial * BCond4 = material->CreateBC(mat, BC4,bcdirichlet, val1, val2);
        TPZMaterial * BCond5 = material->CreateBC(mat, BC5,bcdirichlet, val1, val2);
        
        cmesh->InsertMaterialObject(BCond0);
        cmesh->InsertMaterialObject(BCond1);
        cmesh->InsertMaterialObject(BCond2);
        cmesh->InsertMaterialObject(BCond3);
        cmesh->InsertMaterialObject(BCond4);
        cmesh->InsertMaterialObject(BCond5);
    }
    
    cmesh->SetDimModel(dim);
	cmesh->AutoBuild();//Ajuste da estrutura de dados computacional
    
//#ifdef LOG4CXX
//	if(logdata->isDebugEnabled())
//	{
//        std::stringstream sout;
//        sout<<"\n\n Malha Computacional_1 Fluxo\n ";
//        cmesh->Print(sout);
//        LOGPZ_DEBUG(logdata,sout.str())
//	}
//#endif
	
	return cmesh;
}

TPZCompMesh *CMeshPressure2(int pOrder,TPZGeoMesh *gmesh)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(MatId,dim);
	material->NStateVariables();
    
    //    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDefaultOrder(pOrder);
    
    ///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    if(IsHomogeneo==true)
    {
        TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcdirichlet, val1, val2);
        TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcdirichlet, val1, val2);
        TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcdirichlet, val1, val2);
        TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcdirichlet, val1, val2);
        
        cmesh->InsertMaterialObject(BCond0);
        cmesh->InsertMaterialObject(BCond1);
        cmesh->InsertMaterialObject(BCond2);
        cmesh->InsertMaterialObject(BCond3);
    }else{
        TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcdirichlet, val1, val2);
        TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcdirichlet, val1, val2);
        TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcdirichlet, val1, val2);
        TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcdirichlet, val1, val2);
        TPZMaterial * BCond4 = material->CreateBC(mat, BC4,bcdirichlet, val1, val2);
        TPZMaterial * BCond5 = material->CreateBC(mat, BC5,bcdirichlet, val1, val2);
        
        cmesh->InsertMaterialObject(BCond0);
        cmesh->InsertMaterialObject(BCond1);
        cmesh->InsertMaterialObject(BCond2);
        cmesh->InsertMaterialObject(BCond3);
        cmesh->InsertMaterialObject(BCond4);
        cmesh->InsertMaterialObject(BCond5);
    }
    
    cmesh->SetDimModel(dim);
    
    if(IsContinuou==false)
    {
        cmesh->SetAllCreateFunctionsDiscontinuous();

        //Ajuste da estrutura de dados computacional
        cmesh->AutoBuild();

        int ncon = cmesh->NConnects();
        for(int i=0; i<ncon; i++)
        {
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            //newnod.SetPressure(true);
            newnod.SetLagrangeMultiplier(1);
        }

        int nel = cmesh->NElements();
        for(int i=0; i<nel; i++){
            TPZCompEl *cel = cmesh->ElementVec()[i];
            TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
            celdisc->SetConstC(1.);
            celdisc->SetCenterPoint(0, 0.);
            celdisc->SetCenterPoint(1, 0.);
            celdisc->SetCenterPoint(2, 0.);
            celdisc->SetTrueUseQsiEta();
            if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
            {
                if(fTriang==true) celdisc->SetTotalOrderShape();
                else celdisc->SetTensorialShape();
            }
            
        }


        #ifdef PZDEBUG
        int ncel = cmesh->NElements();
        for(int i =0; i<ncel; i++){
            TPZCompEl * compEl = cmesh->ElementVec()[i];
            if(!compEl) continue;
            TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
            if(facel)DebugStop();
            
        }
        #endif
    }
    else{
        cmesh->SetAllCreateFunctionsContinuous();
        
        //Ajuste da estrutura de dados computacional
        cmesh->AutoBuild();
    }
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_2 pressure\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str());
    //	}
    //#endif
	return cmesh;
}

TPZCompMesh *CMeshMixed2(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim =2;
    TPZMixedPoisson *material = new TPZMixedPoisson(MatId,dim);
    
    //incluindo os dados do problema
    REAL coefk = 1.;
    material->SetPermeability(coefk);
    REAL coefvisc = 1.;
    material->SetViscosity(coefvisc);
    
    //permeabilidade
    if(IsHomogeneo==true){
        TPZFMatrix<REAL> Ktensor(dim,dim,0.);
        TPZFMatrix<REAL> InvK(dim,dim,0.);
        Ktensor(0,0)=1.; Ktensor(1,1)=1.;
        InvK=Ktensor;
        material->SetPermeabilityTensor(Ktensor,InvK);
    }
    else{
        TPZAutoPointer<TPZFunction<STATE> > tensorK;
        tensorK = new TPZDummyFunction<STATE>(PermeabilityTensor, 5);
        material->SetPermeabilityFunction(tensorK);
    }
    
    if(IsStab==true){
        material->SetStabilizedMethod();
        material->SetStabilizationCoeficients(Delta1,Delta2);
	}
    if(IsStab==true && Useh2==true){
        material->SetHdois();
    }
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    if(IsHomogeneo==true){
        solexata = new TPZDummyFunction<STATE>(SolExataMista, 5);
        material->SetForcingFunctionExact(solexata);
    }else{
        solexata = new TPZDummyFunction<STATE>(SolFluxoHeter, 5);
        material->SetForcingFunctionExact(solexata);
    }
    
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    TPZDummyFunction<STATE> *dum;
    if(IsHomogeneo==true){
        dum = new TPZDummyFunction<STATE>(ForcingMista, 5);
    }else{
        dum = new TPZDummyFunction<STATE>(ForcingHeter, 5);
    }
    dum->SetPolynomialOrder(20);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    mphysics->SetDimModel(dim);
    //Criando condicoes de contorno
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    
    if(IsHomogeneo==true)
    {
        TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcneumann, val1, val2);
        TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcneumann, val1, val2);
        TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcneumann, val1, val2);
        TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcdirichlet, val1, val2);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
        bcmatNeumannAcima = new TPZDummyFunction<STATE>(NeumannAcima, 5);
        BCond2->SetForcingFunction(bcmatNeumannAcima);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAbaixo;
        bcmatNeumannAbaixo = new TPZDummyFunction<STATE>(NeumannAbaixo, 5);
        BCond0->SetForcingFunction(bcmatNeumannAbaixo);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatDirichletEsquerda;
        bcmatDirichletEsquerda = new TPZDummyFunction<STATE>(DirichletEsquerda, 5);
        BCond3->SetForcingFunction(bcmatDirichletEsquerda);
        
        mphysics->InsertMaterialObject(BCond0);
        mphysics->InsertMaterialObject(BCond1);
        mphysics->InsertMaterialObject(BCond2);
        mphysics->InsertMaterialObject(BCond3);
    }else{
        TPZMaterial * BCond0 = material->CreateBC(mat, BC0,bcneumann, val1,val2);
        TPZMaterial * BCond1 = material->CreateBC(mat, BC1,bcneumann, val1,val2);
        TPZMaterial * BCond2 = material->CreateBC(mat, BC2,bcdirichlet,val1,val2);
        TPZMaterial * BCond3 = material->CreateBC(mat, BC3,bcneumann, val1, val2);
        TPZMaterial * BCond4 = material->CreateBC(mat, BC4,bcneumann, val1, val2);
        TPZMaterial * BCond5 = material->CreateBC(mat, BC5,bcdirichlet,val1,val2);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAbaixoXMenorZero;
        bcmatNeumannAbaixoXMenorZero = new TPZDummyFunction<STATE>(NeumannAbaixoXMenorZero, 5);
        BCond0->SetForcingFunction(bcmatNeumannAbaixoXMenorZero);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAbaixoXMaiorZero;
        bcmatNeumannAbaixoXMaiorZero = new TPZDummyFunction<STATE>(NeumannAbaixoXMaiorZero, 5);
        BCond1->SetForcingFunction(bcmatNeumannAbaixoXMaiorZero);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatDirichletXIgualUm;
        bcmatDirichletXIgualUm = new TPZDummyFunction<STATE>(DirichletXIgualUm, 5);
        BCond2->SetForcingFunction(bcmatDirichletXIgualUm);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcimaXMaiorZero;
        bcmatNeumannAcimaXMaiorZero = new TPZDummyFunction<STATE>(NeumannAcimaXMaiorZero, 5);
        BCond3->SetForcingFunction(bcmatNeumannAcimaXMaiorZero);
        
        TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcimaXMenorZero;
        bcmatNeumannAcimaXMenorZero = new TPZDummyFunction<STATE>(NeumannAcimaXMenorZero, 5);
        BCond4->SetForcingFunction(bcmatNeumannAcimaXMenorZero);

        TPZAutoPointer<TPZFunction<STATE> > bcmatDirichletXIgualMenosUm;
        bcmatDirichletXIgualMenosUm = new TPZDummyFunction<STATE>(DirichletXIgualMenosUm, 5);
        BCond5->SetForcingFunction(bcmatDirichletXIgualMenosUm);

        mphysics->InsertMaterialObject(BCond0);
        mphysics->InsertMaterialObject(BCond1);
        mphysics->InsertMaterialObject(BCond2);
        mphysics->InsertMaterialObject(BCond3);
        mphysics->InsertMaterialObject(BCond4);
        mphysics->InsertMaterialObject(BCond5);
    }
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    //Fazendo auto build
    mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}

#include "TPZParFrontStructMatrix.h"
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numthreads, bool direct)
{
    
    if (direct)
    {
//        TPZSkylineStructMatrix full(fCmesh); //caso simetrico
//        full.SetNumThreads(numthreads);
//        an.SetStructuralMatrix(full);
//        TPZStepSolver<STATE> step;
//        step.SetDirect(ELDLt); //caso simetrico
//        //  step.SetDirect(ELU);
//        an.SetSolver(step);
//        an.Run();
        
//        TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(fCmesh);
        TPZFrontStructMatrix<TPZFrontSym<STATE> > strmat(fCmesh);
        strmat.SetDecomposeType(ELDLt);
        strmat.SetNumThreads(numthreads);
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt); //caso simetrico
        //	step.SetDirect(ELU);
        an.SetSolver(step);
        an.Run();
    }
    else{
        
        TPZSkylineStructMatrix full(fCmesh); //caso simetrico
        full.SetNumThreads(numthreads);
        an.SetStructuralMatrix(full);
        TPZStepSolver<STATE> step;
        TPZAutoPointer<TPZMatrix<STATE> > Inverse = ComputeInverse(fCmesh);
        TPZStepSolver<STATE> precond(Inverse);
        precond.SetMultiply();
//        step.SetGMRES(10, 40, precond, 1.0e-10, 0);
        step.SetCG(10, precond, 1.0e-10, 0);
        an.SetSolver(step);
        an.Run();
    }
	
	//Saida de Dados: solucao e  grafico no VT
//	ofstream file("Solutout");
//	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

TPZFMatrix<STATE> * ComputeInverse(TPZCompMesh * mphysics)
{
    int neq = mphysics->NEquations();
    TPZFMatrix<STATE> * PreInverse =  new TPZFMatrix<STATE> (neq,neq,0.0);
    TPZSkylineStructMatrix skyl(mphysics);
    std::set<int> matids; // to be computed
    matids.insert(MatId);
    matids.insert(BC0);
    matids.insert(BC1);
    matids.insert(BC2);
    matids.insert(BC3);
    skyl.SetMaterialIds(matids);
    TPZFMatrix<STATE> rhsfrac;
    TPZAutoPointer<TPZGuiInterface> gui = new TPZGuiInterface;
    TPZAutoPointer<TPZMatrix<STATE> > matfrac = skyl.CreateAssemble(rhsfrac, gui);
    TPZFMatrix<STATE> oldmat = *matfrac.operator->();

    
//#ifdef LOG4CXX
//    if(logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout.precision(20);
//        matfrac->Print("K = ",sout,EMathematicaInput);
//        LOGPZ_DEBUG(logger,sout.str());
//    }
//#endif
    
    matfrac->Inverse( * PreInverse,ELDLt);
    
    
//#ifdef LOG4CXX
//        if(logger->isDebugEnabled())
//        {
//            std::stringstream sout;
//            sout.precision(20);
////            TPZFMatrix<> Iden(neq,neq,0.);
////            oldmat.Multiply(*PreInverse, Iden);
//            PreInverse->Print("Kinv = ",sout,EMathematicaInput);
////            Iden.Print("Iden = ",sout,EMathematicaInput);
//            LOGPZ_DEBUG(logger,sout.str());
//        }
//#endif

    
    return PreInverse;
    
}

void PosProcessMultph(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(3), vecnames(4);
	vecnames[0]  = "Flux";
    vecnames[1]  = "GradFluxX";
    vecnames[2]  = "GradFluxY";
    vecnames[3]  = "ExactFlux";
    scalnames[0] = "Pressure";
    scalnames[1] = "DivFlux";
    scalnames[2] = "ExactPressure";

    
	const int dim = 2;
	int div =1;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
//	std::ofstream out("malha.txt");
//	an.Print("nothing",out);
    
}

void SolExataMista(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    double x = pt[0];
    double y = pt[1];
    
    solp[0] = (-2./pi)*cos(pi*x)*exp(y/2.);
    flux(0,0)= -2*sin(pi*x)*exp(y/2.);
    flux(1,0)= (1./pi)*cos(pi*x)*exp(y/2.);
    flux(2,0)= (1.- 4.*pi*pi)/(2.*pi)*(cos(pi*x)*exp(y/2.));
}

void SolExataPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    double x = pt[0];
    double y = pt[1];
    
    solp[0] = (-2./pi)*cos(pi*x)*exp(y/2.);
    flux(0,0)= 2*sin(pi*x)*exp(y/2.);
    flux(1,0)= -(1./pi)*cos(pi*x)*exp(y/2.);
    
}

void ForcingMista(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= (1.- 4.*pi*pi)/(2.*pi)*(cos(pi*x)*exp(y/2.));
}

void NeumannAbaixo(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    double x = pt[0];
    //double y = pt[1];
    disp[0] = -cos(pi*x)/pi;
}

void NeumannAcima(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    double x = pt[0];
    //double y = pt[1];
    disp[0] = 0.4087179842429694*cos(pi*x);
}

void DirichletEsquerda(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
//    TPZFMatrix<STATE> flux;
//    SolExataMista(pt, disp, flux);
    double y=pt[1];
    disp[0]=-(2./pi)*exp(y/2);
}

//meio heterogeneo
void SolFluxoHeter(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    REAL x = pt[0];
    REAL y = pt[1];
    
    if(x < 0.)
    {
        solp[0] = sin(y) + x*cos(y) + 2.*x*sin(y);
        flux(0,0)= -cos(y) - 2.*sin(y);//-KGradP[[1]]
        flux(1,0)= -cos(y) - 2.*x*cos(y) + x*sin(y);//-KGradP[[2]]
        flux(2,0)= x*cos(y) + sin(y) + 2.*x*sin(y);//div[-KgradP]
    }else
    {
        solp[0] = exp(x)*sin(y);
        flux(0,0) = -exp(x)*cos(y) - 2.*exp(x)*sin(y);//-KGradP[[1]]
        flux(1,0) = -2.*exp(x)*cos(y) - exp(x)*sin(y);//-KGradP[[2]]
        flux(2,0)= -2.*exp(x)*cos(y);//div[-KgradP]
    }
}

void SolPressaoHeter(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    REAL x = pt[0];
    REAL y = pt[1];
    
    if(x < 0.)
    {
        solp[0] = sin(y) + x*cos(y) + 2.*x*sin(y);
        flux(0,0)= cos(y) + 2.*sin(y);//dPdx
        flux(1,0)= cos(y) - x*sin(y) + 2.*x*cos(y);//dPdy
    }else
    {
        solp[0] = exp(x)*sin(y);
        flux(0,0) = exp(x)*sin(y);//dPdx
        flux(1,0) = exp(x)*cos(y);//dPdy
    }

}

void ForcingHeter(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
	REAL x = pt[0];
    REAL y = pt[1];
    
    if(x < 0.)
    {
        disp[0] = x*cos(y) + sin(y) + 2.*x*sin(y);
    }else
    {
        disp[0] = -2.*exp(x)*cos(y);
    }
}

void NeumannAbaixoXMenorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    if(x>0.) DebugStop();
    
    TPZManVector<STATE> p(1);
    TPZFNMatrix<10,STATE> fluxo(2,1);
    SolFluxoHeter(pt,p,fluxo);
    
    disp.Resize(1);
    disp[0] = 0.5403023058681398 + 1.922075596544176*x;//fluxo(0,0)*normal[0]+fluxo(1,0)*normal[1];
}


void NeumannAbaixoXMaiorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    if(x<0.) DebugStop();
    
    TPZManVector<STATE> p(1);
    TPZFNMatrix<10,STATE> fluxo(2,1);
    SolFluxoHeter(pt,p,fluxo);
    
    disp.Resize(1);
    disp[0] = 2.*exp(x)*cos(1.) - exp(x)*sin(1.);//fluxo(0,0)*normal[0]+fluxo(1,0)*normal[1];
    
}


void DirichletXIgualUm(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    REAL y = pt[1];
    if(x != 1.) DebugStop();
    
    TPZManVector<STATE> p(1);
    TPZFNMatrix<10,STATE> fluxo(2,1);
    SolFluxoHeter(pt,p,fluxo);
    
    disp.Resize(1);
    disp[0] = 2.718281828459045*sin(y);//p[0];
}

void NeumannAcimaXMaiorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    if(x<0.) DebugStop();
    
    TPZManVector<STATE> p(1);
    TPZFNMatrix<10,STATE> fluxo(2,1);
    SolFluxoHeter(pt,p,fluxo);
    
    disp.Resize(1);
    disp[0] = -2.*exp(x)*cos(1.) - exp(x)*sin(1.);//fluxo(0,0)*normal[0]+fluxo(1,0)*normal[1];
}

void NeumannAcimaXMenorZero(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    if(x>0.) DebugStop();
    
    TPZManVector<STATE> p(1);
    TPZFNMatrix<10,STATE> fluxo(2,1);
    SolFluxoHeter(pt,p,fluxo);
    
    disp.Resize(1);
    disp[0] = -0.23913362692838303*x - cos(1.);//fluxo(0,0)*normal[0]+fluxo(1,0)*normal[1];
}

void DirichletXIgualMenosUm(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    REAL y = pt[1];
    if(x != -1.) DebugStop();
    
    TPZManVector<STATE> p(1);
    TPZFNMatrix<10,STATE> fluxo(2,1);
    SolFluxoHeter(pt,p,fluxo);
    
    disp.Resize(1);
    disp[0] = sin(y) - cos(y) - 2.*sin(y);//p[0];
}


void PermeabilityTensor(const TPZVec<REAL> &pt, TPZVec<STATE> &kabs, TPZFMatrix<STATE> &tensorK){
    
     
    if(IsHomogeneo){
    tensorK.Resize(4,2);

        tensorK(0,0)=1.;
        tensorK(0,1)=0.;
        tensorK(1,0)=0.;
        tensorK(1,1)=1.;
        
        //kinv
        tensorK(2,0)=1.;
        tensorK(2,1)=0.;
        tensorK(3,0)=0.;
        tensorK(3,1)=1.;

    }
    
    else{
    REAL x = pt[0];
    tensorK.Resize(4,2);
    kabs.Resize(1, 0.);
    
    if(x<0){
        //K
        tensorK(0,0)=1.;
        tensorK(0,1)=0.;
        tensorK(1,0)=0.;
        tensorK(1,1)=1.;
        
        //kinv
        tensorK(2,0)=1.;
        tensorK(2,1)=0.;
        tensorK(3,0)=0.;
        tensorK(3,1)=1.;
    }
    else{
        //K
        tensorK(0,0)=2.;
        tensorK(0,1)=1.;
        tensorK(1,0)=1.;
        tensorK(1,1)=2.;
        
        //kinv
        tensorK(2,0)=2./3.;
        tensorK(2,1)=-1./3.;
        tensorK(3,0)=-1./3.;
        tensorK(3,1)=2./3.;
    }
    }
}

void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        if(IsHomogeneo==true){
            cel->EvaluateError(SolExataMista, elerror, 0);
        }else{
            cel->EvaluateError(SolFluxoHeter, elerror, 0);
        }
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
        
    }
   // out << "Errors associated with HDiv space\n";
    out << "L2 Norm for flux = "    << sqrt(globerrors[1]) << endl;
    out << "L2 Norm for divergence = "    << sqrt(globerrors[2])  <<endl;
    out << "Hdiv Norm for flux = "    << sqrt(globerrors[3])  <<endl;

}

void ErrorL22(TPZCompMesh *l2mesh, std::ostream &out)
{
    int64_t nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        if(IsHomogeneo==true){
            cel->EvaluateError(SolExataPressao, elerror, 0);
        }else{
            cel->EvaluateError(SolPressaoHeter, elerror, 0);
        }
        int nerr = elerror.size();
        globerrors.resize(nerr);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "\n";
    out << "Errors associated with L2 or H1 space\n";
    out << "H1 Norm = "    << sqrt(globerrors[0]) << endl;
    out << "L2 Norm = "    << sqrt(globerrors[1]) << endl;
    out << "Semi H1 Norm = "    << sqrt(globerrors[2]) << endl;
}

void ComputeFluxError(TPZCompMesh *cmesh, std::ostream &out){
    
    TPZManVector<REAL,2> errors(3,0.);
    
    int dimmesh = cmesh->Dimension();
    int nel = cmesh->NElements();
    int iel;
    for(iel=0; iel<nel; iel++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if(!cel) continue;
        
        int dimcel = cel->Dimension();
        if(dimcel != dimmesh) continue;
        
        TPZGeoEl *gel = cel->Reference();
        TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, 2);
        TPZManVector<int,3> prevorder(dimcel), maxorder(dimcel,intrule->GetMaxOrder());
        intrule->GetOrder(prevorder);
        intrule->SetOrder(maxorder);
        
        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
        TPZMaterialData data;
        sp->InitMaterialData(data);
        const int npoints = intrule->NPoints();
        TPZManVector<REAL,3> qsi(dimcel), xVec(3);
        
        for(int ip = 0; ip < npoints; ip++)
        {
            REAL weight;
            intrule->Point(ip,qsi,weight);
            sp->ComputeShape(qsi, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
            weight *= fabs(data.detjac);
            sp->ComputeSolution(qsi,data);
            
            TPZManVector<STATE,2> flux(2,0.);
            STATE divfluxo;
            flux[0]=data.sol[0][0];
            flux[1]=data.sol[0][1];
            divfluxo =  data.dsol[0](0,0)+data.dsol[0](1,1);
            
            TPZManVector<STATE> uExato(1);
            TPZFNMatrix<100,STATE> duExato(2,1);
            gel->X(qsi,xVec);
            SolExataMista(xVec, uExato, duExato);
        
            
            TPZManVector<STATE,2> diff(2,0.);
            diff[0] = flux[0]- duExato(0,0);
            diff[1] = flux[1]- duExato(1,0);
            
            //erro L2 do fluxo
            errors[0] += weight*(diff[0]*diff[0] + diff[1]*diff[1]);
            
            //erro L2 do divergente do fluxo
            STATE diffDiv = abs(divfluxo - duExato(2,0));
            errors[1] += weight*diffDiv*diffDiv;
            
            //erro Hdiv para o fluxo
//            errors[2] += errors[0] + errors[1];
        }
        intrule->SetOrder(prevorder);
    }
    
    //erro Hdiv para o fluxo
    errors[2] = errors[0] + errors[1];
    
    errors[0] = sqrt(errors[0]);
    errors[1] = sqrt(errors[1]);
    errors[2] = sqrt(errors[2]);
    
    out << "\n";
    out << "Erros associados ao fluxo na norma L2\n";
    out << "Norma L2  para fluxo = " << errors[0] << endl;
    out << "Norma L2 para divergente = " << errors[1] << endl;
    out << "Norma Hdiv para o fluxo = " << errors[2] << endl;
    
}///method

void ComputePressureError(TPZCompMesh *cmesh, std::ostream &out){
    
    TPZManVector<REAL,2> errors(3,0.);
    errors.Fill(0.);
    
    int dimmesh = cmesh->Dimension();
    int nel = cmesh->NElements();
    int iel;
    for(iel=0; iel<nel; iel++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if(!cel) continue;
        
        int dimcel = cel->Dimension();
        if(dimcel != dimmesh) continue;
        
        TPZGeoEl *gel = cel->Reference();
        TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, 2);
        TPZManVector<int,3> prevorder(dimcel), maxorder(dimcel,intrule->GetMaxOrder());
        intrule->GetOrder(prevorder);
        intrule->SetOrder(maxorder);
        
        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
        TPZMaterialData data;
        sp->InitMaterialData(data);
        const int npoints = intrule->NPoints();
        TPZManVector<REAL,3> qsi(dimcel), xVec(3);
        
        for(int ip = 0; ip < npoints; ip++)
        {
            REAL weight;
            intrule->Point(ip,qsi,weight);
            sp->ComputeShape(qsi, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphi, data.dphix);
            weight *= fabs(data.detjac);
            sp->ComputeSolution(qsi,data);
            
            TPZManVector<REAL,2> gradP(2,0.);
            REAL diffP;
            
            TPZManVector<STATE> uExato(1);
            TPZFNMatrix<100,STATE> duExato(2,1);
            gel->X(qsi,xVec);
            SolExataPressao(xVec, uExato, duExato);
            
            //erro L2 da pressao
            REAL solP = data.sol[0][0];
            diffP = solP - uExato[0];
            errors[1] += weight*(diffP*diffP);
            
            //erro semi H1 da pressao
            TPZFNMatrix<3,REAL> dsoldx;
            TPZFMatrix<REAL> dsoldaxes(2,1);
            dsoldaxes(0,0) = data.dsol[0][0];
            dsoldaxes(1,0) = data.dsol[0][1];
            
            TPZAxesTools<REAL>::Axes2XYZ(dsoldaxes, dsoldx, data.axes);
            TPZManVector<REAL,2> diffGrad(2,0.);
            diffGrad[0] = dsoldx(0,0)-duExato(0,0);
            diffGrad[1] = dsoldx(1,0)-duExato(1,0);
            errors[2] += weight*(diffGrad[0]*diffGrad[0] + diffGrad[1]*diffGrad[1]);
            
            //erro H1 para a pressao
            //errors[0] += errors[1] + errors[2];
        }
        intrule->SetOrder(prevorder);
    }
    
    //erro H1 para a pressao
    errors[0] = errors[1] + errors[2];
    
    errors[0] = sqrt(errors[0]);
    errors[1] = sqrt(errors[1]);
    errors[2] = sqrt(errors[2]);
    
    out << "\n";
    out << "Erros associados a pressao nas normas L2 e H1\n";
    out << "Norma H1 para a pressao = " << errors[0] << endl;
    out << "Norma L2 para a pressao = " << errors[1] << endl;
    out << "Norma semi-H1 para a pressao = " << errors[2] << endl;
    

}///method

void NEquationsCondensed(TPZCompMesh *cmesh, int64_t &neqglob,int64_t &neqcond, bool ismisto){
    
    int64_t ncon = cmesh->NConnects();
    neqglob = 0;
    neqcond = 0;
    for(int i = 0; i< ncon; i++){
        TPZConnect &co  = cmesh->ConnectVec()[i];
        //if(co.HasDependency()) continue;
        if(co.HasDependency()  || co.IsCondensed() || !co.NElConnected() || co.SequenceNumber() == -1) continue;
        
        int nelc = co.NElConnected();
        if (nelc==0) DebugStop();
        
        int dofsize = co.NShape()*co.NState();
        neqglob += dofsize;
        
        //equacoes condensaveis
        if (nelc == 1){
            neqcond += dofsize;
        }
    }
    
    if(ismisto && IsContinuou==false){
        int nel2D =0;
        for(int i=0; i<cmesh->NElements(); i++){
            TPZCompEl *cel = cmesh->ElementVec()[i];
            if(!cel) continue;
            if(cel->Reference()->Dimension() == cmesh->Dimension()){
                nel2D++;
            }
        }
        neqcond = neqcond - nel2D;
    }
}

