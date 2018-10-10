//
//  ToolsOptimal.cpp
//  PZ
//
//  Created by Denise de Siqueira on 6/18/15.
//
//

#include "ToolsOptimal.h"
#include "pzpoisson3d.h"
#include"pzgeoelbc.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgengrid.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "PZMatPoissonControl.h"


#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("HpAdaptivity.main"));

#endif

const int matId = 1;
const int dirichlet = 0;

const int bc0 = -1;
const int bc1 = -2;
const int bc2 = -3;
const int bc3 = -4;


TPZGeoMesh *OptimalGeoMesh(bool ftriang, REAL Lx, REAL Ly){
    
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
    
    if(ftriang==true)
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    else{
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 2;
        TopolQuad[3] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    
	gmesh->BuildConnectivity();
    
//#ifdef LOG4CXX
//    if(Logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout<<"\n\n Malha Geometrica Inicial\n ";
//        gmesh->Print(sout);
//        LOGPZ_DEBUG(Logger,sout.str())
//    }
//#endif
    
	return gmesh;
}

TPZCompMesh *OptimalCompMesh(TPZGeoMesh *gmesh,int porder){

    
    
    REAL dim = 2;
    TPZMatPoissonControl *material = new TPZMatPoissonControl( matId,dim);
//	TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    material->NStateVariables();
	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
	
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	
    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
	
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(StateVar, 5);
    material->SetForcingFunctionExact(solexata);
	
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > force;
    TPZDummyFunction<STATE> *dum;
    dum = new TPZDummyFunction<STATE>(OptForcing, 5);
    dum->SetPolynomialOrder(20);
    force = dum;
    material->SetForcingFunction(force);
	
    
    cmesh->SetDefaultOrder(porder);
    cmesh->SetDimModel(dim);
	
	
    cmesh->SetAllCreateFunctionsContinuous();
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    
    return cmesh;

    
    cmesh->SetName("Optimal CompMesh");
#ifdef LOG4CXX
    {
        std::stringstream sout;
        cmesh->Print(sout);
       // LOGPZ_DEBUG(Logger,sout.str())
    }
#endif

    
    
    return cmesh;
    
}

void StateVar(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    u.Resize(1, 0.);
    du.Resize(3, 1.);
    du(0,0)=du(1,0)=du(2,0)=0.;
    
    const REAL sol = 10.*x*y*(1-x)*(1-y);
    u[0] = sol;
    
    
}

void OptForcing(const TPZVec<REAL> &pt, TPZVec<STATE> &res){
	
    res[0]=0;
    res[1]=0;
	double x = pt[0];
    double y = pt[1];
    res[0]= 0.;
    res[1]=10.*x*y*(1-x)*(1-y);
}

void SolveLUOpt ( TPZAnalysis &an ){
    TPZCompMesh *malha = an.Mesh();
    time_t tempoinicial;
    time_t tempofinal;
    time(&tempoinicial);
    //TPZFrontStructMatrix<TPZFrontNonSym> mat ( malha );// não funciona com método iterativo
    TPZSkylineStructMatrix mat(malha);
    //	TPZFStructMatrix mat( malha );
    //	TPZSpStructMatrix mat( malha );
    TPZStepSolver<STATE> solv;
    
    //solv.SetDirect ( ELU );
    //solv.SetDirect(ECholesky);
    solv.SetDirect(ELDLt);
    
    std::cout << "ELU " << std::endl;
    an.SetSolver ( solv );
    an.SetStructuralMatrix ( mat );
    std::cout << std::endl;
    an.Solution().Redim ( 0,0 );
    std::cout << "Assemble " << std::endl;
    time(&tempofinal);
    std::cout << "Antes do assemble " << tempofinal - tempoinicial << std::endl;
    
    time(&tempoinicial);
    an.Assemble();
    time(&tempofinal);
    std::cout << "  Tempo do assemble " << tempofinal - tempoinicial << std::endl;
    
    time(&tempoinicial);
    an.Solve();
    time(&tempofinal);
    std::cout << "  Tempo do assemble " << tempofinal - tempoinicial << std::endl;
    
    std::cout << std::endl;
    std::cout << "  No equacoes = " << malha->NEquations() << std::endl;
    // cout << "Banda = " << malha->BandWidth() << endl;
}
