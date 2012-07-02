#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "pzanalysiserror.h"
#include "pzanalysis.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"
#include "pzmatrix.h"
#include "TPZCompElDisc.h"
#include "pzfstrmatrix.h"
#include "pzinterpolationspace.h"
#include "pzsubcmesh.h"
#include "pzlog.h"
#include "pzelctemp.h"
#include "pzelchdiv.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzfstrmatrix.h"
#include "pzgengrid.h"
#include "pzbndcond.h"
#include "pzmaterial.h"
#include "pzelmat.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pzlog.h"
#include <cmath>
#include "pzhdivpressure.h"

#include "TPZRefPattern.h"


#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("Steklov.main"));

#endif



void ValFunction(TPZVec<REAL> &loc, TPZFMatrix<STATE> &Val1, TPZVec<STATE> &Val2, int &BCType);
TPZGeoMesh * MalhaGeo(const int h);
TPZGeoMesh * MalhaGeoT(const int h);
/** Resolver o problema do tipo 
 * -Laplac(u) = 0
 * du/dn = lambda u em todo contorno
 */

using namespace std;


TPZCompMesh *CreateMesh2d(TPZGeoMesh &gmesh,int porder);

void PrettyPrint(TPZMaterialData &data, std::ostream &out);

void DrawCommand(std::ostream &out, TPZInterpolationSpace *cel, int shape, int resolution);

int main()
{
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		InitializePZLOG();
		std::stringstream sout;
		sout<< "Problema de Steklov"<<endl;
		LOGPZ_DEBUG(logger, sout.str());
	}
#endif
	int porder = 3;
    int h = 1;
			
    TPZGeoMesh *gmesh2 = MalhaGeo(h);//malha geometrica
			
			
    TPZCompMesh *cmesh = CreateMesh2d(*gmesh2,porder+1);//malha computacional
			
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    cmesh->LoadReferences();//mapeia para a malha geometrica lo
			
    TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
    
    TPZAnalysis analysis(cmesh);
    cmesh->SetName("Malha depois de Analysis-----");
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
	
    TPZCompEl *cel = cmesh->ElementVec()[0];
    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
    TPZMaterialData data;
    intel->InitMaterialData(data);
    TPZManVector<REAL> intpoint(2,0.);
    intel->ComputeRequiredData(data, intpoint);
    PrettyPrint(data, std::cout);
    int resolution = 5;
    int shapeindex = 4;
    DrawCommand(std::cout,intel, shapeindex, resolution);
    
    
	return 0;
}

TPZCompMesh *CreateMesh2d(TPZGeoMesh &gmesh,int porder){
	TPZCompEl::SetgOrder(porder);
	TPZCompMesh *comp = new TPZCompMesh(&gmesh);
	
	
	
	// Criar e inserir os materiais na malha
	TPZMatPoisson3d *mat = new TPZMatPoisson3d(1,2);
	TPZMaterial * automat(mat);
	comp->InsertMaterialObject(automat);
	
	
	// Condicoes de contorno
	TPZFMatrix<STATE> val1(1,1,1.),val2(1,1,0.);
	
	TPZMaterial *bnd = automat->CreateBC (automat,-1,2,val1,val2);//misto tbem
	comp->InsertMaterialObject(bnd);
	
	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-2,2,val1,val2);
	TPZBndCond *bndcond = dynamic_cast<TPZBndCond *> (bnd);
	bndcond->SetValFunction(ValFunction);
	comp->InsertMaterialObject(bnd);
	
	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-3,2,val1,val2);
	comp->InsertMaterialObject(bnd);
	
	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-4,2,val1,val2);
	comp->InsertMaterialObject(bnd);
	
   // comp->SetAllCreateFunctionsHDiv();
		comp->SetAllCreateFunctionsHDivPressure();
	//comp->SetAllCreateFunctionsContinuous();
	
	// Ajuste da estrutura de dados computacional
	comp->AutoBuild();
	
	
	comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
	comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
	comp->SetName("Malha Computacional Original");
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		comp->Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	
	
    return comp;
	
}


void ValFunction(TPZVec<REAL> &loc, TPZFMatrix<STATE> &Val1, TPZVec<STATE> &Val2, int &BCType)
{
	BCType = 2;
	Val1.Redim(1, 1);
	Val1(0,0) = 1.;
	Val2[0] = loc[0];
}


TPZGeoMesh * MalhaGeo(const int h){//malha quadrilatera
	
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	
	//Criar nos
	const int nnode = 4;//AQUI
	const int dim = 2;//AQUI
	
	REAL co[nnode][dim] = {{-1.,-1},{1.,-1},{1.,1.},{-1.,1.}};
	int indices[1][nnode]={{0,1,2,3}};//como serao enumerados os nos
	
	
	
	/*for(int i = 0; i < nnode; i++)
	{
		indices[0][i] = i;
	}
	*/
	
	int nod;
	TPZVec<REAL> coord(dim);
	for(nod=0; nod<nnode; nod++) {
		int nodind = gmesh->NodeVec().AllocateNewElement();
		
		for(int d = 0; d < dim; d++)
		{
			coord[d] = co[nod][d];
		}
		gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
	}
	//Criacao de elementos
	
	
	TPZVec<int> nodind(4);
	for(int i=0; i<4; i++){
		nodind[i] = indices[0][i];
	}
	
	int index;
	TPZGeoEl *elvec = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index); //AQUI
	
	gmesh->BuildConnectivity();
	
	//Cria as condicoes de contorno
	TPZGeoElBC gbc1(elvec,4,-1);// condicao de fronteira tipo -1: (x,y=0)
	TPZGeoElBC gbc2(elvec,5,-2);// condicao de fronteira tipo -2: (x=1,y)
	TPZGeoElBC gbc3(elvec,6,-3);// condicao de fronteira tipo -3: (x,y=1)
	TPZGeoElBC gbc4(elvec,7,-4);// condicao de fronteira tipo -4: (x=0,y)
	
	const std::string nameref;
	
 TPZAutoPointer<TPZRefPattern> ref;
	//gmesh->RefPatternList(EQuadrilateral);
	
	//	Refinamento uniforme
	for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for(int i = 0; i < n; i++){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}		
			
			
		}
		
		
	}
	
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		gmesh->Print(sout);
		LOGPZ_DEBUG(logger, sout.str());
	}
#endif 		 		 
	return gmesh;
} 

TPZGeoMesh * MalhaGeoT(const int h){//malha triangulo
	
	
	
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	
	//Criar ns
	const int nnode = 4;//AQUI
	const int nelem = 2;
	TPZGeoEl *elvec[nelem];	
	const int dim = 2;//AQUI
	
	REAL co[nnode][dim] = {{-1.,-1},{1.,-1},{1.,1.},{-1.,1.}};
	int indices[2][nnode];//como serao enumerados os nos
	
	//el 1
	indices[0][0] = 0;
	indices[0][1] = 1;
	indices[0][2] = 3;
	//el2
	indices[1][0] = 2;
	indices[1][1] = 3;
	indices[1][2] = 1;
	
	
	int nod;
	TPZVec<REAL> coord(dim);
	for(nod=0; nod<nnode; nod++) {
		int nodind = gmesh->NodeVec().AllocateNewElement();
		
		for(int d = 0; d < dim; d++)
		{
			coord[d] = co[nod][d];
		}
		gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
	}
	//Criacao de elementos
	
	
	TPZVec<int> nodind1(3);
	TPZVec<int> nodind2(3);
	for(int i=0; i<3; i++){
		nodind1[i] = indices[0][i];
		nodind2[i] = indices[1][i];
	}
	
	int index;
	elvec[0] = gmesh->CreateGeoElement(ETriangle,nodind1,1,index); //AQUI
	elvec[1] = gmesh->CreateGeoElement(ETriangle,nodind2,1,index); //AQUI
	
	
	gmesh->BuildConnectivity();
	
	//	como usar os padroes de refinamento do caju..nao posso rodar ainda pq tenho o Pz antigo..nao tenho este gRefDBase.
	//{
	//				TPZAutoPointer< TPZRefPattern > refExemplo = gRefDBase->FinfRefPattern("a_Quad_Rib_Side_4_5_6");
	//				if(refExemplo)
	//				{
	//						TPZGeoEl * elExemplo = gmesh->Elementvec()[10];
	//						elExemplo->SetRefPattern(refExemplo);
	//						elEmemplo->Divide(filhos);
	//				}
	//		}
	
	//Cria as condicoes de contorno
	TPZGeoElBC gbc1(elvec[0],3,-1);// condicao de fronteira tipo -1: 
	TPZGeoElBC gbc2(elvec[0],5,-2);// condicao de fronteira tipo -2: 
	
	TPZGeoElBC gbc3(elvec[1],3,-3);// condicao de fronteira tipo -3: 
	TPZGeoElBC gbc4(elvec[1],5,-4);// condicao de fronteira tipo -4: 
	
	const std::string nameref;
	
	TPZAutoPointer<TPZRefPattern> ref;
	
	
	//	Refinamento uniforme
	/*for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for(int i = 0; i < n; i++){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}		
			}
		
	}*/
	for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for(int i = 0; i < n; i++){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}		
			
			
		}
		
		
	}
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		gmesh->Print(sout);
		LOGPZ_DEBUG(logger, sout.str());
	}
#endif 		 		 
	return gmesh;
}

void PrettyPrint(TPZMaterialData &data, std::ostream &out)
{
    int numfunction = data.fVecShapeIndex.size();
    for (int i=0; i<numfunction; i++) {
        int vecindex = data.fVecShapeIndex[i].first;
        int shapeindex = data.fVecShapeIndex[i].second;
        REAL shapeval = data.phi(shapeindex,0);
        TPZManVector<REAL,3> vec(3,0.);
        for (int iv=0; iv<3; iv++) {
            vec[iv] = data.fNormalVec(iv,vecindex);
        }
        out << "function " << i << " val = " << shapeval << " * (" << vec << ")\n";
    }
}

void DrawCommand(std::ostream &out, TPZInterpolationSpace *cel, int shape, int resolution)
{
    
    out << "master = Graphics[Line[{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {-1, -1}}]];\n";
    
    std::stringstream allgraphstr;
    allgraphstr << "AllGraph" << shape;
    std::string allgraph = allgraphstr.str();
    out << allgraph << " = {master};\n";
    
    TPZMaterialData data;
    cel->InitMaterialData(data);
    
    for (int i=0; i< resolution+2; i++) {
        for (int j=0; j< resolution+2; j++) {
            TPZManVector<REAL,2> point(2,0.);
            point[0] = -1. + i*2./(resolution+1.);
            point[1] = -1. + j*2./(resolution+1.);
            cel->ComputeRequiredData(data, point);
            int vecindex = data.fVecShapeIndex[shape].first;
            int shapeindex = data.fVecShapeIndex[shape].second;
            TPZManVector<REAL, 2> vector1(2,0.), vector2(2,0.);
            for (int v=0; v<2; v++) {
                vector2[v] = point[v]+data.phi(shapeindex,0)*data.fNormalVec(v,vecindex);
                if (abs(vector2[v]) < 1.e-4) {
                    vector2[v] = 0.;
                }
            }
            std::stringstream sout;
            sout << "Graphics[{Arrowheads[0.05],Thickness[0.01],Arrow[{{" << point << "},{" << vector2 << "}}]}]";
            out << "AppendTo[" << allgraph << "," << sout.str() << "];\n";
        }
    }
    
}
