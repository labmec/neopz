#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"

#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzinterpolationspace.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzmaterial.h"
#include "pzmat2dlin.h"
#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"

#include "pzanalysis.h"

#include "pzlog.h"

#include "pzelmat.h"


#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("HDivDemo"));

#endif



void ValFunction(TPZVec<REAL> &loc, TPZFMatrix<STATE> &Val1, TPZVec<STATE> &Val2, int &BCType);

TPZGeoMesh * ReadGmsh(std::string filename);
TPZGeoMesh * MalhaGeo(const int h);
TPZGeoMesh * MalhaGeoT(const int h);
/** Resolver o problema do tipo 
 * -Laplac(u) = 0
 * du/dn = lambda u em todo contorno
 */

using namespace std;


TPZCompMesh *CreateHDivMesh2d(TPZGeoMesh &gmesh,int porder);
TPZCompMesh *CreatePressureMesh2d(TPZGeoMesh &gmesh,int porder);
TPZCompMesh *CreateMultiPhysicsMesh(TPZGeoMesh *gmesh, TPZCompMesh *cmeshHDiv, TPZCompMesh *cmeshPressure);

void PrettyPrint(TPZMaterialData &data, std::ostream &out);

void DrawCommand(std::ostream &out, TPZInterpolationSpace *cel, int shape, int resolution);

int main()
{
	
#ifdef LOG4CXX
    InitializePZLOG();
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout<< "Problema de Steklov"<<endl;
		LOGPZ_DEBUG(logger, sout.str());
	}
#endif
	int porder = 1;
			
//    TPZGeoMesh *gmesh2 = MalhaGeo(h);//malha geometrica
    
    // @omar:: Dear Phil, I have prepared a gmsh script which give examples of the following issues:
    // Entities tagging in gmsh
    // Refinement
    // Macro defintions
    // please see "fractured_reservoir.geo" file
    std::string gmsh_geometry = PZSOURCEDIR;
    gmsh_geometry += "/Projects/Hdiv_demo/NaturalFractures/fractured_reservoir.msh";
    TPZGeoMesh *gmesh2 = ReadGmsh(gmsh_geometry);
    DebugStop();// @omar:: Please check the material id of the fractures, boundaries and the reservoir.
    
    TPZCompMesh *cmeshHDiv = CreateHDivMesh2d(*gmesh2,porder);//malha computacional
    TPZCompMesh *cmeshPress = CreatePressureMesh2d(*gmesh2,porder);
    
    TPZCompMesh *multiPhysics = CreateMultiPhysicsMesh(gmesh2, cmeshHDiv, cmeshPress);
			
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        multiPhysics->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    multiPhysics->LoadReferences();//mapeia para a malha geometrica lo
			
    TPZAdmChunkVector<TPZCompEl *> elvec = multiPhysics->ElementVec();
    
    TPZAnalysis analysis(multiPhysics);
    multiPhysics->SetName("Malha depois de Analysis-----");
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        multiPhysics->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
	
    {
        
        TPZCompEl *cel = cmeshHDiv->ElementVec()[0];
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        TPZMaterialData data;
        intel->InitMaterialData(data);
        TPZManVector<REAL> intpoint(2,0.);
        intel->ComputeRequiredData(data, intpoint);
        PrettyPrint(data, std::cout);
        int resolution = 5;
        int shapeindex = 10;
        DrawCommand(std::cout,intel, shapeindex, resolution);
    }
    
    {
        TPZCompEl *cel = multiPhysics->Element(0);
        TPZElementMatrix ek,ef;
        cel->CalcStiff(ek, ef);
        ek.fMat.Print("EK=",std::cout,EMathematicaInput);
    }
    
	return 0;
}

TPZCompMesh *CreateHDivMesh2d(TPZGeoMesh &gmesh,int porder){
	TPZCompEl::SetgOrder(porder);
	TPZCompMesh *comp = new TPZCompMesh(&gmesh);
	
    comp->SetDimModel(2);
	
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
		comp->SetAllCreateFunctionsHDiv();
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

TPZCompMesh *CreatePressureMesh2d(TPZGeoMesh &gmesh,int porder)
{
    TPZMat2dLin *mat = new TPZMat2dLin(1);
    TPZFNMatrix<1,STATE> xk(1,1,1.),xc(1,1,0.),xf(1,1,0.);
    mat->SetMaterial(xk, xc, xf);
    TPZCompMesh *cmesh = new TPZCompMesh(&gmesh);
    cmesh->SetDefaultOrder(porder);
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDimModel(2);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    return cmesh;
}

TPZCompMesh *CreateMultiPhysicsMesh(TPZGeoMesh *gmesh, TPZCompMesh *cmeshHDiv, TPZCompMesh *cmeshPressure)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    TPZMixedPoisson *mixed = new TPZMixedPoisson(1,2);
    cmesh->InsertMaterialObject(mixed);

    // Condicoes de contorno
    TPZFMatrix<STATE> val1(1,1,1.),val2(1,1,0.);
    
    TPZMaterial *bnd = mixed->CreateBC (mixed,-1,2,val1,val2);//misto tbem
    cmesh->InsertMaterialObject(bnd);
    
    // Mixed
    val1(0,0) = 1.;
    val2(0,0)=0.;
    bnd = mixed->CreateBC (mixed,-2,2,val1,val2);
    TPZBndCond *bndcond = dynamic_cast<TPZBndCond *> (bnd);
    bndcond->SetValFunction(ValFunction);
    cmesh->InsertMaterialObject(bnd);
    
    // Mixed
    val1(0,0) = 1.;
    val2(0,0)=0.;
    bnd = mixed->CreateBC (mixed,-3,2,val1,val2);
    cmesh->InsertMaterialObject(bnd);
    
    // Mixed
    val1(0,0) = 1.;
    val2(0,0)=0.;
    bnd = mixed->CreateBC (mixed,-4,2,val1,val2);
    cmesh->InsertMaterialObject(bnd);
    
    // comp->SetAllCreateFunctionsHDiv();
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    //comp->SetAllCreateFunctionsContinuous();
    
    // Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    TPZManVector<TPZCompMesh *,2> meshvector(2);
    meshvector[0] = cmeshHDiv;
    meshvector[1] = cmeshPressure;

    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);
    
    return cmesh;

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
	
	
	TPZVec<long> nodind(4);
	for(int i=0; i<4; i++){
		nodind[i] = indices[0][i];
	}
	
	long index;
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
	
	
	TPZVec<long> nodind1(3);
	TPZVec<long> nodind2(3);
	for(int i=0; i<3; i++){
		nodind1[i] = indices[0][i];
		nodind2[i] = indices[1][i];
	}
	
	long index;
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

TPZGeoMesh * ReadGmsh(std::string filename)
{
    TPZGmshReader readgmsh;
    TPZGeoMesh * geometry = readgmsh.GeometricGmshMesh(filename);
    ofstream textfile("geometry.txt");
    geometry->Print(textfile);
    std::ofstream vtkfile("geometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geometry, vtkfile, true);
    
    return geometry;
}
