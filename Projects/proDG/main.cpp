/**
 * @file
 * @brief Implements the build of a computational mesh as tutorial example of the interpolation NeoPZ module
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
#include "pzbndcond.h"
#include <TPZMatLaplacian.h>
#include "tpzdifureac.h"
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

#include <iostream>
#include <fstream>
using namespace std;

// cria a malha geométrica, dados os comprimentos em x e em y, o argumento bool é para escolher a malha com triangulos ou quadrilateros
TPZGeoMesh * GetMesh(REAL Lx,REAL Ly, bool triang_elements);

//nDiv indica o número de divisões feitas em cada direção
void UniformRefine(TPZGeoMesh* gmesh, int nDiv);

//cria a malha computational
TPZCompMesh *MalhaComp(TPZGeoMesh * gmesh, int pOrder);
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, bool symmetric_matrix);

//pos processamento da solucao
void SaidaSolucao(TPZAnalysis &an, std::string plotfile);
void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC0(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBC3(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &deriv);
//os ids das bc se colocam com iinteiros negativos
int matId = 1;
int bc0 = -1;	
int bc1 = -2; 
int bc2 = -3;
int bc3 = -4;

int dirichlet = 0;
int neumann = 1;

REAL erroL2;

int main() {
    InitializePZLOG();
    
     std::ofstream erro("Calculotaxa.txt");//este archivo guarda os erros e taxas
     TPZVec<REAL> calcErro;
     
     for(int pOrder=1; pOrder<4; pOrder++){
     
        erro<<"Ordem polinomial"<<pOrder<<endl; //escreve no arquivo onde guarda os erros
	
	for(int ndiv =0; ndiv<5	; ndiv++){
	  TPZGeoMesh *gmesh = GetMesh(2,2,true);  //lx, ly, triangles, or not(quadrilaterals)
	  UniformRefine(gmesh, ndiv);//refinamento uniforme
     
	    erro<<"\n Número de divisoes: "<<ndiv<<endl;
	    std::ofstream out2("geomesh1.txt"); //malha geométrica refinada
	    gmesh->Print(out2);
     
	    std::ofstream Dummyfile("GeometricMesh.vtk");
	    TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
     
// ----- Malha computacional -----
     
	    TPZCompMesh *cmesh = MalhaComp(gmesh,pOrder);
//   	    std::ofstream out3("compmesh0.txt"); //malha computacional, vai depois do refinamento
//          cmesh->Print(out3);
     
	    cmesh->ApproxSpace().CreateInterfaces(*cmesh);//cria as interfaces na malha comput.
//	    std::ofstream out5("compmesh1.txt");
//	    cmesh->Print(out5);
     
	 /*** Stiff Matrix ***/

       int64_t neq = cmesh->NEquations();
       TPZFMatrix<STATE> rhs(neq,1);//right hand side
       TPZFStructMatrix full(cmesh);
       TPZMatrix<STATE> *Stiff = full.CreateAssemble(rhs, 0);
       ofstream ArgStiff("Stiff.txt");
       Stiff->Print("Stiff=",ArgStiff,EMathematicaInput); //imprime matriz de rigidez em formato legível por Mathematica
      
	    TPZAnalysis an(cmesh);
	    
// 	    an.LoadSolution();
	    ResolverSistema(an, cmesh, false); //o booleano indica q o sistema a resolver é simétrico ou nao, configurar como true no caso SIPG
     
	    // Arquivo de saida para plotar a solução
	    string plotfile("Solution_aproximada.vtk");
	    SaidaSolucao(an,plotfile);
      
	    //calculo do erro
      
	    an.SetExact(*SolExata);//pasar a solução exata
	    an.PostProcess(calcErro,erro);//clcula erro em norma L2, H1 e seminorma
	    
	    //arquivos com graficos das solucoes
	    //TPZVec<std::string> scalnames(1),vecnames(1);
	    //scalnames[0] ="Solution";
	    //vecnames[0]="Derivative";
	      //    char buf[256] ;//nome do arquivo
              //   sprintf(buf,"GrafSol_pOrder%d_ndiv_%d.vtk",pOrder,ndiv);//arquivo para cada aproximação ser vista em paraview, ordem polinomial e # de divisões
              //   an.DefineGraphMesh(2,scalnames,vecnames,buf);
	      //   an.PostProcess(2); //o parámetro indica a resolução para visualizar
	    
	    
	}
    } 
	return 0;
	
}

TPZGeoMesh *GetMesh (REAL Lx,REAL Ly, bool triang_elements){
  
	int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	  
	TPZVec <int64_t> TopolQuad(4);
       TPZVec <int64_t> TopolTriang(3);
	TPZVec <int64_t> TopolLine(2);
	
	//indice dos nos
	int64_t id = 0;
	REAL valx, valy;
	REAL x00 = -1., y00 = -1.; //Valores mínimos en cada coordenada del domínio
    
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		//valx = -1. + xi*Lx; //quadrilatero [-1,1] x [-1,1]
		//valy = -1.; 
		valx = x00 + xi*Lx;
		valy = y00;		
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0,valx);//coord X
		Node[id].SetCoord(1,valy);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		//valx = Lx - xi*Lx - 1;
		//valy = -1. + Ly;
		valx = Lx - xi*Lx +x00;
		valy = y00+Ly;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0,valx);//coord X
		Node[id].SetCoord(1,valy);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
    
    if(triang_elements==true)
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 3;
        TopolTriang[2] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 0;
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
    else{
        
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 2;
        TopolQuad[3] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
        id++;
        
	//elementos de contorno
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
	return gmesh;
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

TPZCompMesh *MalhaComp(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	
	TPZdifureac *material = new TPZdifureac(matId,dim);
	TPZMaterial * mat(material);
	
	//O problema que resolve o material é "-K Laplaciano(u) + Alpha u = f", com alpha constante e K um múltiplo da matriz identidade
	
	REAL parK = 1.;//0.;
	REAL parAlpha = 0.;//1.;//-1;
	REAL functionf = -1.;//0;//-4.;
	
	REAL sigma = 1.0;//18.0; //6.0 for pOrder=1, 18 for  pO = 2 and 36 for pOrder=36
	REAL beta = 1.0;	//penalty term sigma/\e\^beta, multiplica [u][v]
	
	material->SetPenaltyConstants(sigma,beta);//configurar termos de penalidade
	material->SetParameters(parK,functionf,parAlpha); //configurar parametros, si f é não constante aquí se passa um valor qualquer e depois se insere a função certa
	//material->NotSetSymetricterm(); //termo de simetría 0. Sem termo de simetria deve ter penalidade ou a matriz é singular
	//material->SetSymmetric(); //termo de simetría -1
	material->SetNonSymmetric();//termo de simetría 1
    	material->NStateVariables();
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	
	cmesh->SetDimModel(dim);
    	cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->InsertMaterialObject(mat); //mat2);
    
	//funcao do lado direito da equacao do problema
     TPZAutoPointer<TPZFunction<STATE> > forcef;
     TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>((void (*)(const TPZVec<REAL> &, TPZVec<STATE> &))Forcing, 5);
     dum->SetPolynomialOrder(20); //?
     forcef = dum;
     material->SetForcingFunction(forcef); //no caso de que o termo fonte da equação esteja definido por uma função é preciso fazer isto
     
    TPZAutoPointer<TPZFunction<STATE> > FExact;
    FExact = new TPZDummyFunction<STATE>(SolExata, 5);
    material->SetForcingFunctionExact(FExact);    
     
	///Inserir condicao de contorno
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.); // ?
    REAL uDbc0 =0.;
    val2(0,0) = uDbc0;
    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    
    REAL uDbc2 =0.;
    val2(0,0) = uDbc2;
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    
    REAL uDbc1 =0.;
    val2(0,0) = uDbc1;
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    
    REAL uDbc3 =0.;
    val2(0,0) = uDbc3;
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    
    //insere condições de contorno não constantes
     TPZAutoPointer<TPZFunction<STATE> > fCC0;
       fCC0 = new TPZDummyFunction<STATE>((void (*)(const TPZVec<REAL> &, TPZVec<STATE> &))ForcingBC0, 5);
     TPZAutoPointer<TPZFunction<STATE> > fCC1;
       fCC1 = new TPZDummyFunction<STATE>((void (*)(const TPZVec<REAL> &, TPZVec<STATE> &))ForcingBC1, 5);
     TPZAutoPointer<TPZFunction<STATE> > fCC2;
       fCC2 = new TPZDummyFunction<STATE>((void (*)(const TPZVec<REAL> &, TPZVec<STATE> &))ForcingBC2, 5);
     TPZAutoPointer<TPZFunction<STATE> > fCC3;
       fCC3 = new TPZDummyFunction<STATE>((void (*)(const TPZVec<REAL> &, TPZVec<STATE> &))ForcingBC3, 5);
//          
        BCond0->SetForcingFunction(fCC0);
  	BCond1->SetForcingFunction(fCC1);
  	BCond2->SetForcingFunction(fCC2);
  	BCond3->SetForcingFunction(fCC3);
    	
	cmesh->InsertMaterialObject(BCond0);
	cmesh->InsertMaterialObject(BCond1);
	cmesh->InsertMaterialObject(BCond2);
	cmesh->InsertMaterialObject(BCond3);
    
    	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
      	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
	
	
	return cmesh;
}

void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, bool symmetric_matrix)
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
	
	//Saida de Dados: solucao txt
	ofstream file("Solution.out");
	an.Solution().Print("solution", file);
	
}

void SaidaSolucao(TPZAnalysis &an, std::string plotfile){
    
	TPZManVector<std::string,10> scalnames(2), vecnames(2);
	scalnames[0] = "Solution";
	scalnames[1] = "ExactSolution";
	vecnames[0]= "Derivative";
	vecnames[1]= "ExactDerivative";
	
    
	const int dim = an.Mesh()->Dimension();
	int div = 3;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}
void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE>&disp){ //função do lado direito
    
    double x = pt[0];
    double y = pt[1];
    //disp[0]= 2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);
 disp[0]= (1-4.*y*y)*exp(-x-y*y);   
  //disp[0]=1 + 2.*x + 3.*y +4.*x*y+ 5.*x*x+6*y*y;
  //  disp[0] = 1.+2.*x + 3.*x*x + 4.*x*y + 5*x*x*y + 6.*x*x*y*y + 7.*y + 8.*y*y + 9.*x*y*y; //esta es la misma solucion, teste k=0, alpha = 1
  //disp[0] = -21.+ (-16.*x) -9.*x*x + 4.*x*y + 5*x*x*y + 6.*x*x*y*y -3.*y - 4.*y*y + 9.*x*y*y;  
}
void ForcingBC0(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    double x = pt[0];
    double y = pt[1];
    //disp[0]= 0.5*x*(x-1.);//exp(-x*x-y*y);
    //disp[0] = -(7. + 4.*x + 5.*x*x + 12*x*x*y + 16.*y + 18.*x*y);
    disp[0]=exp(-x-y*y);
   // disp[0]=1 + 2.*x + 3.*y +4.*x*y+ 5.*x*x+6*y*y;
}
void ForcingBC2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];
    //disp[0]= -Pi*cos(Pi*y)*sin(Pi*x);
    //disp[0]= 7.+ 4.*x + 5.*x*x + 12*x*x*y+ 16.*y + 18.*x*y;
    disp[0]=exp(-x-y*y);
  //  disp[0]=1 + 2.*x + 3.*y +4.*x*y+ 5.*x*x+6*y*y;
}

void ForcingBC1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];
    //disp[0]= -Pi*cos(Pi*y)*sin(Pi*x);
   // disp[0] = 2. + 6.*x + 4.*y + 10.*x*y + 12.*x*y*y + 9.*y*y ;
   // disp[0] = 6.+ 16.*y + 23.*y*y;
   disp[0]=exp(-x-y*y);
  //disp[0]=1 + 2.*x + 3.*y +4.*x*y+ 5.*x*x+6*y*y;
}

void ForcingBC3(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    double x = pt[0];
    double y = pt[1];
  //  disp[0] = -(2. + 6.*x + 4.*y + 10.*x*y + 12.*x*y*y + 9.*y*y) ;
    //disp[0]= -Pi*cos(Pi*y)*sin(Pi*x);

    disp[0]=exp(-x-y*y);

  //disp[0]=1 + 2.*x + 3.*y +4.*x*y+ 5.*x*x+6*y*y;
}
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &deriv) {
	double x = pt[0];
	double y = pt[1];
	
	//teste 1
 	p[0]= exp(-x-y*y);
	/*p[0]=1 + 2.*x + 3.*y +4.*x*y+ 5.*x*x+6*y*y;
 	deriv(0,0)=2. + 4*y +10.*x;
 	deriv(1,0)=3. + 4.*x + 12.*y;
 	deriv(2,0)=0.;
	*/
	deriv(0,0)=-p[0];
 	deriv(1,0)=-2.*y*p[0];
 	deriv(2,0)=0.;
	//teste 2
// 	p[0]= 1.+2.*x + 3.*x*x + 4.*x*y + 5*x*x*y + 6.*x*x*y*y + 7.*y + 8.*y*y + 9.*x*y*y;//0.5*x*(x-1.);
// 	deriv(0,0)=2. + 6.*x + 4.*y + 10.*x*y + 12.*x*y*y + 9.*y*y ;//x-0.5; //-2.*x*p[0];
// 	deriv(1,0)=4.*x + 5.*x*x + 12*x*x*y + 7 + 16.*y + 18.*x*y; //0.; //-2.*y*p[0];
// 	deriv(2,0)=0.; // deriv c respeito a z
  
  
}

