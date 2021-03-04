/**
 * @file
 * @brief Simple 2D example simulating the electric potential between the plates of a parallel-plate capacitor given that its plates are located at y=-x (with V=0) and y=8-x ( with V=1).
 * @author Francisco Orlandini & Nathan Shauer
 * @since 2014
 */

///Bibliotecas utilizadas por este programa
#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "TPZMatExSimples2D.h"
#include "TPZParFrontStructMatrix.h"
#include "pzstepsolver.h"

//------------------PROBLEMA MODELO------------------------

/**
 * @brief Funcao para criar a malha geometrica do problema a ser simulado
 * @note A malha sera unidimensional formada por nel elementos de tamanho elsize
 * @param uNDiv number of divisions ortogonal to the plates performed on the domain
 * @param vNDiv number of divisions parallel to the plates performed on the domain
 * @param nel numero de elementos
 * @param elsize tamanho dos elementos
 */
TPZGeoMesh *CreateGMesh(int64_t nel, int uNDiv, int vNDiv);

/**
 * @brief Funcao para criar a malha computacional do problema a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder);

///Adiciona namespace std
using namespace std;

///Funcao principal do programa
int main(int argc, char *argv[])
{
    int dim = 2;//dimensao do problema
    int uNDiv=3,
    vNDiv=4;//numero de divisoes feitas no dominio
	int nel = uNDiv*vNDiv; //numero de elementos a serem utilizados
    int pOrder = 4; //ordem polinomial de aproximacao
	TPZGeoMesh *gmesh = CreateGMesh(nel, uNDiv, vNDiv); //funcao para criar a malha geometrica
	
	TPZCompMesh *cmesh = CMesh(gmesh, pOrder); //funcao para criar a malha computacional

	// Resolvendo o Sistema
    bool optimizeBandwidth = false; //impede a renumeracao das equacoes do problema(para obter o mesmo resultado do Oden)
	TPZAnalysis an(cmesh, optimizeBandwidth); //cria objeto de analise que gerenciaria a analise do problema
    
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh);
    
    an.SetStructuralMatrix(strmat);
    
    TPZStepSolver<STATE> solver;
    solver.SetDirect(ECholesky);
    an.SetSolver(solver);
    
	an.Run();//assembla a matriz de rigidez (e o vetor de carga) global e inverte o sistema de equacoes
    
	TPZFMatrix<STATE> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
	solucao.Print("Sol",cout,EMathematicaInput);//imprime na formatacao do Mathematica
	
    //fazendo pos processamento para paraview
    TPZStack<string> scalnames, vecnames;
    scalnames.Push("State");//setando para imprimir u
    string plotfile = "ModelProblemSol.vtk";//arquivo de saida que estara na pasta debug
    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
	int postProcessResolution = 3;//define resolucao do pos processamento
    an.PostProcess(postProcessResolution);//realiza pos processamento
    
	std::cout << "FINISHED!" << std::endl;
	
	return 0;
}



TPZGeoMesh *CreateGMesh(int64_t nel, int uNDiv, int vNDiv)
{
	TPZGeoMesh * gmesh = new TPZGeoMesh;//Inicializa objeto da classe TPZGeoMesh
	
	int64_t nnodes = (uNDiv+1)*(vNDiv+1); //numero de nos do problema
	gmesh->NodeVec().Resize(nnodes); //Redimensiona o tamanho do vetor de nos da malha geometrica
	int mat1d = 1; //define id para um material(formulacao fraca)
    int bc0 = -1; //define id para um material(cond contorno esq)
    int bc1 = -2; //define id para um material(cond contorno dir)
    int bc2 = -3; //define id para um material(cond contorno inf)
    int bc3 = -4; //define id para um material(cond contorno sup)
	// Colocando nos na malha
	for (int64_t i = 0 ; i < nnodes; i++) 
	{
        TPZManVector <REAL,3> coord(3,0.);
        coord[0]=8*((float)(i/(uNDiv+1))/vNDiv)+((float)(i%(uNDiv+1))/uNDiv)*16-8;
        coord[0]/=2;
        coord[1]=coord[0]+8-((float)(i%(uNDiv+1))/uNDiv)*16;
		gmesh->NodeVec()[i].SetCoord(coord); //seta coordenada de um no no vetor de nos da malha
		gmesh->NodeVec()[i].SetNodeId(i); //atribui identificacao para um no
	}
	
	// Criando Elementos
    TPZManVector<int64_t,4> topolQuad(4,0.); //vetor que será inicializado com o indice dos nos de um elemento quadrilateral
	TPZManVector <int64_t,4> topolLine(2,0.); //vetor que sera inicializado com o indice dos nos de um elemento unidimensional
    TPZVec <int64_t> TopolPoint(1); //vetor que sera inicializado com o indice do no de um elemento zero-dimensional
	int64_t id; //id do elemento que sera preenchido pelo metodo CreateGeoElement
    int64_t column, row;
	for (int64_t iel = 0; iel < nel; iel++) 
	{
        column=iel%(uNDiv);
        row=iel/(uNDiv);
        topolQuad[0]=(row)*(uNDiv+1)+column;//Lower left vertex
        topolQuad[1]=(row)*(uNDiv+1)+column+1;//Lower right vertex
        topolQuad[2]=(row+1)*(uNDiv+1)+column+1;//Upper right vertex
        topolQuad[3]=(row+1)*(uNDiv+1)+column;//Upper left vertex
        
        //column==0 <=> left side
        // Cond Contorno esquerda
        if(column==0)
        {
            topolLine[0] = topolQuad[0];
            topolLine[1] = topolQuad[3];
            gmesh->CreateGeoElement(EOned, topolLine, bc0, id);
        }
        //column==uNDiv-1 <==> right side
        // Cond Contorno direita
        if(column==uNDiv-1)
        {
            topolLine[0] = topolQuad[1];
            topolLine[1] = topolQuad[2];
            gmesh->CreateGeoElement(EOned, topolLine, bc1, id);
        }
        //row==0<==>lower plate(V=0)
        // Cond Contorno inferior
        if(row==0)
        {
            topolLine[0] = topolQuad[0];
            topolLine[1] = topolQuad[1];
            gmesh->CreateGeoElement(EOned, topolLine, bc2, id);
        }
        //row==vNDiv-1 <==> upper plate (V=1)
        // Cond Contorno superior
        if(row==vNDiv-1)
        {
            topolLine[0] = topolQuad[2];
            topolLine[1] = topolQuad[3];
            gmesh->CreateGeoElement(EOned, topolLine, bc3, id);
        }
        gmesh->CreateGeoElement(EQuadrilateral, topolQuad, mat1d, id);//cria elemento quadrilateral
        gmesh->ElementVec()[id];
	}
    gmesh->BuildConnectivity(); //constroi a conectividade de vizinhanca da malha
#ifdef PZDEBUG
    std::ofstream out("geomesh.vtk"), outtxt("gmesh.txt");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out,true);//imprimindo a malha geometrica no vtk
    gmesh->Print(outtxt);
#endif
    
    
    
	return gmesh;
}

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder)
{
	const int dim = 2; //dimensao do problema
	const int matId = 1, bc0 = -1, bc1 = -2, bc2=-3, bc3=-4; //MESMOS ids da malha geometrica
    const int dirichlet = 0, neumann = 1;
//    const int mixed = 2; //tipo da condicao de contorno do problema ->default dirichlet na esquerda e na direita
	
    
	///criar malha computacional
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
	cmesh->SetDimModel(dim);//seta dimensao do modelo
	
    // Criando material
    TPZMatExSimples2D *material = new TPZMatExSimples2D(matId);//criando material que implementa a formulacao fraca do problema modelo
    
	// Inserindo material na malha
	cmesh->InsertMaterialObject(material);
		
	///Inserir condicao de contorno esquerda
	TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(material, bc0, neumann, val1, val2);//cria material que implementa a condicao de contorno da esquerda
	
    cmesh->InsertMaterialObject(BCond0);//insere material na malha
    
	// Condicao de contorno da direita
	TPZMaterial * BCond1 = material->CreateBC(material, bc1, neumann, val1, val2);//cria material que implementa a condicao de contorno da direita
    
    cmesh->InsertMaterialObject(BCond1);//insere material na malha
    
    val2(0,0) = 1.0;//potencial na placa inferior
    // Condicao de contorno da placa inferior
    TPZMaterial * BCond2 = material->CreateBC(material, bc2, dirichlet, val1, val2);//cria material que implementa a condicao de contorno da placa inferior
    
    cmesh->InsertMaterialObject(BCond2);//insere material na malha
    
    val2(0,0) = 1.5;//potencial na placa superior
    // Condicao de contorno da placa superior
    TPZMaterial * BCond3 = material->CreateBC(material, bc3, dirichlet, val1, val2);//cria material que implementa a condicao de contorno da placa superior
	
    cmesh->InsertMaterialObject(BCond3);//insere material na malha
	
	//Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
	cmesh->AutoBuild();
	
	return cmesh;
	
}
