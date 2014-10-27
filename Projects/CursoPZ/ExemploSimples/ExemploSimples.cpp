/**
 * @file
 * @brief Example to Run Model Problem from Oden's Chapter I in a NeoPZ Environment
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
#include "TPZMatModelProblem.h"

//------------------PROBLEMA MODELO------------------------

/**
 * @brief Funcao para criar a malha geometrica do problema a ser simulado
 * @note A malha sera unidimensional formada por nel elementos de tamanho elsize
 * @param nel numero de elementos
 * @param elsize tamanho dos elementos
 */
TPZGeoMesh *CreateGMesh(long nel, REAL elsize);

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
    int dim = 1;//dimensao do problema
	REAL dom = 1.; //comprimento do dominio unidimensional com inicio na origem zero
	int nel = 4; //numero de elementos a serem utilizados
    int pOrder = 1; //ordem polinomial de aproximacao
	REAL elsize = dom/nel; //tamanho de cada elemento
	TPZGeoMesh *gmesh = CreateGMesh(nel, elsize); //funcao para criar a malha geometrica
	
	std::ofstream out("gmesh.vtk"); //define arquivo de saida para impressao da malha no paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true); //imprime a malha no formato vtk
	
	TPZCompMesh *cmesh = CMesh(gmesh, pOrder); //funcao para criar a malha computacional

	// Resolvendo o Sistema
    bool optimizeBandwidth = false; //impede a renumeracao das equacoes do problema(para obter o mesmo resultado do Oden)
	TPZAnalysis an(cmesh, optimizeBandwidth); //cria objeto de analise que gerenciaria a analise do problema
	an.Run();//assembla a matriz de rigidez (e o vetor de carga) global e inverte o sistema de equacoes
    
	TPZFMatrix<REAL> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
	solucao.Print("Sol",cout,EMathematicaInput);//imprime na formatacao do Mathematica
	
    //fazendo pos processamento para paraview
    TPZStack<string> scalnames, vecnames;
    scalnames.Push("State");//setando para imprimir u
    string plotfile= "ModelProblemSol.vtk";//arquivo de saida que estara na pasta debug
    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
	int postProcessResolution = 0;//define resolucao do pos processamento
    an.PostProcess(postProcessResolution);//realiza pos processamento
    
	std::cout << "FINISHED!" << std::endl;
	
	return 0;
}



TPZGeoMesh *CreateGMesh(long nel, REAL elsize)
{
	TPZGeoMesh * gmesh = new TPZGeoMesh;//Inicializa objeto da classe TPZGeoMesh
	
	long nnodes = nel + 1; //numero de nos do problema
	gmesh->NodeVec().Resize(nnodes); //Redimensiona o tamanho do vetor de nos da malha geometrica
	int mat1d = 1; //define id para um material(formulacao fraca)
    int bc0 = -1; //define id para um material(cond contorno esq)
    int bc1 = -2; //define id para um material(cond contorno dir)
	
	// Colocando nos na malha
	for (long i = 0 ; i < nnodes; i++) 
	{
		const REAL pos = i * elsize;
		TPZVec <REAL> coord(3,0.);
		coord[0] = pos;
		gmesh->NodeVec()[i].SetCoord(coord); //seta coordenada de um no no vetor de nos da malha
		gmesh->NodeVec()[i].SetNodeId(i); //atribui identificacao para um no
	}
	
	// Criando Elementos
	TPZVec <long> topol(2); //vetor que sera inicializado com o indice dos nos de um elemento unidimensional
    TPZVec <long> TopolPoint(1); //vetor que sera inicializado com o indice do no de um elemento zero-dimensional
	long id; //id do elemento que sera preenchido pelo metodo CreateGeoElement
	
	for (long iel = 0; iel < nel; iel++) 
	{
		const long ino1 = iel;
		const long ino2 = iel + 1;
		topol[0] = ino1;
		topol[1] = ino2;
        gmesh->CreateGeoElement(EOned, topol, mat1d, id);//cria elemento unidimensional
        gmesh->ElementVec()[id];
	}
	
	// Cond Contorno esquerda
	TopolPoint[0] = 0;
    gmesh->CreateGeoElement(EPoint, TopolPoint, bc0, id);

	// Cond Contorno Direita
	TopolPoint[0] = nnodes-1;
	gmesh->CreateGeoElement(EPoint, TopolPoint, bc1, id);
	
	gmesh->BuildConnectivity(); //constroi a conectividade de vizinhanca da malha
	
	return gmesh;
}

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder)
{
	const int dim = 1; //coeficientes da eq diferencial unidimensional
	const int matId = 1, bc0 = -1, bc1 = -2; //MESMOS ids da malha geometrica
	const int dirichlet = 0, neumann = 1, mixed = 2; //tipo da condicao de contorno do problema ->default dirichlet na esquerda e na direita 
	
	// Criando material
	TPZMatModelProblem *material = new TPZMatModelProblem(matId);//criando material que implementa a formulacao fraca do problema modelo
	
	///criar malha computacional
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
	cmesh->SetDimModel(dim);//seta dimensao do modelo
	
	// Inserindo material na malha
	cmesh->InsertMaterialObject(material);
		
	///Inserir condicao de contorno esquerda
	TPZFMatrix<REAL> val1(1,1,0.), val2(1,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(material, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno da esquerda
	
	// Condicao de contorno da direita
	TPZMaterial * BCond1 = material->CreateBC(material, bc1, dirichlet, val1, val2);//cria material que implementa a condicao de contorno da direita
	
	cmesh->InsertMaterialObject(BCond0);//insere material na malha
	cmesh->InsertMaterialObject(BCond1);//insere material na malha
	
	//Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
	cmesh->AutoBuild();
	
	return cmesh;
	
}