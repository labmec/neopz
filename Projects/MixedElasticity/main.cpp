#include <fstream>

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include <TPZGeoElement.h>

#include "pzmat2dlin.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"

#include <pz_config.h>

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoblend.h"

#include "TPZMatElasticity2D.h"
#include "TPZInterfaceEl.h"
#include "pzdiscgal.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include <tpzarc3d.h>

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzl2projection.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"
#include "meshgen.h"

#include "pzmixedelasmat.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSpStructMatrix.h"
#include "pzlog.h"
#include <iostream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzfunction.h"
#include "TPZReadGIDGrid.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"


#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
#endif

//------------------Problema Elasticidade------------------------


/**
 * @brief Funcao para criar a malha geometrica do problema a ser simulado
 * @note A malha sera unidimensional formada por nel elementos de tamanho elsize
 * @param uNDiv number of divisions ortogonal to the plates performed on the domain
 * @param vNDiv number of divisions parallel to the plates performed on the domain
 * @param nel numero de elementos
 * @param elsize tamanho dos elementos
 */
TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy);

/**
 * @brief Funcao para criar a malha computacional da velocidade a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_S(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional da pressão a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_U(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional da pressão a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_P(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional multi-fisica ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder);


//Função para criar interface entre elmentos:

TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);


void Error(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);

//Variáveis globais do problema:

const int dim = 2; //Dimensão do problema
const int matID = 1; //Materia do elemento volumétrico
const int matBCbott = -1, matBCtop = -2, matBCleft = -3, matBCright = -4; //Materiais das condições de contorno
const int matInterface = 4; //Material do elemento de interface
const int matIntBCbott=-11, matIntBCtop=-12, matIntBCleft=-13, matIntBCright=-14; //Materiais das condições de contorno (elementos de interface)
const int matPoint =-5;//Materia de um ponto
const int dirichlet = 0, neumann = 1, mixed = 2, pointtype=5, dirichletvar=4; //Condições de contorno do problema ->default Dirichlet na esquerda e na direita
const int pointtypex = 3, pointtypey = 4;
const REAL visco=1., permeability=1., theta=-1.; //Coeficientes: viscosidade, permeabilidade, fator simetria
const REAL Pi=M_PI;

const int quadmat1 = 1; //Parte inferior do quadrado
const int quadmat2 = 2; //Parte superior do quadrado
const int quadmat3 = 3; //Material de interface


void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget);


using namespace std;

//teste 0
// definition of f
void f_source(const TPZVec<REAL> & x, TPZVec<STATE>& f, TPZFMatrix<STATE>& gradu);
// definition of v analytic
void v_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f);
// definition of p analytic
void p_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f);
// definition of sol analytic
void sol_exact(const TPZVec<REAL> & x, TPZVec<STATE>& p, TPZFMatrix<STATE>& v);


//Função principal do programa:

int main(int argc, char *argv[])
{
    
    TPZMaterial::gBigNumber = 1.e16;
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    //Dados do problema:
    TElasticityExample1 test;
    TPZManVector<STATE,2> displ(2), force(2);
    TPZFNMatrix<4,STATE> sigma(2,2);
    TPZManVector<REAL,3> x(3,1.);
	TPZManVector<STATE, 3> x_state(3, (STATE)1.);
    test.Sigma(x_state,sigma);
    test.Force(x,force);
    int h_level = 1;
    
    
    double hx=2.,hy=2.; //Dimensões em x e y do domínio
    int nelx=h_level, nely=h_level; //Número de elementos em x e y
    int nx=nelx+1 ,ny=nely+1; //Número de nos em x  y
    int pOrder = 1; //Ordem polinomial de aproximação
    //double elsizex=hx/nelx, elsizey=hy/nely; //Tamanho dos elementos
    //int nel = elsizex*elsizey; //Número de elementos a serem utilizados
    
    //Gerando malha geométrica:
    
    TPZGeoMesh *gmesh = CreateGMesh(nx, ny, hx, hy); //Função para criar a malha geometrica
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Gerando malha computacional:
    
    TPZCompMesh *cmesh_S = CMesh_S(gmesh, pOrder); //Função para criar a malha computacional da tensão
    TPZCompMesh *cmesh_U = CMesh_U(gmesh, pOrder); //Função para criar a malha computacional da deslocamento
    TPZCompMesh *cmesh_P = CMesh_P(gmesh, pOrder); //Função para criar a malha computacional da rotação
    TPZCompMesh *cmesh_m = CMesh_m(gmesh, pOrder); //Função para criar a malha computacional multifísica
    
#ifdef PZDEBUG
    {
        std::ofstream filecS("MalhaC_S.txt"); //Impressão da malha computacional da tensão (formato txt)
        std::ofstream filecU("MalhaC_U.txt"); //Impressão da malha computacional da deslocamento (formato txt)
        std::ofstream filecP("MalhaC_P.txt"); //Impressão da malha computacional da rotação (formato txt)
        cmesh_S->Print(filecS);
        cmesh_U->Print(filecU);
        cmesh_P->Print(filecP);
    }
#endif
    
    TPZManVector<TPZCompMesh *, 2> meshvector(3);
    meshvector[0] = cmesh_S;
    meshvector[1] = cmesh_U;
    meshvector[2] = cmesh_P;
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    cmesh_m->LoadReferences();
    
    //    AddMultiphysicsInterfaces(*cmesh_m,matInterface,matID);
    //    AddMultiphysicsInterfaces(*cmesh_m,matIntBCbott,matBCbott);
    //    AddMultiphysicsInterfaces(*cmesh_m,matIntBCtop,matBCtop);
    //    AddMultiphysicsInterfaces(*cmesh_m,matIntBCleft,matBCleft);
    //    AddMultiphysicsInterfaces(*cmesh_m,matIntBCright,matBCright);
    
#ifdef PZDEBUG
    std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
    gmesh->Print(fileg1);
    
    std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
    cmesh_m->Print(filecm);
#endif
    
    //Resolvendo o Sistema:
    int numthreads = 0;
    
    bool optimizeBandwidth = false; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    TPZFStructMatrix matskl(cmesh_m); //caso nao simetrico ***
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    
    
    std::cout << "Assemble matrix with NDoF = " << cmesh_m->NEquations() << std::endl;
    
    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
    
    
#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    {
        std::ofstream filestiff("stiffness.txt");
        an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
        
        std::ofstream filerhs("rhs.txt");
        an.Rhs().Print("R = ",filerhs,EMathematicaInput);
    }
#endif
    
    std::cout << "Solving Matrix " << std::endl;
    
    an.Solve();
    
    {
        std::ofstream file("file.txt");
        an.Solution().Print("sol=",file,EMathematicaInput);
        
    }
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector,cmesh_m);

    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("SigmaX");
    vecnames.Push("Flux");
    vecnames.Push("displacement");
    vecnames.Push("Stress");
    an.DefineGraphMesh(2, scalnames, vecnames, "mixedelast.vtk");
    an.PostProcess(0);
    
#ifdef PZDEBUG
    //Imprimindo vetor solução:
    {
        TPZFMatrix<STATE> solucao=cmesh_m->Solution();//Pegando o vetor de solução, alphaj
        std::ofstream solout("sol.txt");
        solucao.Print("Sol",solout,EMathematicaInput);//Imprime na formatação do Mathematica
        
        std::ofstream fileAlpha("alpha.txt");
        an.Solution().Print("Alpha = ",fileAlpha,EMathematicaInput);
    }
#endif
    
    
    //    //Calculo do erro
    //    std::cout << "Comuting Error " << std::endl;
    //    TPZManVector<REAL,3> Errors;
    //    ofstream ErroOut("Erro.txt");
    //    an.SetExact(sol_exact);
    //    an.PostProcessError(Errors,ErroOut);
    //
    //
    //
    //    //Pós-processamento (paraview):
    //    std::cout << "Post Processing " << std::endl;
    //    std::string plotfile("ElasticityTest.vtk");
    //    TPZStack<std::string> scalnames, vecnames;
    //    vecnames.Push("Displacement");
    //    vecnames.Push("Stress");
    //    vecnames.Push("Rotation");
    ////    vecnames.Push("V_exact");
    ////    vecnames.Push("P_exact");
    //    //        vecnames.Push("V_exactBC");
    //
    //
    //    int postProcessResolution = 3; //  keep low as possible
    //
    //    int dim = gmesh->Dimension();
    //    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    //    an.PostProcess(postProcessResolution,dim);
    
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
}

// Teste 0
// definition of f
void f_source(const TPZVec<REAL> & x, TPZVec<STATE>& f, TPZFMatrix<STATE>& gradu){
    f.resize(1);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE f_x = 8.0*Pi*Pi*sin(2.0*Pi*xv)*sin(2.0*Pi*yv);
    
    f[0] = f_x;
}

// definition of v analytic
void v_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    
    f.resize(2);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE v_x =  -2.*Pi*cos(2.*Pi*xv)*sin(2.*Pi*yv);
    STATE v_y =  -2.*Pi*cos(2.*Pi*yv)*sin(2.*Pi*xv);
    
    f[0] = v_x; // x direction
    f[1] = v_y; // y direction
}

// definition of p analytic
void p_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    
    f.resize(1);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE p_x = sin(2.*Pi*xv)*sin(2.*Pi*yv);
    
    f[0] = p_x;
}

// BC - Solução analítica - Artigo
void solucao_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    
    f.resize(3);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE v_x =  -2.*Pi*cos(2.*Pi*xv)*sin(2.*Pi*yv);
    STATE v_y =  -2.*Pi*cos(2.*Pi*yv)*sin(2.*Pi*xv);
    STATE p =  sin(2.*Pi*xv)*sin(2.*Pi*yv);
    
    f[0] = v_x; // x direction
    f[1] = v_y; // y direction
    f[2] = p; //
}


// Solução analítica - Artigo
void sol_exact(const TPZVec<REAL> & x, TPZVec<STATE>& sol, TPZFMatrix<STATE>& dsol){
    
    dsol.Resize(3,2);
    sol.Resize(3);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE v_x = -2.*Pi*cos(2.*Pi*xv)*sin(2.*Pi*yv);
    STATE v_y = -2.*Pi*cos(2.*Pi*yv)*sin(2.*Pi*xv);
    STATE pressure= sin(2.*Pi*xv)*sin(2.*Pi*yv);
    
    sol[0]=v_x;
    sol[1]=v_y;
    sol[2]=pressure;
    
    // vx direction
    dsol(0,0)= 4.*Pi*Pi*sin(2.*Pi*xv)*sin(2.*Pi*yv);
    dsol(0,1)= -4.*Pi*Pi*cos(2.*Pi*xv)*cos(2.*Pi*yv);
    
    // vy direction
    dsol(1,0)= -4.*Pi*Pi*cos(2.*Pi*xv)*cos(2.*Pi*yv);
    dsol(1,1)= 4.*Pi*Pi*sin(2.*Pi*xv)*sin(2.*Pi*yv);
    
    // Gradiente pressão
    
    dsol(2,0)= -v_x;
    dsol(2,1)= -v_y;
    
    
}



TPZGeoMesh *CreateGMesh(int nx, int ny, double hx, double hy)
{
    
    int i,j;
    int64_t id, index;
    
    
    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    
    //Vetor auxiliar para armazenar coordenadas:
    
    TPZVec <REAL> coord (3,0.);
    
    
    //Inicialização dos nós:
    
    for(i = 0; i < ny; i++){
        for(j = 0; j < nx; j++){
            id = i*nx + j;
            coord[0] = (j)*hx/(nx - 1)-1;
            coord[1] = (i)*hy/(ny - 1)-1;
            //using the same coordinate x for z
            coord[2] = 0.;
            //cout << coord << endl;
            //Get the index in the mesh nodes vector for the new node
            index = gmesh->NodeVec().AllocateNewElement();
            //Set the value of the node in the mesh nodes vector
            gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
        }
    }
    
    //Ponto 1
    //    TPZVec<int64_t> pointtopology(1);
    //    pointtopology[0] = 0;
    //
    //    gmesh->CreateGeoElement(EPoint,pointtopology,matPoint,id);
    
    
    //Vetor auxiliar para armazenar as conecções entre elementos:
    
    TPZVec <int64_t> connect(4,0);
    
    
    //Conectividade dos elementos:
    
    for(i = 0; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            index = (i)*(nx - 1)+ (j);
            connect[0] = (i)*ny + (j);
            connect[1] = connect[0]+1;
            connect[2] = connect[1]+(nx);
            connect[3] = connect[0]+(nx);
            gmesh->CreateGeoElement(EQuadrilateral,connect,matID,id);
            if(id != index) DebugStop();
        }
    }
    
    
    //Gerando informação da vizinhança:
    
    gmesh->BuildConnectivity();
    
    for(i = 0; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            index = (i)*(nx - 1)+ (j);
            TPZGeoEl *gel = gmesh->Element(index);
            if (j==0) {
                TPZGeoElBC gbcl(gel,7,matBCleft);
            }
            if(j == nx-2)
            {
                TPZGeoElBC gbcr(gel,5,matBCright);
            }
            if (i==0) {
                TPZGeoElBC gbcb(gel,4,matBCbott);
            }
            if (i==ny-2) {
                TPZGeoElBC gbct(gel,6,matBCtop);
            }
        }
    }

    {
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
    }
    
    //Impressão da malha geométrica:
    
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;
    
    
    
}

TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}


TPZCompMesh *CMesh_S(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim);//Insere dimensão do modelo
    
    
    //Definição do espaço de aprximação:
    
    cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
    
    
    //Criando material cujo nSTATE = 2:
    
    REAL E = 1.; //* @param E elasticity modulus
    REAL nu=0.; //* @param nu poisson coefficient
    REAL fx=0.;//* @param fx forcing function \f$ -x = fx \f$
    REAL fy=0.;//* @param fx forcing function \f$ -x = fx \f$
    int plain=0.; //* @param plainstress = 1 \f$ indicates use of plainstress
    
    TPZMaterial * material = new TPZElasticityMaterial(matID,E,nu,fx,fy,plain,dim);
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
    
    //Dimensões do material (para H1 e descontinuo):
    //TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
    //material->SetMaterial(xkin, xcin, xfin);
    
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZFMatrix<STATE> val2s(2,1,0.);
    val2s(0,0) = 10.0; // vx -> 0
    val2s(1,0) = 0.0; // vy -> 0
    
    TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, neumann, val1, val2); //Cria material que implementa a condição de contorno inferior
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, neumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    TPZMaterial * BCond2 = material->CreateBC(material, matBCleft, neumann, val1, val2s); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    TPZMaterial * BCond3 = material->CreateBC(material, matBCright, neumann, val1, val2s); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

TPZCompMesh *CMesh_U(TPZGeoMesh *gmesh, int pOrder)
{
    
    // @omar::
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    
    
    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    //Criando material cujo nSTATE = 2:
    
    REAL E = 1.; //* @param E elasticity modulus
    REAL nu=0.; //* @param nu poisson coefficient
    REAL fx=0.;//* @param fx forcing function \f$ -x = fx \f$
    REAL fy=0.;//* @param fx forcing function \f$ -x = fx \f$
    int plain=0.; //* @param plainstress = 1 \f$ indicates use of plainstress
    
    TPZMaterial * material = new TPZElasticityMaterial(matID,E,nu,fx,fy,plain,dim);
    
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
    //Dimensões do material (para H1 e descontinuo):
    //    TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
    //    material->SetMaterial(xkin, xcin, xfin);
    
    //Condições de contorno:
    
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0; // Ux -> 0
    val2(1,0) = 0.0; // Uy -> 0
    
    TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    TPZMaterial * BCond2 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    TPZMaterial * BCond3 = material->CreateBC(material, matBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    //    //    Ponto de pressao:
    //    //
    //    TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
    //    ////
    //    TPZMaterial * BCPoint = material->CreateBC(material, matPoint, pointtype, val3, val4); //Cria material que implementa um ponto para a pressao
    //    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
    
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    std::set<int> materialids;
    materialids.insert(matID);
    cmesh->AutoBuild(materialids);
    cmesh->LoadReferences();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild();
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    //    cmesh->AdjustBoundaryElements();
    //    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

TPZCompMesh *CMesh_P(TPZGeoMesh *gmesh, int pOrder)
{
    
    // @omar::
    
    //pOrder--; // Space restriction apapapa
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    
    // @omar::
    //cmesh->SetAllCreateFunctionsDiscontinuous();
    
    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    
    //Criando material cujo nSTATE = 1:
    
    TPZMaterial *material = new TPZMatPoisson3d(matID,dim);//criando material que implementa a formulacao fraca do problema modelo
    
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
    //    //Dimensões do material (para H1 e descontínuo):
    //    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    //    material->SetMaterial(xkin, xcin, xfin);
    
    //Condições de contorno:
    
    
    //    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    //
    //        val2(0,0) = 0.0; // px -> 0
    //        val2(1,0) = 0.0; // py -> 0
    //
    //        TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //        cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    //
    //        TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    //        cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    //
    //        TPZMaterial * BCond2 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //        cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    //
    //        TPZMaterial * BCond3 = material->CreateBC(material, matBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //        cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    //    Ponto de pressao:
    //
    //    TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
    //    ////
    //    TPZMaterial * BCPoint = material->CreateBC(material, matPoint, pointtype, val3, val4); //Cria material que implementa um ponto para a pressao
    //    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
    //
    
    
    std::set<int> materialids;
    materialids.insert(matID);
    cmesh->AutoBuild(materialids);
    cmesh->LoadReferences();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild();
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    //    cmesh->AdjustBoundaryElements();
    //    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}


TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    int bc_inte_order = 10;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    TElasticityExample1 example;
    
    // Criando material:
    
    REAL E = 1.; //* @param E elasticity modulus
    REAL nu=0.; //* @param nu poisson coefficient
    
    
    REAL fx=0.;//* @param fx forcing function \f$ -x = fx \f$
    REAL fy=0.;//* @param fx forcing function \f$ -x = fx \f$
    int plain=0.; //* @param plainstress = 1 \f$ indicates use of plainstress
    
    TPZMaterial * material = new TPZElasticityMaterial(matID,E,nu,fx,fy,plain,dim);
    
    
    // Inserindo material na malha
    //    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (f_source);
    //    TPZAutoPointer<TPZFunction<STATE> > pp = new TPZDummyFunction<STATE> (p_exact);
    //    TPZAutoPointer<TPZFunction<STATE> > vp = new TPZDummyFunction<STATE> (v_exact);
    
    //    material->SetForcingFunction(fp);
    //    material->SetForcingFunctionExactPressure(pp);
    //    material->SetForcingFunctionExact(vp);
    
    cmesh->InsertMaterialObject(material);
    
    
    //Condições de contorno:
    

    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZFMatrix<STATE> val2s(2,1,0.);
    val2s(0,0) = 10.0; // vx -> 0
    val2s(1,0) = 0.0; // vy -> 0
    
    TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, neumann, val1, val2); //Cria material que implementa a condição de contorno inferior

    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    //BCond0->SetForcingFunction(solucao_exact,bc_inte_order);
    BCond0->SetForcingFunction(example.ValueFunction());
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    BCond1->SetForcingFunction(example.ValueFunction());
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond1->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    TPZMaterial * BCond2 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    BCond2->SetForcingFunction(example.ValueFunction());
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond2->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    TPZMaterial * BCond3 = material->CreateBC(material, matBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond3->SetForcingFunction(example.ValueFunction());
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond3->SetForcingFunction(solucao_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    //Ponto
    
    //    TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
    //    val4(0,0)=0.0;
    //
    //    TPZMaterial * BCPoint = material->CreateBC(material, matPoint, pointtype, val3, val4); //Cria material que implementa um ponto para a pressão
    //    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
    
    
    
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}


void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->MaterialId() != matfrom) {
            continue;
        }
        
        int nsides= gel->NSides();
        
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        if (celstack.size() != 2) {
            DebugStop();
        }
        gel->SetMaterialId(mattarget);
        int64_t index;
        new TPZMultiphysicsInterfaceElement(cmesh,gel,index,celstack[1],celstack[0]);
    }
    
}

//Função Erro (Não utilizada)
void Error(TPZCompMesh *cmesh, std::ostream &out, int p, int ndiv)
{
    DebugStop();
    int64_t nel = cmesh->NElements();
    //int dim = cmesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(sol_exact, elerror, false);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with HDiv space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
    out << "L2 Norm for flux - L2 Norm for divergence - Hdiv Norm for flux " << endl;
    out <<  setw(16) << sqrt(globalerrors[1]) << setw(25)  << sqrt(globalerrors[2]) << setw(21)  << sqrt(globalerrors[3]) << endl;
    //
    //    out << "L2 Norm for flux = "    << sqrt(globalerrors[1]) << endl;
    //    out << "L2 Norm for divergence = "    << sqrt(globalerrors[2])  <<endl;
    //    out << "Hdiv Norm for flux = "    << sqrt(globalerrors[3])  <<endl;
    //
}


