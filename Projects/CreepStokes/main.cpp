

#include <cmath>
#include <set>
#include <math.h>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "TPZStokesMaterial.h"
#include "TPZDarcyPMaterial.h"
#include "TPZCouplingDSMaterial.h"

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiPhysicsInterfaceEl.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"


//------------------Darcy-Stokes Coupling Simulation------------------------




/**
 * @brief Funcao para criar a malha geometrica do problema a ser simulado
 * @note A malha sera unidimensional formada por nel elementos de tamanho elsize
 * @param uNDiv number of divisions ortogonal to the plates performed on the domain
 * @param vNDiv number of divisions parallel to the plates performed on the domain
 * @param nel numero de elementos
 * @param elsize tamanho dos elementos
 */
TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hs, double hd);

/**
 * @brief Funcao para criar a malha computacional da velocidade a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional da pressão a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder);

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

const int matIdS = 1;
const int matIdD = 2;

const int matSBCbott = -1, matSBCtop = -2, matSBCleft = -3, matSBCright = -4; //Materiais das condições de contorno (Stokes)
const int matDBCbott = -21, matDBCtop = -22, matDBCleft = -23, matDBCright = -24; //Materiais das condições de contorno (Darcy)

const int matSInterface = 4,matDInterface = 14, matInterfaceDS=5; //Material do elemento de interface

const int matIntSBCbott=-11, matIntSBCtop=-12, matIntSBCleft=-13, matIntSBCright=-14; //Materiais das condições de contorno (elementos de interface - Stokes)
const int matIntDBCbott=-31, matIntDBCtop=-32, matIntDBCleft=-33, matIntDBCright=-34; //Materiais das condições de contorno (elementos de interface - Darcy)

const int matPoint =-5, matPoint2 =-6;//Materia de um ponto
const int dirichlet = 0, neumann = 1, penetration = 2, pointtype=5, dirichletvar=4; //Condições de contorno do problema ->default Dirichlet na esquerda e na direita
const int pointtypex = 3, pointtypey = 4;
const REAL visco=1.6, permeability=0.00005420731707, theta=-1.; //Coeficientes: viscosidade, permeabilidade, fator simetria
const REAL Pi=M_PI;

const int quadmat1 = 1; //Parte inferior do quadrado
const int quadmat2 = 2; //Parte superior do quadrado
const int quadmat3 = 3; //Material de interface


void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget);

void AddInterfaceCoupllingDS(TPZGeoMesh *gmesh,TPZCompMesh *cmesh, int matInterfaceDS ,int matleft, int matright);

using namespace std;

// Solução Exata Stokes
void fS_source(const TPZVec<REAL> & x, TPZVec<STATE>& f);
// definition of v analytic
void vS_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f);
// definition of p analytic
void pS_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f);
// definition of sol analytic
void sol_exact(const TPZVec<REAL> & x, TPZVec<STATE>& p, TPZFMatrix<STATE>& v);


// Solução Exata Darcy
void fS_source(const TPZVec<REAL> & x, TPZVec<STATE>& f);
// definition of v analytic
void vS_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f);
// definition of p analytic
void pS_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f);
// definition of sol analytic
//void sols_exact(const TPZVec<REAL> & x, TPZVec<STATE>& p, TPZFMatrix<STATE>& v);



//Função principal do programa:

int main(int argc, char *argv[])
{
    
    TPZMaterial::gBigNumber = 1.e16;
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    //Dados do problema:
    
    
    int h_level = 1;
    
    
    double hx=0.4318,hs=0.7747, hd=0.2667; //Dimensões em x e y do domínio
    int nelx=h_level, nely=h_level*2; //Número de elementos em x e y
    int nx=nelx+1 ,ny=nely+1; //Número de nos em x  y
    int pOrder = 1; //Ordem polinomial de aproximação
    //double elsizex=hx/nelx, elsizey=hy/nely; //Tamanho dos elementos
    //int nel = elsizex*elsizey; //Número de elementos a serem utilizados
    
    //Gerando malha geométrica:
    
    TPZGeoMesh *gmesh = CreateGMesh(nx, ny, hx, hs, hd); //Função para criar a malha geometrica
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Gerando malha computacional:
    
    TPZCompMesh *cmesh_v = CMesh_v(gmesh, pOrder); //Função para criar a malha computacional da velocidade
    TPZCompMesh *cmesh_p = CMesh_p(gmesh, pOrder); //Função para criar a malha computacional da pressão
    TPZCompMesh *cmesh_m = CMesh_m(gmesh, pOrder); //Função para criar a malha computacional multifísica
    
#ifdef PZDEBUG
    {
        std::ofstream filecv("MalhaC_v.txt"); //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecp("MalhaC_p.txt"); //Impressão da malha computacional da pressão (formato txt)
        cmesh_v->Print(filecv);
        cmesh_p->Print(filecp);
    }
#endif
    
    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_v;
    meshvector[1] = cmesh_p;
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    cmesh_m->LoadReferences();
    
    
    AddMultiphysicsInterfaces(*cmesh_m,matSInterface,matIdS);
    AddMultiphysicsInterfaces(*cmesh_m,matDInterface,matIdD);
    
    
    //AddInterfaceCoupllingDS(gmesh,cmesh_m,matInterfaceDS ,matIdS,matIdD);
    AddInterfaceCoupllingDS(gmesh,cmesh_m,matInterfaceDS ,matIdD,matIdS);
    
    //AddMultiphysicsInterfaces(*cmesh_m,matIntSBCbott,matSBCbott);
    AddMultiphysicsInterfaces(*cmesh_m,matIntSBCtop,matSBCtop);
    AddMultiphysicsInterfaces(*cmesh_m,matIntSBCleft,matSBCleft);
    AddMultiphysicsInterfaces(*cmesh_m,matIntSBCright,matSBCright);
    
    AddMultiphysicsInterfaces(*cmesh_m,matIntDBCbott,matDBCbott);
    //AddMultiphysicsInterfaces(*cmesh_m,matIntDBCtop,matDBCtop);
    AddMultiphysicsInterfaces(*cmesh_m,matIntDBCleft,matDBCleft);
    AddMultiphysicsInterfaces(*cmesh_m,matIntDBCright,matDBCright);
    
    
#ifdef PZDEBUG
    std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
    gmesh->Print(fileg1);
    
    std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
    cmesh_m->Print(filecm);
#endif
    
    //Resolvendo o Sistema:
    int numthreads = 0;
    
    bool optimizeBandwidth = true; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    TPZSkylineNSymStructMatrix matskl(cmesh_m); //caso nao simetrico ***
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
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
    
    //    REAL a[] = {-4.94792e-18, 4.16667e-18, -4.94792e-18,-3.19038e-33,4.94792e-18,3.14852e-32,4.94792e-18,4.16667e-18, -4.94792e-18, -4.16667e-18, 4.94792e-18,-4.16667e-18,1.54074e-49, 0.216146, -0.216146, 1.0689e-15,-0.208333,-0.0078125, 0.0078125, 0.208333};
    //    int size = an.Solution().Rows();
    //    for (int i = 0; i < size; i++) {
    //        an.Solution()(i,0) = a[i];
    //    }
    //    an.LoadSolution();
    
    
#ifdef PZDEBUG
    //Imprimindo vetor solução:
    {
        TPZFMatrix<REAL> solucao=cmesh_m->Solution();//Pegando o vetor de solução, alphaj
        std::ofstream solout("sol.txt");
        solucao.Print("Sol",solout,EMathematicaInput);//Imprime na formatação do Mathematica
        
        std::ofstream fileAlpha("alpha.txt");
        an.Solution().Print("Alpha = ",fileAlpha,EMathematicaInput);
    }
#endif
    
    //Calculo do erro
    std::cout << "Comuting Error " << std::endl;
    TPZManVector<REAL,3> Errors;
    ofstream ErroOut("Erro.txt");
    an.SetExact(sol_exact);
    an.PostProcessError(Errors,ErroOut);
    
    
    
    //Pós-processamento (paraview):
    std::cout << "Post Processing " << std::endl;
    std::string plotfile("CouplingDS.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("P");
    vecnames.Push("V");
    vecnames.Push("f");
    vecnames.Push("V_exact");
    vecnames.Push("P_exact");
    //        vecnames.Push("V_exactBC");
    
    
    int postProcessResolution = 0; //  keep low as possible
    
    int dim = gmesh->Dimension();
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(postProcessResolution,dim);
    an.PostProcess(postProcessResolution,dim);
    
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
}

// Teste: Acoplamento Darcy-Stokes

//Função Exata Darcy

// definition of f
void fD_source(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    f.resize(1);
    
    //STATE xv = x[0];
    //STATE yv = x[1];
    
    STATE f_x = 0.;
    
    f[0] = f_x;
}

// definition of v analytic
void vD_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    
    f.resize(2);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE v_x = -(-exp(-yv) + exp(yv))*cos(xv);
    STATE v_y = -(exp(-yv) + exp(yv))*sin(xv);
    
    f[0] = v_x; // x direction
    f[1] = v_y; // y direction
}

// definition of p analytic
void pD_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    
    f.resize(1);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE p_x = (-exp(-yv) + exp(yv))*sin(xv);
    
    f[0] = p_x;
}

// BC - Solução analítica - Artigo
void solucaoD_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    
    f.resize(3);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE v_x = -(-exp(-yv) + exp(yv))*cos(xv);
    STATE v_y = -(exp(-yv) + exp(yv))*sin(xv);
    STATE p = (-exp(-yv) + exp(yv))*sin(xv);
    
    f[0] = -v_x; // x direction
    f[1] = -v_y; // y direction
    f[2] = p; //
}

//Função Exata Stokes

// definition of f
void fS_source(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    
    f.resize(2);
    
    
    f[0] = 0.;
    f[1] = 0.;
    
    
    
    
}

// definition of v analytic
void vS_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    
    f.resize(2);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE v_x = (2.*cos(xv)*cos(Pi*yv)*sin(Pi*yv))/Pi;
    STATE v_y = sin(xv)*(-2. + pow(sin(Pi*yv),2)/(Pi*Pi));
    
    f[0] = v_x; // x direction
    f[1] = v_y; // y direction
}

// definition of p analytic
void pS_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    
    f.resize(1);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE p_x = sin(xv)*sin(yv);
    
    f[0] = p_x;
}

// source on interface
void grad_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    
    f.resize(1);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE gradv_x = 2.*cos(xv);
    //STATE gradv_y = 0.;
    
    f[0] = gradv_x; // x direction
    //f[1] = gradv_y; // y direction
}


// BC - Solução analítica - Artigo
void solucaoS_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    
    f.resize(3);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE v_x =  (2.*cos(xv)*cos(Pi*yv)*sin(Pi*yv))/Pi;
    STATE v_y =  sin(xv)*(-2. + pow(sin(Pi*yv),2)/(Pi*Pi));
    STATE p =  sin(xv)*sin(yv);
    
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
    
    //Stokes
    if(yv>0.){
        
        STATE v_x =  (2.*cos(xv)*cos(Pi*yv)*sin(Pi*yv))/Pi;
        STATE v_y =  sin(xv)*(-2. + pow(sin(Pi*yv),2)/(Pi*Pi));
        STATE p =  sin(xv)*sin(yv);
        
        sol[0]=v_x;
        sol[1]=v_y;
        sol[2]=p;
        
        // vx direction
        dsol(0,0)= -(2.*cos(Pi*yv)*sin(xv)*sin(Pi*yv))/Pi;
        dsol(0,1)= 2.*cos(xv)*cos(Pi*yv)*cos(Pi*yv) + 2.*cos(xv)*sin(Pi*yv)*sin(Pi*yv);
        
        // vy direction
        dsol(1,0)= cos(xv)*(-2. + sin(Pi*yv)*sin(Pi*yv)/(Pi*Pi));
        dsol(1,1)= (2.*cos(Pi*yv)*sin(xv)*sin(Pi*yv))/Pi;
        
        // Gradiente pressão
        
        dsol(2,0)= cos(xv)*sin(yv);
        dsol(2,1)= cos(yv)*sin(xv);
        
        
    }
    //Darcy
    if(yv<0.){
        
        STATE v_x = -(-exp(-yv) + exp(yv))*cos(xv);
        STATE v_y = -(exp(-yv) + exp(yv))*sin(xv);
        STATE p = (-exp(-yv) + exp(yv))*sin(xv);
        
        sol[0]=v_x;
        sol[1]=v_y;
        sol[2]=p;
        
        // vx direction
        dsol(0,0)= -(-exp(-yv) + exp(yv))*sin(xv);
        dsol(0,1)= (exp(-yv) + exp(yv))*cos(xv);
        
        // vy direction
        dsol(1,0)= (exp(-yv) + exp(yv))*cos(xv);
        dsol(1,1)= (-exp(-yv) + exp(yv))*sin(xv);
        
        // Gradiente pressão
        
        dsol(2,0)= v_x;
        dsol(2,1)= v_y;
        
        
    }
    
    
}



TPZGeoMesh *CreateGMesh(int nx, int ny, double hx, double hs, double hd)
{
    
    int i,j;
    long id, index, idno;
    
    
    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    
    //Vetor auxiliar para armazenar coordenadas:
    
    TPZVec <REAL> coord (3,0.);
    
    
    //Inicialização dos nós:
    
    long nnodes = nx*ny; //numero de nos do problema
    gmesh->NodeVec().Resize(nnodes); //Redimensiona o tamanho do vetor de nos da malha geometrica
    
    
    for(i = 0; i < ny; i++){
        for(j = 0; j < nx; j++){
            idno = i*nx + j;
            coord[0] = (j)*hx/(nx - 1);
            if(i<=ny/2){
                coord[1] = (i)*hd*2/(ny - 1);
            }else{
                coord[1] = hd+(i-(ny-1)/2)*hs*2/(ny - 1);
            }
            
            //using the same coordinate x for z
            coord[2] = 0.;
            //cout << coord << endl;
            //gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[idno].SetCoord(coord);
            //Get the index in the mesh nodes vector for the new node
            //Set the value of the node in the mesh nodes vector
            //gmesh->NodeVec()[idno] = TPZGeoNode(id,coord,*gmesh);
            gmesh->NodeVec()[idno].SetNodeId(idno);
        }
    }
    
    //    //    Ponto 1
    //    TPZVec<long> pointtopology(1);
    //    //pointtopology[0] = nx*(ny-1)/2;
    //    pointtopology[0] = 0;
    //
    //    gmesh->CreateGeoElement(EPoint,pointtopology,matPoint,id);
    //
    //    //    Ponto 2
    //    TPZVec<long> pointtopology2(1);
    //    pointtopology2[0] = nx*ny-1;
    //    //pointtopology2[0] = nx*((ny+1)/2);
    //
    //    gmesh->CreateGeoElement(EPoint,pointtopology2,matPoint2,id);
    //
    
    //Vetor auxiliar para armazenar as conecções entre elementos:
    
    TPZVec <long> connect(4,0);
    
    long int idD, idS;
    
    //Conectividade dos elementos (Parte Darcy):
    
    for(i = 0; i < (ny - 1)/2.; i++){
        for(j = 0; j < (nx - 1); j++){
            idD = (i)*(nx - 1)+ (j);
            connect[0] = (i)*nx + (j);
            connect[1] = connect[0]+1;
            connect[2] = connect[1]+(nx);
            connect[3] = connect[0]+(nx);
            gmesh->CreateGeoElement(EQuadrilateral,connect,matIdD,id);
        }
    }
    
    //Conectividade dos elementos (Parte Stokes):
    
    for(i = (ny - 1)/2.; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            idS = (i)*(nx - 1)+ (j);
            connect[0] = (i)*nx + (j);
            connect[1] = connect[0]+1;
            connect[2] = connect[1]+(nx);
            connect[3] = connect[0]+(nx);
            gmesh->CreateGeoElement(EQuadrilateral,connect,matIdS,id);
            
        }
    }
    
    
    //Gerando informação da vizinhança:
    
    gmesh->BuildConnectivity();
    
    {
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
    }
    long el, numelements = gmesh->NElements();
    
    TPZManVector <long> TopolPlate(4);
    
    for (el=0; el<numelements; el++)
    {
        long totalnodes = gmesh->ElementVec()[el]->NNodes();
        TPZGeoEl *plate = gmesh->ElementVec()[el];
        for (int i=0; i<4; i++){
            TopolPlate[i] = plate->NodeIndex(i);
        }
        
        //Colocando as condicoes de contorno:
        TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
        TPZManVector <REAL,3> nodecoord(3);
        
        //Na face x = 1
        TPZVec<long> ncoordzbottVec(0); long sizeOfbottVec = 0;
        TPZVec<long> ncoordzbottVec2(0); long sizeOfbottVec2 = 0;
        
        
        TPZVec<long> ncoordztopVec(0); long sizeOftopVec = 0;
        TPZVec<long> ncoordzleftupVec(0); long sizeOfleftupVec = 0;
        TPZVec<long> ncoordzleftdownVec(0); long sizeOfleftdownVec = 0;
        TPZVec<long> ncoordzrightupVec(0); long sizeOfrightupVec = 0;
        TPZVec<long> ncoordzrightdownVec(0); long sizeOfrightdownVec = 0;
        
        for (long i = 0; i < totalnodes; i++)
        {
            Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (nodecoord[2] == 0. & nodecoord[1] == 0.)
            {
                sizeOfbottVec++;
                ncoordzbottVec.Resize(sizeOfbottVec);
                ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
            }
            
            //            //teste!!
            //            if (nodecoord[2] == 0. & nodecoord[1] == 0.)
            //            {
            //                sizeOfbottVec2++;
            //                ncoordzbottVec2.Resize(sizeOfbottVec2);
            //                ncoordzbottVec2[sizeOfbottVec2-1] = TopolPlate[i];
            //            }
            
            if (nodecoord[2] == 0. & nodecoord[1] == hs+hd)
            {
                sizeOftopVec++;
                ncoordztopVec.Resize(sizeOftopVec);
                ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == 0.& nodecoord[1] <= hd)
            {
                sizeOfleftdownVec++;
                ncoordzleftdownVec.Resize(sizeOfleftdownVec);
                ncoordzleftdownVec[sizeOfleftdownVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == 0.& nodecoord[1] >= hd)
            {
                sizeOfleftupVec++;
                ncoordzleftupVec.Resize(sizeOfleftupVec);
                ncoordzleftupVec[sizeOfleftupVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == hx & nodecoord[1] <= hd)
            {
                sizeOfrightdownVec++;
                ncoordzrightdownVec.Resize(sizeOfrightdownVec);
                ncoordzrightdownVec[sizeOfrightdownVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == hx & nodecoord[1] >= hd)
            {
                sizeOfrightupVec++;
                ncoordzrightupVec.Resize(sizeOfrightupVec);
                ncoordzrightupVec[sizeOfrightupVec-1] = TopolPlate[i];
            }
        }
        
        if (sizeOfbottVec == 2) {
            int sidesbott = plate->WhichSide(ncoordzbottVec);
            TPZGeoElSide platesidebott(plate, sidesbott);
            TPZGeoElBC(platesidebott,matDBCbott);
            TPZGeoElBC(platesidebott,matIntDBCbott);
        }
        
        //        if (sizeOfbottVec2 == 2) {
        //            int sidesbott2 = plate->WhichSide(ncoordzbottVec2);
        //            TPZGeoElSide platesidebott2(plate, sidesbott2);
        //            TPZGeoElBC(platesidebott2,matSBCbott);
        //            TPZGeoElBC(platesidebott2,matIntSBCbott);
        //        }
        
        
        if (sizeOftopVec == 2) {
            int sidestop = plate->WhichSide(ncoordztopVec);
            TPZGeoElSide platesidetop(plate, sidestop);
            TPZGeoElBC(platesidetop,matSBCtop);
            TPZGeoElBC(platesidetop,matIntSBCtop);
        }
        
        if (sizeOfleftdownVec == 2) {
            int sidesleftdown = plate->WhichSide(ncoordzleftdownVec);
            TPZGeoElSide platesideleftdown(plate, sidesleftdown);
            TPZGeoElBC(platesideleftdown,matDBCleft);
            TPZGeoElBC(platesideleftdown,matIntDBCleft);
        }
        
        if (sizeOfleftupVec == 2) {
            int sidesleftup = plate->WhichSide(ncoordzleftupVec);
            TPZGeoElSide platesideleftup(plate, sidesleftup);
            TPZGeoElBC(platesideleftup,matSBCleft);
            TPZGeoElBC(platesideleftup,matIntSBCleft);
        }
        
        if (sizeOfrightdownVec == 2) {
            int sidesrightdown = plate->WhichSide(ncoordzrightdownVec);
            TPZGeoElSide platesiderightdown(plate, sidesrightdown);
            TPZGeoElBC(platesiderightdown,matDBCright);
            TPZGeoElBC(platesiderightdown,matIntDBCright);
        }
        
        if (sizeOfrightupVec == 2) {
            int sidesrightup = plate->WhichSide(ncoordzrightupVec);
            TPZGeoElSide platesiderightup(plate, sidesrightup);
            TPZGeoElBC(platesiderightup,matSBCright);
            TPZGeoElBC(platesiderightup,matIntSBCright);
        }
        
        ncoordzbottVec.Resize(0);
        sizeOfbottVec = 0;
        
        //        ncoordzbottVec2.Resize(0);
        //        sizeOfbottVec2 = 0;
        
        ncoordztopVec.Resize(0);
        sizeOftopVec = 0;
        ncoordzleftdownVec.Resize(0);
        sizeOfleftdownVec = 0;
        ncoordzleftupVec.Resize(0);
        sizeOfleftupVec = 0;
        ncoordzrightdownVec.Resize(0);
        sizeOfrightdownVec = 0;
        ncoordzrightupVec.Resize(0);
        sizeOfrightupVec = 0;
        
    }
    
    // Criando e inserindo elemento de interfação:
    //    TPZVec<long> nodind3(2);
    //
    //    nodind3[0]=1;
    //    nodind3[1]=4;
    //
    //    gmesh->CreateGeoElement(EOned, nodind3, matInterface, index); //Criando elemento de interface (GeoElement)
    
    
    //Criando interface (Geralizado):
    
    TPZVec<long> nodint(2);
    for(i = 0; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            
            if(j>0&&j<(nx-1)){
                nodint[0]=j+nx*i;
                nodint[1]=j+nx*(i+1);
                if(i<(ny-1)/2.){
                    gmesh->CreateGeoElement(EOned, nodint, matDInterface, id); //Criando elemento de interface (GeoElement)
                }
                if(i>=(ny-1)/2.){
                    gmesh->CreateGeoElement(EOned, nodint, matSInterface, id); //Criando elemento de interface (GeoElement)
                }
                
            }
            if(i>0&&i<(ny-1)){
                nodint[0]=j+nx*i;
                nodint[1]=j+nx*i+1;
                
                if(i==(ny-1)/2.){
                    gmesh->CreateGeoElement(EOned, nodint, matInterfaceDS, id); //Criando elemento de interface (GeoElement)
                }
                if(i<(ny-1)/2.){
                    
                    gmesh->CreateGeoElement(EOned, nodint, matDInterface, id); //Criando elemento de interface (GeoElement)
                    
                }
                if(i>(ny-1)/2.){
                    
                    gmesh->CreateGeoElement(EOned, nodint, matSInterface, id); //Criando elemento de interface (GeoElement)
                }
                
                
            }
        }
    }
    
    
    //new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (nodind3,matInterface,*gmesh); //Criando elemento de interface (RefPattern)
    //id++;
    
    gmesh->AddInterfaceMaterial(matIdS, matIdD, matInterfaceDS);
    gmesh->AddInterfaceMaterial(matIdD, matIdS, matInterfaceDS);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();
    
    //Impressão da malha geométrica:
    
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;
    
    
    
}

TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}


TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim);//Insere dimensão do modelo
    
    
    //Definição do espaço de aprximação:
    
    //cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1:
    
    cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
    
    
    //Criando elementos com graus de liberdade differentes para cada elemento (descontínuo):
    
    //cmesh->ApproxSpace().CreateDisconnectedElements(true); //Criando elementos desconectados (descontínuo)
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Darcy:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *materialDarcy = new TPZMat2dLin(matIdD); //Criando material que implementa a formulação fraca do problema modelo
    
    cmesh->InsertMaterialObject(materialDarcy); //Insere material na malha
    
    //Dimensões do material (para H1 e descontinuo):
    //TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
    //materialDarcy->SetMaterial(xkin, xcin, xfin);
    
    //Dimensões do material (para HDiv):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    materialDarcy->SetMaterial(xkin, xcin, xfin);
    
    
    //Condições de contorno:
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial * BCDond0 = materialDarcy->CreateBC(materialDarcy, matDBCbott, neumann, val1, val2); //Cria material que implementa a condição de contorno inferior
    cmesh->InsertMaterialObject(BCDond0); //Insere material na malha
    
    // TPZMaterial * BCDond1 = materialDarcy->CreateBC(materialDarcy, matDBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    // cmesh->InsertMaterialObject(BCDond1); //Insere material na malha
    
    TPZMaterial * BCDond2 = materialDarcy->CreateBC(materialDarcy, matDBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCDond2); //Insere material na malha
    
    TPZMaterial * BCDond3 = materialDarcy->CreateBC(materialDarcy, matDBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCDond3); //Insere material na malha
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Stokes:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *materialStokes = new TPZMat2dLin(matIdS); //Criando material que implementa a formulação fraca do problema modelo
    
    cmesh->InsertMaterialObject(materialStokes); //Insere material na malha
    
    //Dimensões do material (para H1 e descontinuo):
    //TPZFMatrix<STATE> xkin2(2,2,0.), xcin2(2,2,0.), xfin2(2,2,0.);
    //materialStokes->SetMaterial(xkin2, xcin2, xfin2);
    
    //Dimensões do material (para HDiv):
    TPZFMatrix<STATE> xkin2(1,1,0.), xcin2(1,1,0.), xfin2(1,1,0.);
    materialStokes->SetMaterial(xkin2, xcin2, xfin2);
    
    
    //Condições de contorno:
    
    //TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
    //    TPZMaterial * BCSond0 = materialStokes->CreateBC(materialStokes, matSBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //    cmesh->InsertMaterialObject(BCSond0); //Insere material na malha
    
    TPZMaterial * BCSond1 = materialStokes->CreateBC(materialStokes, matSBCtop, penetration, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCSond1); //Insere material na malha
    
    TPZMaterial * BCSond2 = materialStokes->CreateBC(materialStokes, matSBCleft, penetration, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCSond2); //Insere material na malha
    
    TPZMaterial * BCSond3 = materialStokes->CreateBC(materialStokes, matSBCright, penetration, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCSond3); //Insere material na malha
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Coupling:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat1dLin *materialCoupling = new TPZMat1dLin(matInterfaceDS); //Criando material que implementa a formulação fraca do problema modelo
    
    cmesh->InsertMaterialObject(materialCoupling); //Insere material na malha
    
    //Dimensões do material (para H1 e descontinuo):
    //TPZFMatrix<STATE> xkin3(2,2,0.), xcin3(2,2,0.),xbin3(2,2,0.), xfin3(2,2,0.);
    //materialCoupling->SetMaterial(xkin3, xcin3, xbin3, xfin3);
    
    //Dimensões do material (para HDiv):
    TPZFMatrix<STATE> xkin3(1,1,0.), xcin3(1,1,0.), xbin3(1,1,0.), xfin3(1,1,0.);
    materialCoupling->SetMaterial(xkin3, xcin3, xbin3, xfin3);
    
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

TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder)
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
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Darcy:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *materialDarcy = new TPZMat2dLin(matIdD);//criando material que implementa a formulacao fraca do problema modelo
    
    cmesh->InsertMaterialObject(materialDarcy); //Insere material na malha
    
    //Dimensões do material (para H1 e descontínuo):
    TPZFMatrix<STATE> xkin1(1,1,0.), xcin1(1,1,0.), xfin1(1,1,0.);
    materialDarcy->SetMaterial(xkin1, xcin1, xfin1);
    
    //Condições de contorno:
    
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
    //    val2(0,0) = 0.0; // px -> 0
    //    val2(1,0) = 0.0; // py -> 0
    //
    //    TPZMaterial * BCond0 = material->CreateBC(material, matDBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    //
    //    TPZMaterial * BCond1 = material->CreateBC(material, matDBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    //    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    //
    //    TPZMaterial * BCond2 = material->CreateBC(material, matDBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    //
    //    TPZMaterial * BCond3 = material->CreateBC(material, matDBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    //    Ponto de pressao:
    //
    //    TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
    //
    //    TPZMaterial * BCPointD = materialDarcy->CreateBC(materialDarcy, matPoint, pointtype, val3, val4); //Cria material que implementa um ponto para a pressao
    //    cmesh->InsertMaterialObject(BCPointD); //Insere material na malha
    //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Stokes:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *materialStokes = new TPZMat2dLin(matIdS);//criando material que implementa a formulacao fraca do problema modelo
    
    cmesh->InsertMaterialObject(materialStokes); //Insere material na malha
    
    //Dimensões do material (para H1 e descontínuo):
    TPZFMatrix<STATE> xkin2(1,1,0.), xcin2(1,1,0.), xfin2(1,1,0.);
    materialStokes->SetMaterial(xkin2, xcin2, xfin2);
    
    //Condições de contorno:
    
    
    //TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
    //    val2(0,0) = 0.0; // px -> 0
    //    val2(1,0) = 0.0; // py -> 0
    //
    //    TPZMaterial * BCond0 = material->CreateBC(material, matSBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    //
    //    TPZMaterial * BCond1 = material->CreateBC(material, matSBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    //    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    //
    //    TPZMaterial * BCond2 = material->CreateBC(material, matSBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    //
    //    TPZMaterial * BCond3 = material->CreateBC(material, matSBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    //Ponto de pressao:
    //
    //    TPZFMatrix<REAL> val5(1,1,0.), val6(1,1,0.);
    //
    //    TPZMaterial * BCPointS = materialStokes->CreateBC(materialStokes, matPoint2, pointtype, val5, val6); //Cria material que implementa um ponto para a pressao
    //    cmesh->InsertMaterialObject(BCPointS); //Insere material na malha
    
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Coupling:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat1dLin *materialCoupling = new TPZMat1dLin(matInterfaceDS); //Criando material que implementa a formulação fraca do problema modelo
    
    cmesh->InsertMaterialObject(materialCoupling); //Insere material na malha
    
    //Dimensões do material (para H1 e descontinuo):
    TPZFMatrix<STATE> xkin3(1,1,0.), xcin3(1,1,0.), xbin3(1,1,0.), xfin3(1,1,0.);
    materialCoupling->SetMaterial(xkin3, xcin3, xbin3, xfin3);
    
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    std::set<int> materialids;
    materialids.insert(matIdD);
    materialids.insert(matIdS);
    materialids.insert(matInterfaceDS);
    //materialids.insert(matInterfaceDS);
    cmesh->AutoBuild(materialids);
    cmesh->LoadReferences();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    std::set<int> materialids2;
    materialids2.insert(matPoint);
    //materialids2.insert(matpoint2);
    cmesh->AutoBuild(materialids2);
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    
    //   cmesh->AdjustBoundaryElements();
    //   cmesh->CleanUpUnconnectedNodes();
    
    
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
    
    // Criando material Darcy:
    
    TPZDarcyPMaterial *materialDarcy = new TPZDarcyPMaterial(matIdD,dim,permeability,theta);//criando material que implementa a formulacao fraca do problema modelo
    // Inserindo material na malha
    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (fD_source);
    TPZAutoPointer<TPZFunction<STATE> > pp = new TPZDummyFunction<STATE> (pD_exact);
    TPZAutoPointer<TPZFunction<STATE> > vp = new TPZDummyFunction<STATE> (vD_exact);
    
    materialDarcy->SetForcingFunction(fp);
  //  materialDarcy->SetForcingFunctionExactPressure(pp);
    materialDarcy->SetForcingFunctionExact(vp);
    cmesh->InsertMaterialObject(materialDarcy);
    
    
    //Condições de contorno:
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0
    
    TPZFMatrix<REAL> val2tp(2,1,0.);
    
    val2tp(0,0) = 0.0; // vx -> 0
    val2tp(1,0) = -0.015; // vy -> 0
    
    TPZMaterial * BCDond0 = materialDarcy->CreateBC(materialDarcy, matDBCbott, neumann, val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    //BCDond0->SetForcingFunction(solucaoD_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCDond0); //Insere material na malha
    
    //    TPZMaterial * BCDond1 = materialDarcy->CreateBC(materialDarcy, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    //    BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //    BCDond1->SetForcingFunction(solucaoD_exact,bc_inte_order);
    //    cmesh->InsertMaterialObject(BCDond1); //Insere material na malha
    
    TPZMaterial * BCDond2 = materialDarcy->CreateBC(materialDarcy, matDBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    //BCDond2->SetForcingFunction(solucaoD_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCDond2); //Insere material na malha
    
    TPZMaterial * BCDond3 = materialDarcy->CreateBC(materialDarcy, matDBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    //BCDond3->SetForcingFunction(solucaoD_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCDond3); //Insere material na malha
    
    //Ponto
    
    //    TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
    //    val4(0,0)=0.0;
    //
    //    TPZMaterial * BCPointD = materialDarcy->CreateBC(materialDarcy, matPoint, pointtype, val3, val4); //Cria material que implementa um ponto para a pressão
    //    cmesh->InsertMaterialObject(BCPointD); //Insere material na malha
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Stokes:
    
    TPZStokesMaterial *materialStokes = new TPZStokesMaterial(matIdS,dim,visco,theta);//criando material que implementa a formulacao fraca do problema modelo
    // Inserindo material na malha
    TPZAutoPointer<TPZFunction<STATE> > fp2 = new TPZDummyFunction<STATE> (fS_source);
    TPZAutoPointer<TPZFunction<STATE> > pp2 = new TPZDummyFunction<STATE> (pS_exact);
    TPZAutoPointer<TPZFunction<STATE> > vp2 = new TPZDummyFunction<STATE> (vS_exact);
    
    materialStokes->SetForcingFunction(fp2);
   // materialStokes->SetForcingFunctionExactPressure(pp2);
    materialStokes->SetForcingFunctionExact(vp2);
    cmesh->InsertMaterialObject(materialStokes);
    
    
    //Condições de contorno:
    
    //    TPZMaterial * BCSond0 = materialStokes->CreateBC(materialStokes, matSBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    //    BCSond0->SetForcingFunction(solucaoS_exact,bc_inte_order);
    //    cmesh->InsertMaterialObject(BCSond0); //Insere material na malha
    
    TPZMaterial * BCSond1 = materialStokes->CreateBC(materialStokes, matSBCtop, penetration, val1, val2tp); //Cria material que implementa a condicao de contorno superior
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //BCSond1->SetForcingFunction(solucaoS_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCSond1); //Insere material na malha
    
    TPZMaterial * BCSond2 = materialStokes->CreateBC(materialStokes, matSBCleft, penetration, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    //BCSond2->SetForcingFunction(solucaoS_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCSond2); //Insere material na malha
    
    TPZMaterial * BCSond3 = materialStokes->CreateBC(materialStokes, matSBCright, penetration, val1, val2); //Cria material que implementa a condicao de contorno direita
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    //BCSond3->SetForcingFunction(solucaoS_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCSond3); //Insere material na malha
    
    //Ponto
    
    //    TPZFMatrix<REAL> val5(1,1,0.), val6(1,1,0.);
    //    val6(0,0)=0.0;
    //
    //    TPZMaterial * BCPointS = materialStokes->CreateBC(materialStokes, matPoint2, pointtype, val3, val4); //Cria material que implementa um ponto para a pressão
    //    cmesh->InsertMaterialObject(BCPointS); //Insere material na malha
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Acoplamento:
    
    TPZCouplingDSMaterial *materialCoupling = new TPZCouplingDSMaterial(matInterfaceDS,dim,visco,permeability,theta);
    cmesh->InsertMaterialObject(materialCoupling);
    
    TPZAutoPointer<TPZFunction<STATE> > fp3 = new TPZDummyFunction<STATE> (grad_exact);
    TPZAutoPointer<TPZFunction<STATE> > pp3 = new TPZDummyFunction<STATE> (pS_exact);
    TPZAutoPointer<TPZFunction<STATE> > vp3 = new TPZDummyFunction<STATE> (vS_exact);
    
    materialCoupling->SetForcingFunction(fp3);
   // materialCoupling->SetForcingFunctionExactPressure(pp3);
    materialCoupling->SetForcingFunctionExact(vp3);
    cmesh->InsertMaterialObject(materialCoupling);
    
    
    
#ifdef PZDEBUG
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}


void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    long nel = gmesh->NElements();
    for (long el = 0; el<nel; el++) {
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
        long index;
        new TPZMultiphysicsInterfaceElement(cmesh,gel,index,celstack[1],celstack[0]);
    }
    
}

void AddInterfaceCoupllingDS(TPZGeoMesh *gmesh,TPZCompMesh *cmesh, int matInterfaceDS ,int matleft, int matright){
    
    //gmesh->ResetReference();
    //cmesh->LoadReferences();
    
    const int nel = gmesh->NElements();
    long index;
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        const int gelMatId = gel->MaterialId();
        if (gelMatId != matInterfaceDS)
            continue;
        if (gel->HasSubElement())
            continue;
        TPZCompElSide left, right;
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZGeoElSide neigh = gelside.Neighbour();
        
        int debug = 0;
        while (gelside != neigh) {
            if (neigh.Element()->MaterialId() == matleft) {
                TPZCompEl *cel = neigh.Element()->Reference();
                if (!cel) {
                    DebugStop();
                }
                left = TPZCompElSide(cel, 6);
                TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace*>(cel);
                debug++;
            }
            else if (neigh.Element()->MaterialId() == matright) {
                TPZCompEl *cel = neigh.Element()->Reference();
                if (!cel) {
                    DebugStop();
                }
                right = TPZCompElSide(cel, 4);
                debug++;
            }
            neigh = neigh.Neighbour();
        }
        if (debug != 2)
            DebugStop();
        new TPZMultiphysicsInterfaceElement(*cmesh, gel, index, left, right);
    }
    
}


//Função Erro (Não utilizada)
void Error(TPZCompMesh *cmesh, std::ostream &out, int p, int ndiv)
{
    DebugStop();
    long nel = cmesh->NElements();
    //int dim = cmesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        TPZManVector<STATE,10> elerror(10,0.);
        cel->EvaluateError(sol_exact, elerror, NULL);
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


