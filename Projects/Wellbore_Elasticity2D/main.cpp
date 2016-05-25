#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzmat2dlin.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzfunction.h"
#include "TPZVTKGeoMesh.h"

#include <time.h>
#include <stdio.h>
#include <fstream>

#include "pzlog.h"
#include "pzgmesh.h"
#include "TPZMatElasticity2D.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"

#include "pzstepsolver.h"
#include "math.h"
#include "TPZSkylineNSymStructMatrix.h"

#include <cmath>
#include <set>

#include <pzgengrid.h>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.Wellbore_Elasticity2D"));
static LoggerPtr loggeradap(Logger::getLogger("pz.adaptivity"));
static LoggerPtr loggerconv(Logger::getLogger("pz.adaptivity.conv"));
static LoggerPtr loggerpoint(Logger::getLogger("pz.adaptivity.points"));
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>

using namespace std;


TPZGeoMesh *GetMesh (REAL rwb, REAL re, int ncirc, int nrad, REAL DrDcirc);
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder);


int main(int argc, char *argv[])
{
    
    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    std::string FileName = dirname;
    FileName = dirname + "/Projects/Wellbore_Elasticity2D/";
    FileName += "Wellbore_Elasticity2DLog.cfg";
    InitializePZLOG(FileName);
#endif
    
    
    
    //******** Configura malha geometrica ***************/
    // rw = raio do poco (metros)
    // rext = raio externo do contorno (metros)
    // ncircle = nro elementos em 1/4 da parede do poco
    // nradial = nro de elementos da parede do poco ate o raio externo
    // drdcirc = proporcao do primeiro elemento
    REAL rw = 0.1;
    REAL rext = 1.0;
    int ncircle = 30;
    int nradial = 30;
    REAL drdcirc = 1.5;
    
    TPZGeoMesh *gmesh = GetMesh(rw, rext, ncircle, nradial, drdcirc); //funcao para criar a malha geometrica
    const std::string nm("line");
    gmesh->SetName(nm);
    std::ofstream outtxt("gmesh.txt"); //define arquivo de saida para impressao dos dados da malha
    gmesh->Print(outtxt);
    std::ofstream out("gmesh.vtk"); //define arquivo de saida para impressao da malha no paraview
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true); //imprime a malha no formato vtk
    
   
    
    //******** Configura malha Computacional ***************/
    
    int p = 3;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CMesh(gmesh, p);
    
    // Resolvendo o Sistema
    bool optimizeBandwidth = false; //impede a renumeracao das equacoes do problema(para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh, optimizeBandwidth); //cria objeto de analise que gerenciaria a analise do problema
    
    
//    // Solving linear equations
//    // Initial steps
//    TPZAnalysis an (cmesh);
//    TPZSkylineStructMatrix strskyl(cmesh);
//    an.SetStructuralMatrix(strskyl);
//    // Solver (is your choose)
//    TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
//    direct->SetDirect(ECholesky);
//    an.SetSolver(*direct);
//    delete direct;
//    direct = 0;
//    
//    an.Run();
    
    
    //************ Para visualizar K e Rhs ******************************//
    int numofThreads = 0;
    TPZSkylineNSymStructMatrix skylnsym(cmesh);
    TPZStepSolver<STATE> step;
    skylnsym.SetNumThreads(numofThreads);
    step.SetDirect(ELU);
    an.SetStructuralMatrix(skylnsym);
    an.SetSolver(step);
    
    an.Assemble();
    an.Rhs() ;
    
    TPZAutoPointer< TPZMatrix<REAL> > KGlobal;
    TPZFMatrix<STATE> FGlobal;
    KGlobal =   an.Solver().Matrix();
    FGlobal =   an.Rhs();
    
#ifdef PZDEBUG
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        KGlobal->Print("k = ", sout,EMathematicaInput);
        FGlobal.Print("r = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
#endif
    
    
    an.Solve();//assembla a matriz de rigidez (e o vetor de carga) global e inverte o sistema de equacoes
    
    TPZFMatrix<REAL> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
    solucao.Print("Sol",cout,EMathematicaInput);//imprime na formatacao do Mathematica

    
    
    // Post processing
    TPZManVector<std::string> scalarnames(3), vecnames(1);
    scalarnames[0] = "SigmaX";
    scalarnames[1] = "SigmaY";
    scalarnames[2] = "Pressure";
    vecnames[0] = "displacement";
    //vecnames[1] = "";
    an.DefineGraphMesh(2,scalarnames,vecnames,"ElasticitySolutions.vtk");
    
    an.PostProcess(6);
    

    
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
    
    
#ifdef PZDEBUG
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        
        std::stringstream sout;
        sout << " Geometry out ... " << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
#endif
    
    

    

}





// *********** Cria malha Geometrica *************/

// Para gerar a malha eh necessario:
// rwb -> raio do poco
// re -> raio do contorno
// ncirc -> nro elem ao longo 1/4 do poco
// nrad -> nro elem da parede do poco ate contorno externo
// DrDcirc -> proporcao dos elementos da parede do poco

TPZGeoMesh *GetMesh (REAL rwb, REAL re, int ncirc, int nrad, REAL DrDcirc) {
    

    // calcula comprimento radial do primeiro elemento
    REAL szmin;
    REAL Pi = 3.14159;
    szmin = (Pi/2)*(rwb/ncirc)*(DrDcirc);
    
    // calcula comprimento radial da parede do poco ate contorno
    REAL radiallenght;
    radiallenght = re - rwb;
    
    // definindo variacao do angulo theta ao redor do poco
    // em rads!!!!!
    TPZVec<REAL> theta;
    theta.Resize(ncirc+1);
    REAL firsttheta = (Pi/2) / (ncirc);
    for (int k = 0; k<ncirc+1; k++) {
        REAL sumtheta = 0.;
        sumtheta += firsttheta * k;
        theta[k] = sumtheta;
    }

    
    // *******Imprime variacao dos angulos (em rads)
    std::cout<< "elementos de theta: " << endl;
        for (int t=0; t<ncirc+1; t++) {
            std::cout<< "Theta[" << t << "] :" << theta[t] << endl;
            }
            std::cout<< "Print theta " << endl;
            theta.Print();
            std::cout << endl;
    
    
    // nx = number of nodes in x direction
    // ny = number of nodes in y direction
    int nx,ny;
    nx = nrad+1;
    ny = ncirc+1;
    
   
    // Geometric Progression of the elements
    REAL q;
    q = TPZGenGrid::GeometricProgression(szmin, radiallenght, nrad);
    std::cout<< "valor de q " << q << endl; // imprime razao da PG
    
    //Creates the geometric mesh... The nodes and elements
    //will be inserted into mesh object during initilize process
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    long i,j;
    long id, index;
    
    
//*********************** Malha Circunferencial (1/4) *********************//
    
    //vector to store a coordinate
    TPZVec <REAL> coord (2,0.);
    
    // aloca valor acumulado dos raios
    REAL rsum = 0.;
    
    //Nodes initialization
    for(i = 1; i < nx+1; i++){
        for(j = 1; j < ny+1; j++){
          
            // aloca coordenadas em cada no
            coord[0] = (rwb + rsum)* cos(theta[j-1]);
            coord[1] = (rwb + rsum)* sin(theta[j-1]);
            
            // id do elemento
            id = (i) * ny + (j);
            
            //Get the index in the mesh nodes vector for the new node
            index = gmesh->NodeVec().AllocateNewElement();
            //Set the value of the node in the mesh nodes vector
            gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
            
            
            // print
            
            std::cout << "*****Iteracao nro: " << j << endl;
            std::cout << "rsum: " << rsum << endl;
            std::cout << "cos" << "[" << theta[j-1] << "]" <<": " ;
            std::cout << cos(theta[j-1]);
            std::cout << endl;
            std::cout << "sin" << "[" << theta[j-1] << "]" << ": ";
            std::cout << sin(theta[j-1]);
            std::cout << endl;
            std::cout << "Coord x: " << coord[0] << ";" << " Coord y: " << coord[1] << ";" << " Rad: " << theta[j-1] << endl;
            std::cout << endl;
            
            
        }
                    rsum += (szmin * (pow(q,i-1))); //valor acumulado dos raios
    }
    
    
    
    //vector to store a element connectivities
    TPZVec <long> connect(4,0);
    
    
    //Element connectivities
    for(i = 0; i < (nx -1); i++){
        for(j = 0; j < (ny -1); j++){
           // index = (i)*(ny - 1)+ (j);
            connect[0] = (i) * ny + (j);
            connect[1] = connect[0]+(ny);
            connect[2] = connect[1]+1;
            connect[3] = connect[0]+1;
            //Allocates and define the geometric element
            gmesh->CreateGeoElement(EQuadrilateral,connect,1,id);
            std::cout << "connect: " << connect << endl;

            gmesh->ElementVec()[id];
        }
    }
    //Generate neighborhod information
    gmesh->BuildConnectivity();
    
    
    
    //********** Criando Geo de BC, verificar como criar GeoElBC para cada elemento e definir seu lado que sera atribuido a BC ********//
    
    // bc = -1 -> Normal Pressure condition
    for (int i = 0; i<ncirc; i++ ) {
        gmesh->ElementVec()[i]->CreateBCGeoEl(7, -1);
       //TPZGeoElBC gbc1(&elvec[i],7,-1);
    }

    // bc = -2 -> Neumann condition contorno externo bottom
    for (int i = 0; i<nrad; i++ ) {
        gmesh->ElementVec()[ncirc*i]->CreateBCGeoEl(4, -2);
        //TPZGeoElBC gbc2(&elvec[ncirc*i],4,-2);
    }
    
    // bc = -3 -> Neumann condition contorno externo upper
    for (int i = 1; i<nrad+1; i++ ) {
        gmesh->ElementVec()[(ncirc*i)-1]->CreateBCGeoEl(6, -3);
        //TPZGeoElBC gbc3(&elvec[(ncirc*i)-1],6,-3);
    }
    
    // bc = -4 -> Neumann condition arco externo do farfield
    for (int i = 1; i<ncirc+1; i++ ) {
        gmesh->ElementVec()[(ncirc*nrad)-i]->CreateBCGeoEl(5, -4);
        //TPZGeoElBC gbc4(&elvec[(ncirc*nrad)-i],5,-4);
    }
    
    
    return gmesh;
    
}

    



// Cria malha COmputacional

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder)
{
    int matId = 1;
    const int dim = 2; //dimensao do problema


    // Plane strain assumption
   // int linestrain = 1;
    
    //**************** Criando material  ********************************
    TPZMatElasticity2D *material = new TPZMatElasticity2D(matId);//criando material que implementa a formulacao fraca do problema modelo
    
    
    // Setting up paremeters
    REAL Eyoung = 1.0 , ni = 2.0, fbx = 0., fby = 0.;
    material->SetElasticity(Eyoung, ni, fbx, fby);
    
    
    /******* Calculating Inicial Stresses *******/
    // direction = direction/azimuth
    // inclination = wellbore inclination
    // problem assumption, inclined wellbore state = 1
    double direction = 0., inclination = 0.; //rad
    int inclinedwellbore = 1;
    
    // Tensoes in Situ, horizontais e vertical
    REAL SigmaVV = 0., Sigmahh = 0., SigmaHH = 0.;
    
  
    // Seta os parametros do poco 
    material->SetInclinedWellboreParameters(SigmaHH, Sigmahh, SigmaVV, direction, inclination, inclinedwellbore);
    
    
    //Obtem tensor de tensoes iniciais
    REAL SigmaX = 0., SigmaXY = 0., SigmaY = 0., SigmaZ = 0.;
    material->SetPreStress(SigmaX, SigmaXY, SigmaY, SigmaZ);
    
    
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmesh->SetDimModel(dim);//seta dimensao do modelo
    
    // Inserindo material na malha
    cmesh->InsertMaterialObject(material);
    
    
    // cond contorno
    const int bc0 = -1, bc1 = -2, bc2 = -3, bc3 = -4; // ids igual da malha geometrica
    const int neumann = 1, normalpressure = 6, stressfield = 4; // tipo de condicao de contorno
    
    ///Inserir condicao de contorno parede do poco
    TPZFMatrix<REAL> val1(2,1,0.), val2(2,1,0.);
    val1(0,0) = 0.0;
    val1(1,0) = 0.0;
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond0 = material->CreateBC(material, bc0, normalpressure, val1, val2);//cria material que implementa a condicao de contorno da parede do poco
    
    //Condicao de contorno externo bottom
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond1 = material->CreateBC(material, bc1, neumann, val1, val2);//cria material que implementa a condicao de contorno da parede do poco
    
    //Condicao de contorno externo bottom
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond2 = material->CreateBC(material, bc2, neumann, val1, val2);//cria material que implementa a condicao de contorno da parede do poco
    
    //Condicao de contorno arco externo do farfield
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond3 = material->CreateBC(material, bc3, stressfield, val1, val2);//cria material que implementa a condicao de contorno da parede do poco
   
   
    cmesh->InsertMaterialObject(BCond0);//insere material na malha
    cmesh->InsertMaterialObject(BCond1);//insere material na malha
    cmesh->InsertMaterialObject(BCond2);//insere material na malha
    cmesh->InsertMaterialObject(BCond3);//insere material na malha
    
    cmesh->SetAllCreateFunctionsContinuous();
    
    //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    
   return cmesh;
    
}




