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
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
//static LoggerPtr loggeradap(Logger::getLogger("pz.adaptivity"));
//static LoggerPtr loggerconv(Logger::getLogger("pz.adaptivity.conv"));
//static LoggerPtr loggerpoint(Logger::getLogger("pz.adaptivity.points"));
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>

using namespace std;


TPZGeoMesh *GetMesh (REAL rwb, REAL re, int ncirc, int nrad, REAL DrDcirc);
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder);
TPZGeoMesh *CircularGeoMesh (REAL rwb, REAL re, int ncirc, int nrad, REAL DrDcirc);
TPZCompMesh *CircularCMesh(TPZGeoMesh *gmesh, int pOrder);


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
    REAL rext = 0.5;
    int ncircle = 4;
    int nradial = 4;
    REAL drdcirc = 1.5;
    
    
    TPZGeoMesh *gmesh = CircularGeoMesh (rw, rext, ncircle, nradial, drdcirc); //funcao para criar a malha GEOMETRICA de todo o poco
    //TPZGeoMesh *gmesh = GetMesh(rw, rext, ncircle, nradial, drdcirc); //funcao para criar a malha GEOMETRICA de 1/4 do poco
    const std::string nm("line");
    gmesh->SetName(nm);
    std::ofstream outtxt("gmesh.txt"); //define arquivo de saida para impressao dos dados da malha
    gmesh->Print(outtxt);
    std::ofstream out("gmesh.vtk"); //define arquivo de saida para impressao da malha no paraview
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true); //imprime a malha no formato vtk
    

    
    //******** Configura malha Computacional ***************/
    
    int p = 1;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CircularCMesh(gmesh, p); //funcao para criar a malha COMPUTACIONAL de todo o poco
    //TPZCompMesh *cmesh = CMesh(gmesh, p); //funcao para criar a malha COMPUTACIONAL de 1/4 do poco

    // Solving linear equations
    // Initial steps
    TPZAnalysis an (cmesh);
    TPZSkylineStructMatrix strskyl(cmesh);
    an.SetStructuralMatrix(strskyl);
    // Solver (is your choose)
    TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
    direct->SetDirect(ECholesky);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    
    an.Run();
//
//    
////    //************ Para visualizar K e Rhs ******************************//
////    int numofThreads = 0;
////    TPZSkylineNSymStructMatrix skylnsym(cmesh);
////    TPZStepSolver<STATE> step;
////    skylnsym.SetNumThreads(numofThreads);
////    step.SetDirect(ELU);
////    an.SetStructuralMatrix(skylnsym);
////    an.SetSolver(step);
////    
    an.Assemble();
    an.Rhs() ;
    
    TPZAutoPointer< TPZMatrix<REAL> > KGlobal;
    TPZFMatrix<STATE> FGlobal;
    KGlobal =   an.Solver().Matrix();
    FGlobal =   an.Rhs();
    
#ifdef PZDEBUG
//    #ifdef LOG4CXX
//            if(logger->isDebugEnabled())
//            {
//                std::stringstream sout;
//                KGlobal->Print("k = ", sout,EMathematicaInput);
//                FGlobal.Print("r = ", sout,EMathematicaInput);
//                LOGPZ_DEBUG(logger,sout.str())
//            }
//    #endif
#endif
    
////
    an.Solve();//assembla a matriz de rigidez (e o vetor de carga) global e inverte o sistema de equacoes
    
    TPZFMatrix<REAL> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
    solucao.Print("Sol",cout,EMathematicaInput);//imprime na formatacao do Mathematica

    
    
    // Post processing
    int ndiv = 2;
    TPZManVector<std::string> scalarnames(3), vecnames(1);
    scalarnames[0] = "SigmaX";
    scalarnames[1] = "SigmaY";
    scalarnames[2] = "SolidPressure";
    vecnames[0] = "Displacement";
    //vecnames[1] = "";
    an.DefineGraphMesh(2,scalarnames,vecnames,"ElasticitySolutions.vtk");
    
    an.PostProcess(ndiv);
//
//
//    
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
//
//    
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
//
//    
//
//    

}






// *********** Cria malha Geometrica Circular (360 graus) *************/

// Para gerar a malha eh necessario:
// rwb -> raio do poco
// re -> raio do contorno
// ncirc -> nro elem ao longo 1/4 do poco
// nrad -> nro elem da parede do poco ate contorno externo
// DrDcirc -> proporcao dos elementos da parede do poco

TPZGeoMesh *CircularGeoMesh (REAL rwb, REAL re, int ncirc, int nrad, REAL DrDcirc) {
    
    
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
    REAL firsttheta = (2*Pi) / (ncirc);
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
    
    
    //*********************** Malha Circunferencial (360 graus) *********************//
    
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
    
    // bc = -1 -> Normal Pressure condition na parede do poco
    for (int i = 0; i<ncirc; i++ ) {
        //        TPZGeoElBC gbc(gmesh->ElementVec()[i],7,-1);
        gmesh->ElementVec()[i]->CreateBCGeoEl(7, -1);
        //TPZGeoElBC gbc1(&elvec[i],7,-1);
    }
    
    
    // bc = -2 -> StressField condition circunferencia externa do farfield
    for (int i = 1; i<ncirc+1; i++ ) {
        gmesh->ElementVec()[(ncirc*nrad)-i]->CreateBCGeoEl(5, -2);
        //TPZGeoElBC gbc4(&elvec[(ncirc*nrad)-i],5,-4);
    }
    
    
   
    return gmesh;
    
}




//****************** Cria malha Computacional para malha 360 graus ***************************************//

TPZCompMesh *CircularCMesh(TPZGeoMesh *gmesh, int pOrder)
{
    int matId = 1;
    const int dim = 2; //dimensao do problema
    
    
    // Plane strain assumption
    // int linestrain = 1;
    
    //**************** Criando material  ********************************
    TPZMatElasticity2D *material = new TPZMatElasticity2D(matId);//criando material que implementa a formulacao fraca do problema modelo
    
    
    // Setting up paremeters
    REAL Eyoung = 1.0 , ni = 0.25, fbx = 0., fby = 0.;
    material->SetElasticity(Eyoung, ni, fbx, fby);
    
    
    /******* Calculating Inicial Stresses *******/
    // direction = direction/azimuth
    // inclination = wellbore inclination
    // problem assumption, inclined wellbore state = 1
    // Pwb = pressao da lama em MPa
    REAL Pi = 3.14159;
    REAL direction = 90., inclination = 30.; // graus
    REAL directionT   = direction*(Pi/180); // rad
    REAL inclinationT = inclination*(Pi/180); // rad
    int inclinedwellbore = 1;
    REAL Pwb = 38.71; // MPa
    
    // Tensoes in Situ, horizontais e vertical em MPa
    REAL SigmaVV = 51.94, Sigmahh = 34.79, SigmaHH = 42.14;
    // Seta os parametros do poco
    material->SetInclinedWellboreParameters(SigmaHH, Sigmahh, SigmaVV, directionT, inclinationT, inclinedwellbore);
    
    
    //Tensor de tensoes iniciais
    REAL SigmaX = 0., SigmaXY = 0., SigmaY = 0., SigmaZ = 0.;
    material->SetInclinedWellborePreStress(SigmaX, SigmaXY, SigmaY, SigmaZ);

#ifdef PZDEBUG
    #ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            
            std::stringstream out;
            out << " Stress rotation " << std::endl;
            out << " SigmaX     = " << SigmaX << std::endl;
            out << " SigmaXY    = " << SigmaXY << std::endl;
            out << " SigmaY     = " << SigmaY << std::endl;
            out << " SigmaZ     = " << SigmaZ << std::endl;
            LOGPZ_DEBUG(logger,out.str())
        }
    #endif
#endif
    
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmesh->SetDimModel(dim);//seta dimensao do modelo
    
    // Inserindo material na malha
    cmesh->InsertMaterialObject(material);
    
    
    // cond contorno
    const int bc0 = -1, bc1 = -2; // ids igual da malha geometrica
    //bc2 = -3, bc3 = -4, bc4 = -5, bc5 = -6;
    const int normalpressure = 6, stressfield = 4; // tipo de condicao de contorno
    //neumann = 1, mixed = 2;
    material->GetPreStress(SigmaX, SigmaXY, SigmaY, SigmaZ); // obtem tensoes iniciais
    
    ///Inserir condicao de contorno parede do poco
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    val1(0,0) = 1;
    val1(1,0) = 0;
    val1(0,1) = 0;
    val1(1,1) = 1;
    val2(0,0) = Pwb;
    val2(1,0) = Pwb;
    TPZMaterial * BCond0 = material->CreateBC(material, bc0, normalpressure, val1, val2);//cria material
    
    ///Inserir condicao de contorno circunferencia externa
    val1(0,0) = SigmaX;
    val1(1,0) = SigmaXY;
    val1(0,1) = SigmaXY;
    val1(1,1) = SigmaY;
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZMaterial * BCond1 = material->CreateBC(material, bc1, stressfield, val1, val2);//cria material
    
    
    cmesh->InsertMaterialObject(BCond0);//insere material na malha
    cmesh->InsertMaterialObject(BCond1);//insere material na malha
    
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



// *********** Cria malha Geometrica 1/4 do Poco *************/

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
//        TPZGeoElBC gbc(gmesh->ElementVec()[i],7,-1);
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


    // bc -5 -> Mixed
    TPZGeoElBC gbc1(gmesh->ElementVec()[0],0,-5);

    // bc -5 -> Mixed
    TPZGeoElBC gbc2(gmesh->ElementVec()[(ncirc*nrad)-1],1,-6);


    return gmesh;

}



// *********** Cria malha Computacional para 1/4 do Poco *************/


TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder)
{
    int matId = 1;
    const int dim = 2; //dimensao do problema


    // Plane strain assumption
   // int linestrain = 1;
    
    //**************** Criando material  ********************************
    TPZMatElasticity2D *material = new TPZMatElasticity2D(matId);//criando material que implementa a formulacao fraca do problema modelo
    
    
    // Setting up paremeters
    REAL Eyoung = 1.0 , ni = 0.25, fbx = 0., fby = 0.;
    material->SetElasticity(Eyoung, ni, fbx, fby);
    
    
    /******* Calculating Inicial Stresses *******/
    // direction = direction/azimuth
    // inclination = wellbore inclination
    // problem assumption, inclined wellbore state = 1
    REAL Pi = 3.14159;
    REAL direction = 90., inclination = 30.; //graus
    REAL directionT   = direction*(Pi/180); // rad
    REAL inclinationT = inclination*(Pi/180); // rad
    int inclinedwellbore = 1;
    
    // Tensoes in Situ, horizontais e vertical em MPa
    REAL SigmaVV = 51.94, Sigmahh = 34.79, SigmaHH = 42.14;
    // Seta os parametros do poco
    material->SetInclinedWellboreParameters(SigmaHH, Sigmahh, SigmaVV, directionT, inclinationT, inclinedwellbore);
    
    
    //Tensor de tensoes iniciais
    REAL SigmaX = 0., SigmaXY = 0., SigmaY = 0., SigmaZ = 0.;
    material->SetPreStress(SigmaX, SigmaXY, SigmaY, SigmaZ);
   
    
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmesh->SetDimModel(dim);//seta dimensao do modelo
    
    // Inserindo material na malha
    cmesh->InsertMaterialObject(material);
    
    
    // cond contorno
    const int bc0 = -1, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5, bc5 = -6; // ids igual da malha geometrica
    const int neumann = 1, normalpressure = 6, stressfield = 4, mixed = 2; // tipo de condicao de contorno
    material->GetPreStress(SigmaX, SigmaXY, SigmaY, SigmaZ); // obtem tensoes iniciais
    
    ///Inserir condicao de contorno parede do poco
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    val1(0,0) = SigmaX;
    val1(1,0) = SigmaXY;
    val1(0,1) = SigmaXY;
    val1(1,1) = SigmaY;
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZMaterial * BCond0 = material->CreateBC(material, bc0, stressfield, val1, val2);//cria material
    
    ///Inserir condicao de contorno externo bottom
    val1(0,0) = SigmaX;
    val1(1,0) = SigmaXY;
    val1(0,1) = SigmaXY;
    val1(1,1) = SigmaY;
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZMaterial * BCond1 = material->CreateBC(material, bc1, stressfield, val1, val2);//cria material
    
    ///Inserir condicao de contorno externo upper
    val1(0,0) = SigmaX;
    val1(1,0) = SigmaXY;
    val1(0,1) = SigmaXY;
    val1(1,1) = SigmaY;
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZMaterial * BCond2 = material->CreateBC(material, bc2, stressfield, val1, val2);//cria material
    
    ///Inserir condicao de contorno arco externo
    val1(0,0) = SigmaX;
    val1(1,0) = SigmaXY;
    val1(0,1) = SigmaXY;
    val1(1,1) = SigmaY;
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZMaterial * BCond3 = material->CreateBC(material, bc3, stressfield, val1, val2);//cria material
    
    ///Inserir condicao de contorno ponto parede do poco
    val1(0,0) = 1.;
    val1(1,0) = 0.;
    val1(0,1) = 0.;
    val1(1,1) = 1;
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZMaterial * BCond4 = material->CreateBC(material, bc4, mixed, val1, val2);//cria material
    
    ///Inserir condicao de contorno ponto arco externo do poco
    val1(0,0) = 0;
    val1(1,0) = 0;
    val1(0,1) = 0;
    val1(1,1) = 1;
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZMaterial * BCond5 = material->CreateBC(material, bc5, mixed, val1, val2);//cria material que implementa a condicao de contorno da parede do poco

    
    
//    //Condicao de contorno externo bottom
//    val2(0,0) = SigmaX;
//    val2(1,0) = SigmaY;
//    TPZMaterial * BCond1 = material->CreateBC(material, bc1, neumann, val1, val2);
//    
//    //Condicao de contorno externo bottom
//    val2(0,0) = SigmaX;
//    val2(1,0) = SigmaY;
//    TPZMaterial * BCond2 = material->CreateBC(material, bc2, neumann, val1, val2);
//    
//    //Condicao de contorno arco externo do farfield
//    val2(0,0) = SigmaX;
//    val2(1,0) = SigmaY;
//    TPZMaterial * BCond3 = material->CreateBC(material, bc3, stressfield, val1, val2);
//   
   
    cmesh->InsertMaterialObject(BCond0);//insere material na malha
    cmesh->InsertMaterialObject(BCond1);//insere material na malha
    cmesh->InsertMaterialObject(BCond2);//insere material na malha
    cmesh->InsertMaterialObject(BCond3);//insere material na malha
    cmesh->InsertMaterialObject(BCond4);//insere material na malha
    cmesh->InsertMaterialObject(BCond5);//insere material na malha
    
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




