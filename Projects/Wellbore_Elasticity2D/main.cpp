#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
#endif
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
    
    
    //int dim = 2;//dimensao do problema
    //REAL dom = 10.0; //comprimento do dominio unidimensional com inicio na origem zero
    //int nel = 2; //numero de elementos a serem utilizados
    //int pOrder = 1; //ordem polinomial de aproximacao
    //REAL elsize = dom/nel; //tamanho de cada elemento
    
    
    
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
        }
    }
    //Generate neighborhod information
    gmesh->BuildConnectivity();
    
    
    return gmesh;
    
}

    



// Cria malha COmputacional

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder)
{
    int matId = 1;
    const int dim = 2; //dimensao do problema


    // Plane strain assumption
    int linestrain = 1;
    
    //**************** Criando material  ********************************
    TPZMatElasticity2D *material = new TPZMatElasticity2D(matId);//criando material que implementa a formulacao fraca do problema modelo
    
    
    // Setting up paremeters
    REAL Eyoung = 1.0 , ni = 2.0, fbx = 0., fby = 0.;
    material->SetElasticity(Eyoung, ni, fbx, fby);
    
    
    /******* Calculating Inicial Stresses *******/
    // direction = direction/azimuth
    // inclination = wellbore inclination
    // inclined wellbore state = 1
    double direction = 0., inclination = 0.; //degrees
    int inclinedwellbore = 1;
    
    // Tensoes in Situ, horizontais e vertical
    REAL SigmaVV = 0., Sigmahh = 0., SigmaHH = 0.;
    
  
    // Seta as tensoes iniciais
    material->SetInclinedWellboreParameters(SigmaHH, Sigmahh, SigmaVV, direction, inclination, inclinedwellbore);
    
    ///criar malha computacional
    //TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    //cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    //cmesh->SetDimModel(dim);//seta dimensao do modelo
    
    
 //  return cmesh;
    
}




