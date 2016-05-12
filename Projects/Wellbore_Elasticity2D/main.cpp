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

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.Wellbore_Elasticity2D"));
#endif
using namespace std;


TPZGeoMesh *CreateGMesh(long nel, REAL elsize);
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
    
    
    /******** Configurar malha geometrica ***************/


    
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





// Cria malha Geometrica

TPZGeoMesh *CreateGMesh(long nel, REAL elsize)
{
    TPZGeoMesh * gmesh = new TPZGeoMesh;//Inicializa objeto da classe TPZGeoMesh
    
   
    
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




