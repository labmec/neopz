

#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "DarcyPTest.h"
#include "StokesTest.h"
#include "CoupledTest.h"

#include "TPZCouplingDSMaterial.h"
#include "TPZStokesMaterial.h"
#include "TPZDarcyPMaterial.h"
#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
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


const int SpaceHDiv = 1; //Velocidade em subespaço de H(div)
const int SpaceContinuous = 2; //Velocidade em subespaço de [H1]ˆ2
const int SpaceDiscontinuous = 3; //Velociadade em subespaço de H(Ph) - Ph: partição
const REAL visco=1., permeability=1., theta=-1.; //Coeficientes: viscosidade, permeabilidade, fator simetria
const REAL Pi=M_PI;

bool DarcyDomain = false, StokesDomain = true, CoupledDomain = false;

//Função principal do programa:

int main(int argc, char *argv[])
{
    
    TPZMaterial::gBigNumber = 1.e16;
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    //Dados do problema:
    
    
    int h_level = 8;
    
    double hx=1.,hy=1.; //Dimensões em x e y do domínio
    //double hx=Pi,hy=2.; //Dimensões em x e y do domínio (acoplamento)
    int nelx=h_level, nely=h_level; //Número de elementos em x e y
    int nx=nelx+1 ,ny=nely+1; //Número de nos em x  y
    int pOrder = 3; //Ordem polinomial de aproximação
    
  
    if (DarcyDomain) {
        DarcyPTest  * Test1 = new DarcyPTest();
        Test1->Run(SpaceHDiv, pOrder, nx, ny, hx, hy,visco,permeability,theta);
    }
    else if (StokesDomain)
    {
        StokesTest  * Test2 = new StokesTest();
        Test2->Run(SpaceHDiv, pOrder, nx, ny, hx, hy,visco,theta);
    }
    else  if(CoupledDomain)
    {
        CoupledTest  * Test3 = new CoupledTest();
        Test3->Run(SpaceHDiv, pOrder, nx, ny, hx, hy,visco,permeability,theta);
    }
    
    return 0;
}





