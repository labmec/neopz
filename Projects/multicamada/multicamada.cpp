//oocalc &
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzplaca.h"
#include "pzintel.h"

#include "TPZGeoElement.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzelgc3d.h"
#include "pzbndcond.h"
#include "pztempmat.h"
#include "pzcompel.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzskylmat.h"
#include "pzstepsolver.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzstack.h"
#include "pzvec.h"
#include "pzsolve.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzelct2d.h"
#include "pzelcc3d.h"
#include "pzelgt3d.h"
#include "pzelct3d.h"
#include "pzelgpi3d.h"
#include "pzelcpi3d.h"
#include "pzelgpr3d.h"
#include "pzelcpr3d.h"
#include "pzelmat.h"
#include "pzelasmat.h"
#include "pzmattest.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzpoisson3d.h"
#include "pzmaterial.h"
#include "TPZConservationLaw.h"
#include "TPZConsLawTest.h"
#include "TPZEulerConsLaw.h"
#include "TPZDiffusionConsLaw.h"
#include "TPZCompElDisc.h"
#include "TPZShapeDisc.h"
#include "TPZInterfaceEl.h"
#include "pzreal.h"
#include "pzdxmesh.h"
#include "pzmatorthotropic.h"
#include "TPZMulticamadaOrtho.h"
#include "TPZPlacaOrthotropic.h"


#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <ostream>
#include <string.h>
#include <math.h>
using namespace std;


/**
 * -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <-
 */

int main(){

  //criando objetos das classes computacional e geométrica
  REAL zmin = 0.;
  REAL dx = 1.0;
  REAL dy = 1.0;
  int nelx = 1;
  int nely = 1;
  cout << "multicamada::main\n"
       << "dx " << dx <<  "\n"
       << "dy " << dy << "\n"
       << "nelx " << nelx << "\n"
       << "nely " << nely << "\n";
  //cin >> dx >> dy >> nelx >> nely;

  TPZMulticamadaOrthotropic *multcam = new TPZMulticamadaOrthotropic(zmin,dx,dy,nelx,nely);
  TPZMatOrthotropic *orto;

  //criando objeto para a classe material
  //REAL E = 200000.0, G = 173913.0, v = 0.15;
  //REAL E = 10000.0, G = 4347.83, v = 0.15;//teste
  //REAL E = 15.133e6, G = 2.59398e7, v = 0.471;//Pinus Pinaster ait.
  //REAL E = 69.0e6, G = 5.14378e6, v = 0.330;//Aluminio
  //REAL eppx=E, eppy=E, eppz=E, vxy=v, vyz=v,vzx=v, gxy=G, gyz=G, gzx=G;
  int numat = 1;
  TPZFMatrix naxes(3,3);
  naxes(0,0) =  1.0;    naxes(0,1) =  0.0;    naxes(0,2) =  0.0;
  naxes(1,0) =  0.0;    naxes(1,1) =  1.0;    naxes(1,2) =  0.0;
  naxes(2,0) =  0.0;    naxes(2,1) =  0.0;    naxes(2,2) =  1.0;
  REAL Elast = 100.0;
  REAL nu = 0.2;
  REAL eppx = Elast;
  REAL eppy = Elast;
  REAL eppz = Elast;
  REAL vxy = nu;
  REAL vyz = nu;
  REAL vzx = nu;
  REAL gxy = Elast/(2.*(1.+nu));
  REAL gyz = Elast/(2.*(1.+nu));
  REAL gzx = Elast/(2.*(1.+nu));

  orto = new TPZMatOrthotropic(numat, naxes, eppx, eppy, eppz, vxy, vyz, vzx, gxy, gyz, gzx);
  multcam->AddPlacaOrtho(orto,0.2);
  if(0){
    Elast = 200.;
    eppx = Elast;
    eppy = Elast;
    eppz = Elast;
    vxy = nu;
    vyz = nu;
    vzx = nu;
    gxy = Elast/(2.*(1.+nu));
    gyz = Elast/(2.*(1.+nu));
    gzx = Elast/(2.*(1.+nu));
    orto = new TPZMatOrthotropic(numat+1, naxes, eppx, eppy, eppz, vxy, vyz, vzx, gxy, gyz, gzx);
    multcam->AddPlacaOrtho(orto,0.1);
  }

  multcam->GenerateMesh();

  ofstream out("QX2.out");

  multcam->SetQX(1.);

  //    multcam->GeoMesh()->Print(out);
  //    multcam->CompMesh()->Print(out);
  out.flush();

  int niter = 20;
  int it;
  cout << "Numero de iterações" << endl;
  cin >> niter;

  TPZFMatrix tensin(5,9,0.),tensout(5,9);

  for(it=0; it<niter; it++) {

    multcam->ComputeSolution(out,0);
    multcam->ComputeCenterForces();
//     multcam->PrintCenterForces(out);
//     multcam->PrintTensors(out);
    multcam->PrintTensors(out,tensin,tensout);
    tensin = tensout;
  }

  out.close();

  return 0;
}

