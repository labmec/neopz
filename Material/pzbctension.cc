// $Id: pzbctension.cc,v 1.3 2003-11-04 00:26:33 phil Exp $

#include "pzbctension.h"
#include "pzadmchunk.h"
#include "pzintel.h"
#include "TPZMulticamadaOrtho.h"
//#include "TPZPlacaOrthotropic.h"

TPZBCTension::TPZBCTension(TPZMaterial *material,int id,int type,
			   TPZFMatrix &val1,TPZFMatrix &val2, TPZMultiCamadaOrthotropic *mult) :
  TPZBndCond(material,id,type,val1,val2) {
  fMx = 0.;
  fQxy = 0.;


  cout << "TPZBCTension::TPZBCTensio FALTA DEFINIR\n";  
}



void TPZBCTension::Tensor(TPZVec<REAL> &Point, int placa, int var, TPZManVector<REAL> &SolOut){

  if (var < 10 || var > 12)
    PZError << "TPZBCTension::Tensor variable error, var = " << var << endl;
  fMultCam->RPlacaOrtho()[placa]->InterpolateEl()->Solution(Point, var, SolOut);

}  

// void  TPZBCTension::FeedBCTension(){

  
// }

REAL TPZBCTension::MX(int placa,int face){

  TPZPlacaOrthotropic *cam =   fMultCam->CamadaN(placa);
  TPZInterpolatedElement *intel = cam->InterpolateEl();
  REAL f = fMultCam->FZ();
  REAL mx = intel->SolicitantEffort(f,cam->FH(),12,face);
  return mx;
}

REAL TPZBCTension::HIGHT(int placa){

  int i,j;
  REAL coordponto[3];
  cout << "digite a coordenada x  do ponto desejado" << " : " << endl;
  cin  >> coordponto [0];
  cout << "digite a coordenada y  do ponto desejado" << " : " << endl;
  cin  >> coordponto [1];

// void  TPZBCTension::FeedBCTension(){

  

  cout << "digite a coordenada z  do ponto desejado" << " : " << endl;
  cin >> coordponto [2];
  cout << endl;
  
  //  if () //módulo da altura da placa > módulo da altura do ponto


    }

REAL TPZBCTension::QXY(int placa, int face){


  TPZPlacaOrthotropic *cam =   fMultCam->CamadaN(placa);
  TPZInterpolatedElement *intel = cam->InterpolateEl();
  REAL f = fMultCam->FZ();
  REAL qxy = intel->SolicitantEffort(f,cam->FH(),12,face);
  return qxy;


}


REAL TPZBCTension::LXLYDIM(int placa){

}

REAL TPZBCTension::MY(int placa, int face){


  TPZPlacaOrthotropic *cam =   fMultCam->CamadaN(placa);
  TPZInterpolatedElement *intel = cam->InterpolateEl();
  REAL f = fMultCam->FZ();
  REAL my = intel->SolicitantEffort(f,cam->FH(),12,face);
  return my;
}


REAL TPZBCTension::MXY(int placa, int face){


  TPZPlacaOrthotropic *cam =   fMultCam->CamadaN(placa);
  TPZInterpolatedElement *intel = cam->InterpolateEl();
  REAL f = fMultCam->FZ();
  REAL mxy = intel->SolicitantEffort(f,cam->FH(),12,face);
  return mxy;
}



REAL TPZBCTension::NX(int placa, int face){


  TPZPlacaOrthotropic *cam =   fMultCam->CamadaN(placa);
  TPZInterpolatedElement *intel = cam->InterpolateEl();
  REAL f = fMultCam->FZ();
  REAL nx = intel->SolicitantEffort(f,cam->FH(),12,face);
  return nx;
}


REAL TPZBCTension::NY(int placa, int face){


  TPZPlacaOrthotropic *cam =   fMultCam->CamadaN(placa);
  TPZInterpolatedElement *intel = cam->InterpolateEl();
  REAL f = fMultCam->FZ();
  REAL ny = intel->SolicitantEffort(f,cam->FH(),12,face);
  return ny;
}

REAL TPZBCTension::NXY(int placa, int face){


  TPZPlacaOrthotropic *cam =   fMultCam->CamadaN(placa);
  TPZInterpolatedElement *intel = cam->InterpolateEl();
  REAL f = fMultCam->FZ();
  REAL nxy = intel->SolicitantEffort(f,cam->FH(),12,face);
  return nxy;
}


REAL TPZBCTension::QX(int placa, int face){


  TPZPlacaOrthotropic *cam =   fMultCam->CamadaN(placa);
  TPZInterpolatedElement *intel = cam->InterpolateEl();
  REAL f = fMultCam->FZ();
  REAL qx = intel->SolicitantEffort(f,cam->FH(),12,face);
  return qx;
}


REAL TPZBCTension::QY(int placa, int face){


  TPZPlacaOrthotropic *cam =   fMultCam->CamadaN(placa);
  TPZInterpolatedElement *intel = cam->InterpolateEl();
  REAL f = fMultCam->FZ();
  REAL qy = intel->SolicitantEffort(f,cam->FH(),12,face);
  return qy;
}
