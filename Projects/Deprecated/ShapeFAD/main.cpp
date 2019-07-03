#include "pzreal.h"
#include "pzshapecube.h"
#include "pzshapequad.h"
#include <iostream>
#include "error.h"
#include "pzcompel.h"
#include <math.h>

void error(char * string)
{
  if(string)cerr << endl << string << endl;

}

int run(){

const int order_ = 20;
int numshape = (order_+1)*(order_+1)*(order_+1);

  TPZVec<REAL> point(3,0.);
  TPZVec<int> id(8);
  int i;

  for(i = 0; i< 8; i ++)
  {
  	id[i] = i;
  }

  TPZVec<int> order(19);
  for(i = 0; i< 19; i ++)
  {
        order[i] = order_;
  }

  TPZCompEl::SetgOrder(order_);

  TPZVec<FADREAL> phi(numshape);
  TPZFMatrix OldPhi(numshape,1), OldDPhi(3,numshape);
  TPZFMatrix DiffPhi(numshape,1), DiffDPhi(3,numshape);

  TPZShapeCube::ShapeCube(point, id, order, phi);
  TPZShapeCube::ShapeCube(point, id, order, OldPhi, OldDPhi);

  /*cout << "Calculated by Fad" << phi;
  cout << "Old derivative method (phi)\n" << OldPhi;
  cout << "Old derivative method (dPhi)\n" << OldDPhi;

  shapeFAD::ExplodeDerivatives(phi, DiffPhi, DiffDPhi);
  DiffPhi-=OldPhi;
  DiffDPhi-=OldDPhi;*/
  //cout << "FAD derivative method (phi)\n" << /*TPZFMatrix (OldPhi -*/ DiffPhi;
  //cout << "FAD derivative method (dPhi)\n" <</* TPZFMatrix (OldDPhi -*/ DiffDPhi;
  return 0;
}

int main()
{
 for (int i = 0; i < 1000 ; i++)run();
 return 0;
}
