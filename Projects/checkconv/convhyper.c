
#include "convhyper.h"

TPZConvHyper::TPZConvHyper(int numnod, int num, REAL E, REAL mu, REAL nu, REAL lambda, REAL coef1, REAL coef2, REAL coef3) : TPZMatHyperElastic(num,E,mu,nu,lambda,coef1,coef2,coef3), fState(3,3,0.),
						       fPhi(numnod,1,0.),fDphi(3,numnod,0.),fAxes(3,3,0.),fSol(3,0.),fX(3,0.),fNumNod(numnod) {
  int i,n;
  for(i=0;i<3;i++) {
    fAxes(i,i) = 1.;
  }
  for(i=0; i<3; i++) {
    for(n=0; n<fNumNod; n++) {
      fDphi(i,n) = rand()%100;
    }
  }

}

int TPZConvHyper::NumCases() {
  return 1;
}

void TPZConvHyper::LoadSolution(TPZFMatrix &state) {
  int i,j,n;
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++) {
      fState(i,j) = 0.;
      for(n=0; n<fNumNod; n++) {
	fState(i,j) += fDphi(i,n)*state(j+n*3,0);
      }
    }
  }
  //fState.Print("state");
}

void TPZConvHyper::Residual(TPZFMatrix &residual, int icase) {
  TPZFMatrix ek(3*fNumNod,3*fNumNod);
  residual.Redim(3*fNumNod,1);
  int i;
  //for(i=0; i<3; i++) residual(i,0) = fState(0,i)*fState(0,i);
  Contribute(fX,fSol,fState,1.,fAxes,fPhi,fDphi,ek,residual);
  residual *= -1.;
  //residual.Print("residual");
}

void TPZConvHyper::ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &coefs, int icase) {
  TPZFMatrix ef(3*fNumNod,1);
  tangent.Redim(3*fNumNod,3*fNumNod);
  int i;
  //for(i=0; i<3; i++) tangent(i,i) = 2.*fState(0,i);
  Contribute(fX,fSol,fState,1.,fAxes,fPhi,fDphi,tangent,ef);
  //tangent.Print("tangent");
}
