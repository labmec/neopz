#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzartdiff.h"
#include "pzdiffmatrix.h"

int main()
{
#ifdef _AUTODIFF

#ifdef FAD
  const int dim = 3;
  const int nstate = dim+2;
  const int nphi = 6;

  // emulating state variables
  TPZVec<TPZDiffMatrix<REAL> > Ai, Tau;
  TPZVec<TPZVec<REAL> > TauDiv;
  TPZVec<TPZDiffMatrix<REAL> > TaudDiv;
  /*char ArtDiff[32]="LS";*/
  TPZArtDiffType artDiffType = LeastSquares_AD;
  TPZFMatrix dsol(dim,nstate);
  TPZFMatrix dphi(dim,nphi);
  TPZVec<REAL> phi(nphi);
  TPZVec<REAL> sol(nstate);


  TPZVec<FADREAL> FADsol(nstate);
  TPZVec<FADREAL> FADdsol(nstate*dim);
  TPZVec<TPZVec<FADREAL> > FADTauDiv;
  // generating data

  TPZArtDiff ArtDiff(artDiffType, 1.4);

  //solution
  sol[0]=2.;
  sol[1]=3.;
  sol[2]=5.;
  sol[3]=7.;
  sol[4]=11.;

  //phi
  phi[0]=2.3;
  phi[1]=3.5;
  phi[2]=5.7;
  phi[3]=7.11;
  phi[4]=11.13;
  phi[5]=13.17;

  //dphi
  int i;
  int j;
  for(i=0;i<dim;i++)
     for(j=0;j<nphi;j++)
        dphi(i,j)=45.8*i-3.2*j; // any choice

  //dsol
  for(i=0;i<dim;i++)
     for(j=0;j<nstate;j++)
        dsol(i,j)=49.8*i-3.1*j; // any choice

  int k, l;

  //FADsol
  for(i=0;i<nstate;i++)
     {
        FADsol[i]=sol[i];
        FADsol[i].diff(0,nphi*nstate);
	FADsol[i].fastAccessDx(0)=0.;
	for(j=0;j<nphi;j++)FADsol[i].fastAccessDx(i+j*nstate)=phi[j];
     }

/*  FADREAL teste;
  for(i=0;i<FADsol.NElements();i++)teste+=FADsol[i];
  cout << "\n\nFADsol\n" << teste;

  cout << "\n\nphi\n" << phi;
*/

  //FADdSol
  for(k=0;k<dim;k++)
     for(i=0;i<nstate;i++)
        {
            FADdsol[k+i*dim]=dsol(k,i);
	    FADdsol[k+i*dim].diff(0,nphi*nstate);
	    FADdsol[k+i*dim].fastAccessDx(0)=0.;
	    for(j=0;j<nphi;j++)
	    	    FADdsol[k+i*dim].fastAccessDx(i+j*nstate)=dphi(k,j);
	}
/*
	cout << "\n\nFADdsol\n";
teste=0;
for(i=0;i<FADdsol.NElements();i+=3)teste+=FADdsol[i];

cout << teste << endl;

teste=0;
for(i=1;i<FADdsol.NElements();i+=3)teste+=FADdsol[i];

cout << teste << endl;

teste=0;
for(i=2;i<FADdsol.NElements();i+=3)teste+=FADdsol[i];

cout << teste << endl;

//cout << "\n\nFADdsol\n" << FADdsol;

cout << "\n\ndphi\n" << dphi;
*/
  TPZFMatrix Jac(dim,dim,1.);
  ArtDiff.PrepareFastDiff(dim,Jac, sol, dsol, dphi, TauDiv, &TaudDiv);
  ArtDiff.PrepareFastDiff(dim, Jac, FADsol, FADdsol, FADTauDiv);


  cout << "\n\nFADTauDiv\n" << FADTauDiv;

  cout << "\n\nTauDiv\n" << TauDiv;

  cout << "\n\nTaudDiv\n" << TaudDiv;
#endif
	
#endif
  return 0;
}
