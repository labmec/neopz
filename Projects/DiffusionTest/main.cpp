#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzartdiff.h"
#include "pzdiffmatrix.h"

void error(char * teste)
{

}

int main()
{

  const int dim = 3;
  const int nstate = dim+2;
  const int nphi = 6;

  // emulating state variables
  TPZVec<TPZDiffMatrix<REAL> > Ai, Tau;
  TPZVec<TPZVec<REAL> > TauDiv;
  TPZVec<TPZDiffMatrix<REAL> > TaudDiv;
  char ArtDiff[32]="LS";
  TPZFMatrix dsol(dim,nstate);
  TPZFMatrix dphi(dim,nphi);
  TPZVec<REAL> phi(nphi);
  TPZVec<REAL> sol(nstate);

  TPZVec<FADREAL> FADsol(nstate);
  TPZVec<FADREAL> FADdsol(nstate*dim);
  TPZVec<TPZVec<FADREAL> > FADTauDiv;
  // generating data

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

  cout << "\n\nFADsol\n" << FADsol;

  cout << "\n\nphi\n" << phi;

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

cout << "\n\nFADdsol\n" << FADdsol;

cout << "\n\ndphi\n" << dphi;

  TPZArtDiff::PrepareFastDiff(dim ,ArtDiff, sol, dsol, dphi, TauDiv, &TaudDiv);
  TPZArtDiff::PrepareFastDiff(dim ,ArtDiff, FADsol, FADdsol, FADTauDiv);


  cout << "\n\nFADTauDiv\n" << FADTauDiv;

  cout << "\n\nTauDiv\n" << TauDiv;

  cout << "\n\nTaudDiv\n" << TaudDiv;

  return 0;
}
