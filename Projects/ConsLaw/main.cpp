#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzartdiff.h"
#include "pzdiffmatrix.h"
#include "pzdiffmatrix.h"
#include "pzeulerconslaw.h"
#include "stdlib.h"
#include "pztempmat.h"


void error(char * teste)
{

}

void Flatten(TPZFMatrix & coeff,
	     TPZFMatrix &phi,
	     TPZFMatrix &dphi,
	     TPZVec<REAL> &sol,
	     TPZFMatrix & dsol);

void CheckConv(const double step,
	     TPZFMatrix & coeff,
	     TPZFMatrix &phi,
	     TPZFMatrix &dphi,
	     TPZEulerConsLaw2 & MatTest);

void CheckJacobFlux(
	     TPZEulerConsLaw2 & MatTest);

int main()
{
  TPZEulerConsLaw2 MatTest(0/*nummat*/,
		 0.3/*timeStep*/,
		 1.4 /*gama*/,
		 3 /*dim */,
		 LeastSquares_AD);

  MatTest.SetTimeDiscr(Implicit_TD, Implicit_TD, Implicit_TD/*Implicit_TD, Implicit_TD*/);
  MatTest.SetContributionTime(Advanced_CT);
  MatTest.Print(cout);

  const int dim = 3;
  const int nstate = dim+2;
  const int nphi = 6;

  // emulating state variables
  TPZFMatrix dsol(dim,nstate);
  TPZFMatrix dphi(dim,nphi);
  TPZFMatrix phi(nphi,1);
  TPZVec<REAL> sol(nstate);

  TPZFMatrix ef(nphi * nstate, 1);
  TPZFMatrix ek(nphi * nstate, nphi * nstate);

  TPZVec<FADREAL> FADsol(nstate);
  TPZVec<FADREAL> FADdsol(nstate*dim);

  // generating data

  //phi
  phi(0)=2.3;
  phi(1)=3.5;
  phi(2)=5.7;
  phi(3)=7.11;
  phi(4)=11.13;
  phi(5)=13.17;

  //dphi
  int i;
  int j;
  for(i=0;i<dim;i++)
     for(j=0;j<nphi;j++)
        dphi(i,j)=45.8*i-3.2*j+2.; // any choice

  // individual solution coefficients
  TPZFMatrix u(nstate * nphi,1);
  for(i = 0; i < nstate * nphi; i++)
  {
    u(i,0) = 10. * (double)(i+1) / (double)(i+2);
    if(i == 4) u(i,0) +=100.;// increasing the energy to avoid negative pressure
  }

   cout << "\nu" << u;

  //flattening solution
  Flatten(u, phi, dphi, sol, dsol);

  cout << "\nsol" << sol;

  cout << "\ndsol" << dsol;

  TPZFMatrix jacinv(dim, dim);
  TPZVec<REAL> x(3);

  MatTest.Contribute(x, jacinv, sol, dsol, 13, jacinv, phi, dphi, ek, ef);

cout << "\nef\n" << ef;

cout << "\nek\n" << ek;

  CheckJacobFlux(MatTest);

cout.flush();

  CheckConv(.01, u, phi, dphi, MatTest);

  return 0;
}

void Flatten(TPZFMatrix & coeff,
	     TPZFMatrix &phi,
	     TPZFMatrix &dphi,
	     TPZVec<REAL> &sol,
	     TPZFMatrix & dsol)
{
   int nshape = phi.Rows();
   int dim = dphi.Rows();
   int nstate = coeff.Rows()/nshape;

   sol.Resize(nstate);
   dsol.Redim(dim, nstate);

   int i_dim, i_shape, i_state;

   // flattening solution and dsol
   for(i_state = 0; i_state < nstate; i_state++)
   {
      sol[i_state] = 0;
      //for(i_dim = 0; i_dim < dim; i_dim++)
      //   dsol(i_dim, i_state)=0;
      for(i_shape = 0; i_shape < nshape; i_shape ++)
      {
         sol[i_state] += coeff(i_shape * nstate + i_state,0)
	            * phi(i_shape,0);
	 for(i_dim = 0; i_dim < dim; i_dim++)
	     dsol(i_dim, i_state) +=
	          dphi(i_dim, i_shape) *
		  coeff(i_shape * nstate + i_state,0);
      }
   }

}

void CheckConv(const double step,
	     TPZFMatrix & coeff,
	     TPZFMatrix &phi,
	     TPZFMatrix &dphi,
	     TPZEulerConsLaw2 & MatTest)
{

   cout << "\nCheckConv matrix\n";

   TPZVec<REAL> sol;
   TPZFMatrix dsol;

   int nCoeff = coeff.Rows(), i, j;

   TPZFMatrix deltaCoeff1(nCoeff),
              deltaCoeff2(nCoeff),
	      updatedCoeff(nCoeff);

   TPZFMatrix Tangent(nCoeff, nCoeff),
              TrashTangent(nCoeff, nCoeff);
   TPZFMatrix F0(nCoeff,1),
   	      F1(nCoeff,1),
   	      F2(nCoeff,1);

   int dim = dphi.Rows();
   TPZFMatrix jacinv(dim, dim);
   TPZVec<REAL> x(3);

   Flatten(coeff, phi, dphi, sol, dsol);

   MatTest.Contribute(x, jacinv, sol, dsol, 13, jacinv, phi, dphi, Tangent, F0);

/* To test the explicit last state contributions, remember to comment out
the contributions to T1 and T2;
   // testing explicit contributions
   MatTest.SetTimeDiscr(Explicit_TD, Explicit_TD, Explicit_TD);
   MatTest.SetContributionTime(Last_CT);
*/
   for(i = 0; i < nCoeff; i++)
   {
      deltaCoeff1.Zero();
      deltaCoeff2.Zero();

      REAL dC1, dC2;
      dC1 = (REAL)rand()/(REAL)RAND_MAX * step * coeff(i,0);
      dC2 = (REAL)rand()/(REAL)RAND_MAX * step * coeff(i,0);
      deltaCoeff1(i,0) = dC1;
      deltaCoeff2(i,0) = dC2;

      F1.Zero();
      F2.Zero();

      updatedCoeff = coeff + deltaCoeff1;
      Flatten(updatedCoeff, phi, dphi, sol, dsol);
      MatTest.Contribute(x, jacinv, sol, dsol, 13, jacinv, phi, dphi, TrashTangent, F1);

      updatedCoeff = coeff + deltaCoeff2;
      Flatten(updatedCoeff, phi, dphi, sol, dsol);
      MatTest.Contribute(x, jacinv, sol, dsol, 13, jacinv, phi, dphi, TrashTangent, F2);

      F1 -= F0;
      F2 -= F0;

      F1 -= Tangent * deltaCoeff1;
      F2 -= Tangent * deltaCoeff2;

      for(j = 0; j < nCoeff; j ++)
      {
         REAL diff1, diff2, lowestValue;
	 lowestValue = fabs(F0(j) * 1e-10);
	 diff1 = F1(j);
	 diff2 = F2(j);
	 if(fabs(diff1) < lowestValue && fabs(diff2) < lowestValue)
	 {
             cout << "\texact";
	 }else
	 {
	    double convRate = log10(diff1/diff2)/log10(dC1/dC2);
	    cout << "\t" << convRate;
	 }

      }
      cout << endl;

   }

}

// verifies if jacobian of Flux is right
void CheckJacobFlux(
	     TPZEulerConsLaw2 & MatTest)
{
   int dim = 3;

   cout << "\nCheckJacobFlux on Flux Jacobian matrix\n";

   TPZVec<REAL> sol;
   TPZFMatrix dsol;

   int nState = dim + 2, i, j, k;

   TPZVec< REAL > U(nState);
   TPZVec< FADREAL > FadU(nState);

   for(i = 0; i < nState; i++)
   {
      U[i] = 59./((double)i+5.5);
      if(i==4) U[i] +=20.;
      FadU[i].diff(i,nState);
      FadU[i].val() = U[i];
   }

   TPZVec< TPZDiffMatrix<REAL> > Tangent(dim),
                           FadTangent(dim),
			   Difference(dim);
   TPZVec< TPZVec< FADREAL > > FadFlux(3);

   MatTest.JacobFlux(U, Tangent);
   MatTest.Flux(FadU, FadFlux[0], FadFlux[1], FadFlux[2]);

cout << "aqui";

   for(k = 0; k < dim; k++)
   {
      FadTangent[k].Redim(nState, nState);
      Difference[k].Redim(nState, nState);
      for(i = 0; i < nState; i ++)
         for(j = 0; j < nState; j++)
	 {
	    FadTangent[k](i,j) = FadFlux[k][i].d(j);
	    Difference[k](i,j) =
	       FadTangent[k](i,j) -
	       Tangent[k](i,j);
	 }
   }

   cout << "\nTangent" << Tangent;

   cout << "\nFadTangent" << FadTangent;

   cout << "\nDifference" << Difference;
}

