#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzartdiff.h"
#include "pzdiffmatrix.h"
#include "pzdiffmatrix.h"
#include "pzeulerconslaw.h"
#include "stdlib.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

int gDebug;
using namespace std;

void Flatten(TPZFMatrix<REAL> & coeff,
			 TPZFMatrix<REAL> &phi,
			 TPZFMatrix<REAL> &dphi,
			 TPZVec<REAL> &sol,
			 TPZFMatrix<REAL> & dsol);

void CheckConv(const double step,
			   TPZFMatrix<REAL> & coeff,
			   TPZFMatrix<REAL> &phi,
			   TPZFMatrix<REAL> &dphi,
			   TPZEulerConsLaw & MatTest);

#ifdef _AUTODIFF
void CheckJacobFlux(
					TPZEulerConsLaw & MatTest);
#endif

int main()
{
	const int dim = 3;
	const int nstate = dim+2;
	const int nphi = 6;
	
	gDebug = 0;
	
	TPZEulerConsLaw MatTest(0/*nummat*/,
							0.3/*timeStep*/,
							1.4 /*gama*/,
							dim,
							Bornhaus_AD);
	
	MatTest.SetTimeDiscr(Implicit_TD, Implicit_TD, Implicit_TD/*Implicit_TD, Implicit_TD*/);
	MatTest.SetContributionTime(Advanced_CT);
	MatTest.Print(cout);
	
    TPZMaterialData data;
	
	// emulating state variables
	data.dsol[0].Redim(dim,nstate);
	data.dphix.Redim(dim,nphi);
	data.phi.Redim(nphi,1);
	data.sol[0].Resize(nstate);
	
	TPZFMatrix<REAL> ef(nphi * nstate, 1, 0.);
	TPZFMatrix<REAL> ek(nphi * nstate, nphi * nstate, 0.);
	
	// generating data
	
	//phi
	data.phi(0)=2.3;
	data.phi(1)=3.5;
	data.phi(2)=5.7;
	data.phi(3)=7.11;
	data.phi(4)=11.13;
	data.phi(5)=13.17;
	
	//dphi
	int i;
	int j;
	for(i=0;i<dim;i++)
		for(j=0;j<nphi;j++)
			data.dphix(i,j)=45.8*i-3.2*j+2.; // any choice
	
	// individual solution coefficients
	TPZFMatrix<REAL> u(nstate * nphi,1, 0.);
	for(i = 0; i < nstate * nphi; i++)
	{
		u(i,0) = 10. * (double)(i+1) / (double)(i+2);
		if(i == nstate * nphi - 1) u(i,0) +=100.;// increasing the energy to avoid negative pressure
	}
	
	cout << "\nu" << u;
	
	//flattening solution
	Flatten(u, data.phi, data.dphix, data.sol[0], data.dsol[0]);
	
	cout << "\nsol" << data.sol[0];
	
	cout << "\ndsol" << data.dsol[0];
	
	data.jacinv.Redim(dim, dim);
	data.x.Resize(3);
	
	
	MatTest.Contribute(data,13., ek, ef);
	
	cout << "\nef\n" << ef;
	
	cout << "\nek\n" << ek;
	
#ifdef _AUTODIFF
	CheckJacobFlux(MatTest);
#endif
	
	cout.flush();
	
	CheckConv(.0001, u, data.phi, data.dphix, MatTest);
	
	return 0;
}

void Flatten(TPZFMatrix<REAL> & coeff,
			 TPZFMatrix<REAL> &phi,
			 TPZFMatrix<REAL> &dphi,
			 TPZVec<REAL> &sol,
			 TPZFMatrix<REAL> & dsol)
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
			   TPZFMatrix<REAL> & coeff,
			   TPZFMatrix<REAL> &phi,
			   TPZFMatrix<REAL> &dphi,
			   TPZEulerConsLaw & MatTest)
{
	
	cout << "\nCheckConv matrix\n";
	
	
	int nCoeff = coeff.Rows(), i, j;
	
	TPZFMatrix<REAL> deltaCoeff1(nCoeff, 1, 0.),
	deltaCoeff2(nCoeff, 1, 0.),
	updatedCoeff(nCoeff, 1, 0.);
	
	TPZFMatrix<REAL> Tangent(nCoeff, nCoeff, 0.),
	TrashTangent(nCoeff, nCoeff, 0.);
	TPZFMatrix<REAL> F0(nCoeff,1, 0.),
	F1(nCoeff,1, 0.),
	F2(nCoeff,1, 0.);
	
	int dim = dphi.Rows();
    TPZMaterialData data;
    data.phi = phi;
    data.dphix = dphi;
    data.jacinv.Resize(dim,dim);
	//   TPZFMatrix<REAL> jacinv(dim, dim, 7.);
	data.x.Resize(3);
	
	Flatten(coeff, data.phi, data.dphix, data.sol[0], data.dsol[0]);
	
	MatTest.Contribute(data, 13., Tangent, F0);
	
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
		Flatten(updatedCoeff, data.phi, data.dphix, data.sol[0], data.dsol[0]);
		MatTest.Contribute(data, 13., TrashTangent, F1);
		
		updatedCoeff = coeff + deltaCoeff2;
		Flatten(updatedCoeff, data.phi, data.dphix, data.sol[0], data.dsol[0]);
		MatTest.Contribute(data, 13., TrashTangent, F2);
		
		F1 -= F0;
		F2 -= F0;
		
		// matrix is already negative
		F1 += Tangent * deltaCoeff1;
		F2 += Tangent * deltaCoeff2;
		
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

#ifdef _AUTODIFF
// verifies if jacobian of Flux is right
void CheckJacobFlux(
					TPZEulerConsLaw & MatTest)
{
	int dim = 3;
	
	cout << "\nCheckJacobFlux on Flux Jacobian matrix\n";
	
	TPZVec<REAL> sol;
	TPZFMatrix<REAL> dsol;
	
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
	
	MatTest.JacobFlux(1.4, dim, U, Tangent);
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
#endif
