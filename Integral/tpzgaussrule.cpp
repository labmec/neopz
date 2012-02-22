/**
 * @file
 * @brief Contains the implementation of the TPZGaussRule methods. 
 */

#include <cmath>

#include "tpzline.h"

#include "tpzgaussrule.h"
#include "tpzintrulelist.h"
#include "pzerror.h"

#include "pzvec.h"

using namespace std;

//***************************************
TPZGaussRule::TPZGaussRule(int order,int type,long double alpha,long double beta) {
	// Cleaning storage
	fLocation.Resize(0);
	fWeight.Resize(0);
	
	if(alpha <= -1.0L) alpha = 0.0L;  // it is required to right computation
	if(beta <= -1.0L) beta = 0.0L;    // it is required to right computation
	fAlpha = alpha;
	fBeta = beta;

	fType = type;
	
	if(order < 0) {
		PZError << "TPZGaussRule criation, order = " << order << " not allowable\n";
		fNumInt = 0;
		fLocation.Resize(0);
		fWeight.Resize(0);
		return;
	}
	switch (fType) {
		case 1:   // Lobatto
		{
			// Computing the points and weights for the symmetric Gaussian quadrature (using Legendre polynomials) 
			ComputingGaussLobattoQuadrature(order);
		}
			break;
		case 2:   // Jacobi
			if(IsZero(fAlpha) && IsZero(fBeta)) fAlpha = fBeta = 1.0L;   // Because fAlpha = fBeta = 0.0 is same as Gauss Legendre case
			// Computing the points and weights for the symmetric Gaussian quadrature (using Legendre polynomials) 
			ComputingGaussJacobiQuadrature(order,fAlpha,fBeta);
			break;
		case 3:  // Chebyshev
			fAlpha = fBeta = -0.5L;
			ComputingGaussChebyshevQuadrature(order);			
			break;
		default:     // Legendre (default)
			// Computing the points and weights for the symmetric Gaussian quadrature (using Legendre polynomials) 
			ComputingGaussLegendreQuadrature(order);
			break;
	}
	// Checks if the cubature rule is right
	if(!CheckCubatureRule())
		PZError << "TPZGaussRule had bad construction: order " << order << " type " << type << std::endl;

}

//***************************************
//***************************************
TPZGaussRule::~TPZGaussRule() {
	fNumInt = 0;
	fLocation.Resize(0);
	fWeight.Resize(0);	
}

//***************************************
//***************************************
long double TPZGaussRule::Loc(int i) const {
	if (i>=0 && i<fNumInt)
		return fLocation[i];
	else {
		PZError << "ERROR(TPZGaussRule::loc) Out of bounds!!\n";
		return 0.0L;
	}
}

//***************************************
//***************************************
long double TPZGaussRule::W(int i) const {
	if (i>=0 && i<fNumInt)
		return fWeight[i];
	else {
		PZError << "ERROR(TPZGaussRule::w) Out of bounds!!\n";
		return 0.0L;
	}
}

void TPZGaussRule::Print(std::ostream &out) {
	int i;
	out << "\nGaussian Quadrature ";
	if(!fType) out << "(Legendre) : " << std::endl;
	else if(fType == 1) out << "(Lobatto) : " << std::endl;
	else if(fType == 2) out << "(Jacobi) : " << std::endl;
	else if(fType == 3) out << "(Chebyshev) : " << std::endl;
	else out << "(Unkown) : " << std::endl;
	out << "\nNumber of points : " << fNumInt  << std::endl;
	for(i=0;i<fNumInt;i++) {
		out << i << " :\t" << fLocation[i] << " \t" << fWeight[i] << std::endl;
	}
	out << std::endl;
}

void TPZGaussRule::SetType(int &type,int order) {
	if(order < 0) order = 1;
	if(type < 0 || type > 3) type = 0;
	if(fType != type) {
		fLocation.Resize(0);
		fWeight.Resize(0);
		fType = type;
		switch(type) {
			case 1:   // Lobatto
			{
				// Computing the points and weights for the symmetric Gaussian quadrature (using Legendre polynomials) 
				ComputingGaussLobattoQuadrature(order);
			}
				break;
			case 2:   // Jacobi
			{
				// Computing the points and weights for the symmetric Gaussian quadrature (using Legendre polynomials) 
				ComputingGaussJacobiQuadrature(order,fAlpha,fBeta);
			}
				break;
			case 3:  // Chebyshev
			{
				ComputingGaussChebyshevQuadrature(order);			
			}
				break;
			default:     // Legendre (default)
			{
				// Computing the points and weights for the symmetric Gaussian quadrature (using Legendre polynomials) 
				ComputingGaussLegendreQuadrature(order);
			}
				break;
		}
	}
}

/** GAUSS LEGENDRE QUADRATURE */

// Compute the points and weights for Gauss Legendre Quadrature over the parametric 1D element [-1.0, 1.0] - This quadrature is symmetric
void TPZGaussRule::ComputingGaussLegendreQuadrature(int order) {
	// Computing number of points appropriated to the wished order = 2*npoints - 1. Note: If odd we need increment one point more
	fNumInt = (int)(0.51*(order+2));
	if(fNumInt < 1)
		fNumInt = 1;
	if(fNumInt > TPZGaussRule::NRULESLEGENDRE_ORDER) 
		fNumInt = TPZGaussRule::NRULESLEGENDRE_ORDER-1;

	
	long double tol = 1.e-19L;
	long double z1, z, pp, p3, p2, p1, dif, den;
	int i, j;
	long iteration;
	
	// Cleaning vector to storage
	fLocation.Resize(0);
	fWeight.Resize(0);
	fLocation.Resize(fNumInt,0.0L);
	fWeight.Resize(fNumInt,0.0L);
	
	int m = (fNumInt+1)*0.5;

	for(i=0;i<m;i++) {
		p1 = (long double)(i + 0.75L);
		p2 = (long double)(fNumInt+0.5L);
		z = cosl((M_PI*p1)/p2);   // p2 never is zero
		iteration = 0L;
		do {
			iteration++;
			p1 = 1.0L;
			p2 = 0.0L;
			for(j=0;j<fNumInt;j++) {
				p3 = p2;
				p2 = p1;
				p1 = (((2.0L)*((long double)j)+(1.0L))*z*p2-(((long double)(j))*p3))/(((long double)(j))+(1.0L));   // denominator never will be zero
			}
			den = z*z - (1.0L);
			if(fabs(den)<1.e-16L)
				z = 0.5L;
			pp = ((long double)fNumInt)*(z*p1-p2)/den;
			z1 = z;
			if(fabs(pp)<1.e-16L)
				z = 0.5L;
			else z = z1-p1/pp;
			dif = fabs(z - z1);
		}while(fabs(dif) > tol && iteration < PZINTEGRAL_MAXITERATIONS_ALLOWED);

		// If the maxime number of iterations is reached then the computing is stopped
		if(iteration == PZINTEGRAL_MAXITERATIONS_ALLOWED) {
			std::cout << "Reached maxime number of iterations for NPOINTS = " << fNumInt << ".\n";
			fNumInt = 0;
			fLocation.Resize(0);
			fWeight.Resize(0);
			break;
		}
		
		fLocation[i] = -z;
		fLocation[fNumInt-1-i] = z;
		fWeight[i] = 2.0L/((1.0L-z*z)*pp*pp);
		fWeight[fNumInt-1-i] = fWeight[i];
	}
}
// Compute the points and weights for Gauss Legendre Quadrature over the parametric 1D element [-1.0, 1.0] - This quadrature is symmetric
void TPZGaussRule::ComputingGaussLegendreQuadrature(int npoints,TPZVec<long double> &Location,TPZVec<long double> &Weight) {
	long double tol = 1.e-19L;
	long double z1, z, pp, p3, p2, p1, dif, den;
	int i, j;
	long iteration;
	
	// Cleaning vector to storage
	Location.Resize(0);
	Weight.Resize(0);
	Location.Resize(npoints,0.0L);
	Weight.Resize(npoints,0.0L);
	
	int m = (npoints+1)*0.5;
	
	for(i=0;i<m;i++) {
		p1 = ((long double)(i))+(0.75L);
		p2 = ((long double)(npoints))+0.5L;
		z = cosl((M_PI*p1)/p2);   // p2 never is zero
		iteration = 0L;
		do {
			iteration++;
			p1 = 1.0L;
			p2 = 0.0L;
			for(j=0;j<npoints;j++) {
				p3 = p2;
				p2 = p1;
				p1 = (((2.0L)*((long double)(j))+(1.0L))*z*p2-(((long double)(j))*p3))/(((long double)(j))+(1.0L));   // denominator never will be zero
			}
			den = (z*z)-(1.0L);
			if(fabs(den)<1.e-16L)
				z = 0.5L;
			pp = ((long double)npoints)*(z*p1-p2)/den;
			z1 = z;
			if(fabs(pp)<1.e-16L)
				z = 0.5L;
			else z = z1-p1/pp;
			dif = fabs(z - z1);
		}while(fabs(dif) > tol && iteration < PZINTEGRAL_MAXITERATIONS_ALLOWED);
		
		// If the maxime number of iterations is reached then the computing is stopped
		if(iteration == PZINTEGRAL_MAXITERATIONS_ALLOWED) {
			std::cout << "Reached maxime number of iterations for NPOINTS = " << npoints << ".\n";
			npoints = 0;
			Location.Resize(0);
			Weight.Resize(0);
			break;
		}
		
		Location[i] = -z;
		Location[npoints-1-i] = z;
		Weight[i] = 2.0L/((1.0L-z*z)*pp*pp);
		Weight[npoints-1-i] = Weight[i];
	}
}

void TPZGaussRule::ComputingGaussChebyshevQuadrature(int order) {
	// Computing number of points appropriated to the wished order = 2*npoints - 1. Note: If odd we need increment one point more
	fNumInt = (int)(0.51*(order+2));
	if(fNumInt < 0)
		fNumInt = 1;

	// Cleaning vector to storage
	fLocation.Resize(0);
	fWeight.Resize(0);
	fLocation.Resize(fNumInt, 0.0L);
	fWeight.Resize(fNumInt, 0.0L);
	
	for(int i=0;i<fNumInt;i++) {
		fLocation[i] = cosl(M_PI*(long double)(i+0.5L)/((long double)(fNumInt)));
		fWeight[i] = M_PI/((long double)(fNumInt));
	}
}

// To compute the polinomial Jacobi on x point, n is the order, alpha and beta are the exponents as knowed
long double TPZGaussRule::JacobiPolinomial(long double x,int alpha,int beta,unsigned int n) {
	// the Jacobi polynomial is evaluated
	// using a recursion formula.
	TPZVec<long double> p(n+1);
	int v, a1, a2, a3, a4;
	
	// initial values P_0(x), P_1(x):
	p[0] = 1.0L;
	if(n==0) return p[0];
	p[1] = (((long double)(alpha+beta+2.0L))*x + ((long double)(alpha-beta)))/2.0L;
	if(n==1) return p[1];
	
	for(unsigned int i=1; i<=(n-1); ++i) {
		v  = 2*i + alpha + beta;
		a1 = 2*(i+1)*(i + alpha + beta + 1)*v;
		a2 = (v + 1)*(alpha*alpha - beta*beta);
		a3 = v*(v + 1)*(v + 2);
		a4 = 2*(i+alpha)*(i+beta)*(v + 2);
		
		p[i+1] = static_cast<long double>( (((long double)(a2)) + ((long double)(a3))*x)*p[i] - ((long double)(a4))*p[i-1])/((long double)(a1));
    } // for
	return p[n];
}

/** Computes the points and weights for Gauss Lobatto quadrature */
void TPZGaussRule::ComputingGaussLobattoQuadrature(int order) {
	// Computing number of points appropriated to the wished order = 2*npoints - 1. Note: If odd we need increment one point more
	fNumInt = (int)(0.51*(order+2));
	fNumInt++;
	
	if(fNumInt < 2)
		fNumInt = 2;
	if(fNumInt > NRULESLOBATTO_ORDER) 
		fNumInt = NRULESLOBATTO_ORDER-1;
	
	// Cleaning vector to storage
	fLocation.Resize(0);
	fWeight.Resize(0);
	fLocation.Resize(fNumInt, 0.0L);
	fWeight.Resize(fNumInt, 0.0L);
	
	// Process to compute inner integration points at [-1.,1.]
	const unsigned int m = fNumInt - 2; // no. of inner points
	
	unsigned int counter = 0;
	// compute quadrature points with a Newton algorithm.
	
	// Set tolerance.
	const long double eps = 1.09e-19L;
    long double deps = 2.23e-16L;
	
	// check whether long double is more accurate than double, and
	// set tolerances accordingly
	const long double tol = ((1.L + eps) != 1.L ? std::max(deps/(100.0L), eps*5.0L) : deps*5.0L);
	
	// it take the zeros of the Chebyshev polynomial (alpha=beta=-0.5) as initial values:
	unsigned int i, k;
	for(i=0; i<m; ++i)
		fLocation[i] = - cosl(((2.0L*((long double)(i))+1.0L)/(2.0L*((long double)(m))))*M_PI);
	
	long double r, s, J_x, f, delta;
	
	for(k=0; k<m; ++k) {
		r = fLocation[k];
		if(k>0) r = (r + fLocation[k-1])/2.0L;
		do {
			counter++;
			s = 0.0L;
			for (i=0; i<k; ++i)
				s += 1.0L/(r - fLocation[i]);
			
			J_x   =  0.5L*((long double)(m + 3.0L))*JacobiPolinomial(r,2,2, m-1);
			f     = JacobiPolinomial(r,1,1, m);
			delta = f/(f*s- J_x);
			r += delta;
		}
		while (fabs(delta) >= tol && counter < PZINTEGRAL_MAXITERATIONS_ALLOWED);
		
		// If maxime iterations is reached, clean all points
		if(counter == PZINTEGRAL_MAXITERATIONS_ALLOWED) {
			fNumInt = 0;
			fWeight.Resize(0);
			fLocation.Resize(0);
			cout << "Maxime number of iterations allowed in quadrature GaussLobatto. \n" ;
			return;
		}
		
		fLocation[k] = r;
    }
	
	int ii;
	for(ii=m-1;ii>-1;ii--)
		fLocation[ii+1] = fLocation[ii];
	
	fLocation[0] = -1.0L;
	fLocation[fNumInt-1] = 1.0L;
	
	// Process to compute the weights to corresponding integration points
	s = 0.L;
	i = (unsigned int) fNumInt;
	f = gamma(i);
	r = gamma(i+1);
	const long double factor = 2.L*f*f / (((long double)(fNumInt-1))*f*r);
	
	fWeight.Resize(fNumInt,0.0L);
	
	for(ii=0; ii < fNumInt; ++ii)
    {
		s = JacobiPolinomial(fLocation[ii],0,0, fNumInt-1);
		fWeight[ii] = factor/(s*s);
    }
}

/** Gauss Jacobi quadrature */
void TPZGaussRule::ComputingGaussJacobiQuadrature(int order,long double alpha, long double beta) {
	// Computing number of points appropriated to the wished order = 2*npoints - 1. Note: If odd we need increment one point more
	fNumInt = (int)(0.51*(order+2));
	if(fNumInt < 0)
		fNumInt = 1;
	long double an;
	TPZVec<long double> b;
	long double bn;
	TPZVec<long double> c;
	long double cc, delta, dp2;
	int i;
	long double p1, prod, r1, r2, r3, temp, x0;
	
	fLocation.Resize(fNumInt,0.0L);
	fWeight.Resize(fNumInt,0.0L);
	
	b.Resize(fNumInt,0.0L);
	c.Resize(fNumInt,0.0L);

	//  This method permit only alpha > -1.0 and beta > -1.0
	if(alpha <= -1.0L) {
		std::cerr << "\nJACOBI_COMPUTE - Fatal error!\n  -1.0 < ALPHA is required.\n";
		std::exit ( 1 );
	}
	if(beta <= -1.0L) {
		std::cerr << "\nJACOBI_COMPUTE - Fatal error!\n  -1.0 < BETA is required.\n";
		std::exit ( 1 );
	}

	//  Set the recursion coefficients.
	for(i=1;i<=fNumInt;i++) {
		if(alpha + beta == 0.0L || beta - alpha == 0.0L) {
			b[i-1] = 0.0L;
		}
		else {
			b[i-1] = ( alpha + beta ) * ( beta - alpha ) / 
			( ( alpha + beta + ( 2.0L * i ) ) 
             * ( alpha + beta + ( 2.0L * i - 2.0L ) ) );
		}
		if(i==1) {
			c[i-1] = 0.0L;
		}
		else {
			c[i-1] = 4.0L * ( i - 1.0L ) * ( alpha + ( i - 1.0L ) ) * ( beta + ( i - 1.0L ) ) 
            * ( alpha + beta + ( i - 1.0L ) ) / 
            ( ( alpha + beta + ( 2.0L * i - 1.0L ) ) * (std::pow(alpha + beta + (2.0L*i - 2.0L),2.0L)) 
			 * ( alpha + beta + ( 2.0L * i - 3.0L ) ) );
		}
	}
	
	delta = (gamma(alpha + 1.0L) * gamma(beta + 1.0L)) / gamma(alpha + beta + 2.0L);
	
	prod = 1.0L;
	for(i=2;i <= fNumInt;i++)
		prod = prod * c[i-1];
	cc = delta * std::pow( 2.0L, alpha + beta + 1.0L ) * prod;
	
	for(i=1;i<=fNumInt;i++) {
		if(i==1) {
			an = alpha / (long double)(fNumInt);
			bn = beta / (long double)(fNumInt);
			
			r1 = ( 1.0L + alpha ) * ( 2.78L / ( 4.0L + (long double)(fNumInt*fNumInt)) + 0.768L * an / (long double)(fNumInt));
			r2 = 1.0L + 1.48L * an + 0.96L * bn + 0.452L * an * an + 0.83L * an * bn;
			
			x0 = ( r2 - r1 ) / r2;
		}
		else if(i==2) {
			r1 = ( 4.1L + alpha ) / ( ( 1.0L + alpha ) * ( 1.0L + 0.156L * alpha ) );
			
			r2 = 1.0L + 0.06L * ((long double)(fNumInt) - 8.0L) * (1.0L + 0.12L * alpha) / (long double)(fNumInt);
			r3 = 1.0L + 0.012L * beta * (1.0L + 0.25L * fabsl(alpha)) / (long double)(fNumInt);

			x0 = x0 - r1*r2*r3*(1.0L - x0);
		}
		else if(i==3) {
			r1 = ( 1.67L + 0.28L * alpha ) / ( 1.0L + 0.37L * alpha );
			r2 = 1.0L + 0.22L * ( (long double)(fNumInt) - 8.0L) / (long double)(fNumInt);
			r3 = 1.0L + 8.0L * beta / ( ( 6.28L + beta ) * (long double)(fNumInt*fNumInt));
			
			x0 = x0 - r1 * r2 * r3 * ( fLocation[0] - x0 );
		}
		else if(i<fNumInt-1) {
			x0 = 3.0L * fLocation[i-2] - 3.0L * fLocation[i-3] + fLocation[i-4];
		}
		else if(i==fNumInt-1) {
			r1 = ( 1.0L + 0.235L * beta ) / (0.766L + 0.119L * beta);
			r2 = 1.0L / ( 1.0L + 0.639L * ((long double)(fNumInt) - 4.0L) 
						 / ( 1.0L + 0.71L * ((long double)(fNumInt) - 4.0L)));
			
			r3 = 1.0L / ( 1.0L + 20.0L * alpha / ( ( 7.5L + alpha ) * 
												  (long double)(fNumInt*fNumInt)));
			
			x0 = x0 + r1 * r2 * r3 * ( x0 - fLocation[i-3] );
		}
		else if(i==fNumInt) {
			r1 = ( 1.0L + 0.37L * beta ) / ( 1.67L + 0.28L * beta );
			
			r2 = 1.0L / ( 1.0L + 0.22L * ( (long double)(fNumInt) - 8.0L ) 
			 / (long double)(fNumInt));
			
			r3 = 1.0L / ( 1.0L + 8.0L * alpha / 
						 ( ( 6.28L + alpha ) * (long double)(fNumInt*fNumInt)));
			
			x0 = x0 + r1 * r2 * r3 * ( x0 - fLocation[i-3] );
		}
		// Improving the root of the Jacobi polinomial
		{
			long double d;
			long double precision;
			int itera, itera_max = 10;
			
			precision = machinePrecision();
			
			for(itera=1;itera<=itera_max;itera++) {
				d = JacobiPolinomial(x0,fNumInt,alpha,beta,b,c,&dp2,&p1);
				d /= dp2;
				x0 = x0 - d;
				if(fabsl(d) <= precision*(fabsl(x0) + 1.0L)) {
					break;
				}
			}
		}
		fLocation[i-1] = x0;
		fWeight[i-1] = cc/(dp2*p1);
	}
	// inverting the position of the points
	for(i=1;i<=fNumInt/2;i++) {
		temp = fLocation[i-1];
		fLocation[i-1] = fLocation[fNumInt-i];
		fLocation[fNumInt-i] = temp;
	}
	// inverting the position of the weights
	for(i=1;i<=fNumInt/2;i++) {
		temp = fWeight[i-1];
		fWeight[i-1] = fWeight[fNumInt-i];
		fWeight[fNumInt-i] = temp;
	}
}
/** Gauss Jacobi quadrature */
void TPZGaussRule::ComputingGaussJacobiQuadrature(int npoints,long double alpha, long double beta,TPZVec<long double> &Location,TPZVec<long double> &Weight) {
	// Computing number of points appropriated to the wished order = 2*npoints - 1. Note: If odd we need increment one point more
	long double an, bn;
	TPZVec<long double> b;
	TPZVec<long double> c;
	long double cc, delta, dp2;
	int i;
	long double p1, prod, r1, r2, r3, temp, x0;
	
	Location.Resize(npoints,0.0L);
	Weight.Resize(npoints,0.0L);
	
	b.Resize(npoints,0.0L);
	c.Resize(npoints,0.0L);
	
	//  This method permit only alpha > -1.0 and beta > -1.0
	if(alpha <= -1.0L) {
		std::cerr << "\nJACOBI_COMPUTE - Fatal error!\n  -1.0 < ALPHA is required.\n";
		std::exit ( 1 );
	}
	if(beta <= -1.0L) {
		std::cerr << "\nJACOBI_COMPUTE - Fatal error!\n  -1.0 < BETA is required.\n";
		std::exit ( 1 );
	}
	
	//  Set the recursion coefficients.
	for(i=1;i<=npoints;i++) {
		if(alpha + beta == 0.0L || beta - alpha == 0.0L) {
			b[i-1] = 0.0L;
		}
		else {
			b[i-1] = ( alpha + beta ) * ( beta - alpha ) / 
			( ( alpha + beta + ( 2.0L * ((long double)(i)) ) ) 
             * ( alpha + beta + ( 2.0L * ((long double)(i)) - 2.0L ) ) );
		}
		if(i==1) {
			c[i-1] = 0.0L;
		}
		else {
			c[i-1] = 4.0L * ( ((long double)(i)) - 1.0L ) * ( alpha + ( ((long double)(i)) - 1.0L ) ) * ( beta + ( ((long double)(i)) - 1.0L ) ) 
            * ( alpha + beta + ( ((long double)(i)) - 1.0L ) ) / 
            ( ( alpha + beta + ( 2.0L * ((long double)(i)) - 1.0L ) ) * (std::pow(alpha + beta + (2.0L*((long double)(i)) - 2.0L),2.0L)) 
			 * ( alpha + beta + ( 2.0L * ((long double)(i)) - 3.0L ) ) );
		}
	}
	
	delta = (gamma(alpha + 1.0L) * gamma(beta + 1.0L)) / gamma(alpha + beta + 2.0L);
	
	prod = 1.0L;
	for(i=2;i <= npoints;i++)
		prod = prod * c[i-1];
	cc = delta * std::pow( 2.0L, alpha + beta + 1.0L ) * prod;
	
	for(i=1;i<=npoints;i++) {
		if(i==1) {
			an = alpha / (long double)(npoints);
			bn = beta / (long double)(npoints);
			
			r1 = ( 1.0L + alpha ) * ( 2.78L / ( 4.0L + (long double)(npoints*npoints)) + 0.768L * an / (long double)(npoints));
			r2 = 1.0L + 1.48L * an + 0.96L * bn + 0.452L * an * an + 0.83L * an * bn;
			
			x0 = ( r2 - r1 ) / r2;
		}
		else if(i==2) {
			r1 = ( 4.1L + alpha ) / ( ( 1.0L + alpha ) * ( 1.0L + 0.156L * alpha ) );
			
			r2 = 1.0L + 0.06L * ((long double)(npoints) - 8.0L) * (1.0L + 0.12L * alpha) / (long double)(npoints);
			r3 = 1.0L + 0.012L * beta * (1.0L + 0.25L * fabsl(alpha)) / (long double)(npoints);
			
			x0 = x0 - r1*r2*r3*(1.0L - x0);
		}
		else if(i==3) {
			r1 = ( 1.67L + 0.28L * alpha ) / ( 1.0L + 0.37L * alpha );
			r2 = 1.0L + 0.22L * ( (long double)(npoints) - 8.0L) / (long double)(npoints);
			r3 = 1.0L + 8.0L * beta / ( ( 6.28L + beta ) * (long double)(npoints*npoints));
			
			x0 = x0 - r1 * r2 * r3 * ( Location[0] - x0 );
		}
		else if(i<npoints-1) {
			x0 = 3.0L * Location[i-2] - 3.0L * Location[i-3] + Location[i-4];
		}
		else if(i==npoints-1) {
			r1 = ( 1.0L + 0.235L * beta ) / (0.766L + 0.119L * beta);
			r2 = 1.0L / ( 1.0L + 0.639L * ((long double)(npoints) - 4.0L) 
						 / ( 1.0L + 0.71L * ((long double)(npoints) - 4.0L)));
			
			r3 = 1.0L / ( 1.0L + 20.0L * alpha / ( ( 7.5L + alpha ) * 
												  (long double)(npoints*npoints)));
			
			x0 = x0 + r1 * r2 * r3 * ( x0 - Location[i-3] );
		}
		else if(i==npoints) {
			r1 = ( 1.0L + 0.37L * beta ) / ( 1.67L + 0.28L * beta );
			
			r2 = 1.0L / ( 1.0L + 0.22L * ( (long double)(npoints) - 8.0L ) 
						 / (long double)(npoints));
			
			r3 = 1.0L / ( 1.0L + 8.0L * alpha / 
						 ( ( 6.28L + alpha ) * (long double)(npoints*npoints)));
			
			x0 = x0 + r1 * r2 * r3 * ( x0 - Location[i-3] );
		}
		
		// Improving the root of the Jacobi polinomial
		long double d;
		long double precision;
		int itera, itera_max = 10;
			
		precision = machinePrecision();
		
		for(itera=1;itera<=itera_max;itera++) {
			d = JacobiPolinomial(x0,npoints,alpha,beta,b,c,&dp2,&p1);
			d /= dp2;
			x0 = x0 - d;
			if(fabsl(d) <= precision*(fabsl(x0) + 1.0L)) {
				break;
			}
		}
		Location[i-1] = x0;
		Weight[i-1] = cc/(dp2*p1);
	}
	// inverting the position of the points
	for(i=1;i<=npoints/2;i++) {
		temp = Location[i-1];
		Location[i-1] = Location[npoints-i];
		Location[npoints-i] = temp;
	}
	// inverting the position of the weights
	for(i=1;i<=npoints/2;i++) {
		temp = Weight[i-1];
		Weight[i-1] = Weight[npoints-i];
		Weight[npoints-i] = temp;
	}
}

/** Evaluate the Jacobi polinomial with alpha and beta real numbers */
long double TPZGaussRule::JacobiPolinomial(long double x, int order,long double alpha, long double beta, 
							TPZVec<long double> &b, TPZVec<long double> &c,long double *dp2, long double *p1 )
{
	long double p2, p0, dp0, dp1;
	int i;
	
	*p1 = 1.0L;
	dp1 = 0.0L;	
	p2 = x + (alpha - beta)/(alpha + beta + 2.0L);
	*dp2 = 1.0L;
	
	for(i=2;i<=order;i++) {
		p0 = *p1;
		dp0 = dp1;
		
		*p1 = p2;
		dp1 = *dp2;
		
		p2 = (x - b[i-1])*(*p1) - c[i-1]*p0;
		*dp2 = (x - b[i-1])*dp1 + (*p1) - c[i-1]*dp0;
	}
	return p2;
}

/** Returns the machine precision of the computer's arithmetic as long double number */
long double machinePrecision() {
	long double precision = 1.0L;
	int iterations = 0;
	
	while ( 1.0L < (long double) (1.0L + precision)  && iterations < PZINTEGRAL_MAXITERATIONS_ALLOWED) {
		precision = precision / 2.0L;
		if(iterations == PZINTEGRAL_MAXITERATIONS_ALLOWED) {
			cout << "\nMaxime iterations reached - epsilon() \n";
			return 2.0*precision;
		}
	}
	return 2.0L*precision;
}

/** Evaluates the Gamma function for integers - Factorial function */
long double gamma(unsigned int n) {
	long double fatorial = (long double)(n - 1);
	for (int i=n-2; i>1; --i)
		fatorial *= (long double)i;
	return fatorial;
}

/** Evaluates the Gamma function for a real number x, Gamma(x) */
long double gamma(long double x) {
	//  Coefficients for minimax approximation over (see reference).
	long double c[7] = { -1.910444077728E-03L, 8.4171387781295E-04L, -5.952379913043012E-04L,
		7.93650793500350248E-04L, -2.777777777777681622553E-03L, 8.333333333333333331554247E-02L,
		5.7083835261E-03L };
	long double precision = 2.22E-16L;
	long double fact, half = 0.5L;
	int i, n;
	long double p[8] = { -1.71618513886549492533811E+00L, 2.47656508055759199108314E+01L, 
		-3.79804256470945635097577E+02L, 6.29331155312818442661052E+02L, 8.66966202790413211295064E+02L, 
		-3.14512729688483675254357E+04L, -3.61444134186911729807069E+04L, 6.64561438202405440627855E+04L };
	bool parity;
	long double q[8] = { -3.08402300119738975254353E+01L, 3.15350626979604161529144E+02L, 
		-1.01515636749021914166146E+03L, -3.10777167157231109440444E+03L, 2.25381184209801510330112E+04L, 
		4.75584627752788110767815E+03L, -1.34659959864969306392456E+05L, -1.15132259675553483497211E+05L };
	long double res, xnum, y, y1, ysq, z;
	long double sqrtpi = 0.9189385332046727417803297L;
	long double sum, value, xbig = 171.624L;
	long double xden, xinf = 1.79E+308L;
	long double xminin = 2.23E-308L;

	parity = false;
	fact = 1.0L;
	n = 0;
	y = x;

	// When the argument is negative
	if(y <= 0.0) {
		y = - x;
		y1 = (long double)((int)(y));   // take the integer part
		res = y - y1;
		if(res != 0.0) {
			if( y1 != (long double)(((int)(y1*half))*2.0L)) {
				parity = true;
			}
			fact = - M_PI/sinl(M_PI*res);
			y = y + 1.0L;
		}
		else {
			res = xinf;
			value = res;
			return value;
		}
	}

	// When the argument is positive but lower than precision
	if(y < precision) {
		if(xminin <= y) res = 1.0L/y;
		else {
			res = xinf;
			value = res;
			return value;
		}
	}
	// When the argument is positive but lower than 12.0L
	else if(y < 12.0L) {
		y1 = y;
		if(y < 1.0L) {
			z = y;
			y = y + 1.0L;
		}
		else {      // note the value is not higher than 12.0L 
			n = (int)(y) - 1;
			y = y - (long double)(n);
			z = y - 1.0L;
		}
		xnum = 0.0;
		xden = 1.0L;
		for(i=0;i<8;i++) {
			xnum = (xnum + p[i])*z;
			xden = xden*z + q[i];
		}
		res = xnum/xden + 1.0L;
		if(y1 < y) {
			res = res/y1;
		}
		else if(y < y1) {
			for(i=1;i<=n;i++) {
				res = res*y;
				y = y + 1.0L;
			}
		}
	}
	// When the argument is higher than 12.0L
	else {
		if(y <= xbig) {
			ysq = y * y;
			sum = c[6];
			for(i=0;i<6;i++) {
				sum = sum / ysq + c[i];
			}
			sum = sum / y - y + sqrtpi;
			sum = sum + (y-half)*log(y);
			res = exp(sum);
		}
		else {
			res = xinf;
			value = res;
			return value;
		}
	}
	// Adjusting
	if(parity) res = - res;
	if(fact != 1.0L) res = fact/res;

	return res;
}

bool TPZGaussRule::CheckCubatureRule() {
	int i;
	TPZVec<REAL> point(3,0.0L);
	long double sum = 0.0L;

	for(i=0;i<fNumInt;i++) {
		// Checking integration point belong at the master element
		point[0] = fLocation[i];
		if(!pztopology::TPZLine::IsInParametricDomain(point))
			break;
		sum += fWeight[i];
	}
	if(i==fNumInt) {
		if(IsZero((REAL)(sum) - pztopology::TPZLine::RefElVolume())) return true;
	}
	return false;   // because any integration point is outside of the master element
}