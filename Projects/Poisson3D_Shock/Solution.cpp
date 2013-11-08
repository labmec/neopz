#include "pzreal.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzmaterial.h"
#include "pzintel.h"

#define sgn(a) (a<0) ? -1 : 1

/** PROBLEM WITH HIGH GRADIENT ON CIRCUNFERENCE  ---  DATA */
extern STATE ValueK;
extern REAL GlobalMaxError;
extern STATE F;
// Circunference with high gradient - data
extern TPZManVector<REAL,3> CCircle;
extern REAL RCircle ;

REAL VarTimesVarMinusOne(int var,int dim,const TPZVec<REAL> &x) {
	if(var < dim)
		return (x[var]*(x[var]-1.));
	return 1.;
}

// Computes ||x-x0||^2
REAL SquareNorm(int dim,const TPZVec<REAL> &x,TPZVec<REAL> &x0) {
	REAL norm = 0.0;
	for(int i=0;i<dim;i++)
		norm += (x[i] - x0[i])*(x[i] - x0[i]);
	return norm;
}
REAL ArgumentArc(int dim,REAL Radius,const TPZVec<REAL> &x,TPZVec<REAL> &x0) {
	return (F*(Radius*Radius - SquareNorm(dim,x,x0)));
}
REAL PolinomicValue(int var,int dim,const TPZVec<REAL> &x,TPZVec<REAL> &x0) {
	if(var<dim)
		return (5*x[var]*x[var]-3.*x[var]+2*x0[var]*(1.-2*x[var]));
	return 1.;
}
void ExactSolutionSphere(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol, TPZVec<STATE> &ddsol) {
	int dim = dsol.Rows();
	dsol.Zero();
	ddsol.Resize(9);
	REAL Coeff, B;
	if(dim==1)
		Coeff = -2.;
	else if(dim==2)
		Coeff = 8.;
	else
		Coeff = -32.;
	B = Coeff/M_PI;
	// Argument value (arc) to compute ArcTangent( arc )
    REAL arc = ArgumentArc(dim,RCircle,x,CCircle);
    REAL Prod, temp, den;
	REAL prodx = VarTimesVarMinusOne(0,dim,x);
	REAL prody = VarTimesVarMinusOne(1,dim,x);
	REAL prodz = VarTimesVarMinusOne(2,dim,x);

    Prod = prodx*prody*prodz;
    temp = M_PI + 2.*atan(arc);
    sol[0] = B*Prod*temp;
	den = 1. + arc*arc;
	
	// Derivatives of the first order
    dsol(0,0) = B*prody*prodz*((2*x[0]-1.)*temp - (4*F*(x[0]-CCircle[0])*(prodx/den)));
    if(dim==2) {
        dsol(1,0) = B*prodx*prodz*((2*x[1]-1.)*temp - (4*F*(x[1]-CCircle[1])*(prody/den)));
    }
    else if(dim==3) {
        dsol(2,0) = B*prodx*prody*((2*x[2]-1.)*temp - (4*F*(x[2]-CCircle[2])*(prodz/den)));
    }

	//Derivatives of the second order
	// ddsol = {d2u/dx2, d2u/dxdy, d2u/dxdz, d2u/dy2, d2u/dydx, d2u/dydz, d2u/dz2, d2u/dzdx, d2u/dzdy}
	for(int i=0;i<9;i++)
		ddsol[i] = 0.0;
	REAL poli = PolinomicValue(0,dim,x,CCircle);
	REAL sum2 = 4*F*(poli/den);
	REAL sum3 = 16*F*F*prodx*(x[0]-CCircle[0])*(x[0]-CCircle[0])*(arc/(den*den));
	ddsol[0] = B*prody*prodz*(2*temp-sum2-sum3);
	if(dim > 1) {
		poli = PolinomicValue(1,dim,x,CCircle);
		sum2 = 4*F*(poli/(den));
		sum3 = 16*F*F*prody*(x[1]-CCircle[1])*(x[1]-CCircle[1])*(arc/(den*den));
		ddsol[3] = B*prodx*prodz*(2*temp-sum2-sum3);
		poli = (2*x[0] - 1)*(2*x[1] - 1.)*temp;
		sum2 = prodx*(2*x[1]-1)*(x[0]-CCircle[0])+prody*(2*x[0]-1)*(x[1]-CCircle[1]);
		sum3 = prodx*prody*(x[0]-CCircle[0])*(x[1]-CCircle[1])*arc;
		ddsol[4] = B*prodz*(poli-((4*F*sum2)/den)-((16*F*F*sum3)/(den*den)));
		ddsol[1] = ddsol[4];
	}
	if(dim==3) {
		poli = PolinomicValue(2,dim,x,CCircle);
		sum2 = 4*F*(poli/(den));
		sum3 = 16*F*F*prodz*(x[2]-CCircle[2])*(x[2]-CCircle[2])*(arc/(den*den));
		ddsol[6] = B*prodx*prody*(2*temp-sum2-sum3);
		poli = (2*x[0] - 1)*(2*x[2] - 1.)*temp;
		sum2 = prodx*(2*x[2]-1)*(x[0]-CCircle[0])+prodz*(2*x[0]-1)*(x[2]-CCircle[2]);
		sum3 = prodx*prodz*(x[0]-CCircle[0])*(x[2]-CCircle[2])*arc;
		ddsol[2] = B*prody*(poli-((4*F*sum2)/den)-((16*F*F*sum3)/(den*den)));
		ddsol[7] = ddsol[2];
		poli = (2*x[1] - 1)*(2*x[2] - 1.)*temp;
		sum2 = prody*(2*x[2]-1)*(x[1]-CCircle[1])+prodz*(2*x[1]-1)*(x[2]-CCircle[2]);
		sum3 = prody*prodz*(x[1]-CCircle[1])*(x[2]-CCircle[2])*arc;
		ddsol[5] = B*prodx*(poli-((4*F*sum2)/den)-((16*F*F*sum3)/(den*den)));
		ddsol[8] = ddsol[5];
	}
}
void ExactSolutionSphere(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) {
	int dim = dsol.Rows();
	dsol.Zero();
	REAL Coeff, B;
	if(dim==1)
		Coeff = -2.;
	else if(dim==2)
		Coeff = 8.;
	else
		Coeff = -32.;
	B = Coeff/M_PI;
	// Argument value (arc) to compute ArcTangent( arc )
    REAL arc = ArgumentArc(dim,RCircle,x,CCircle);
    REAL Prod, temp, den;
	REAL prodx = VarTimesVarMinusOne(0,dim,x);
	REAL prody = VarTimesVarMinusOne(1,dim,x);
	REAL prodz = VarTimesVarMinusOne(2,dim,x);

    Prod = prodx*prody*prodz;
    temp = M_PI + 2.*atan(arc);
    sol[0] = B*Prod*temp;
	den = 1. + arc*arc;
	
	// Derivatives of the first order
    dsol(0,0) = B*prody*prodz*((2*x[0]-1.)*temp - (4*F*(x[0]-CCircle[0])*(prodx/den)));
    if(dim==2) {
        dsol(1,0) = B*prodx*prodz*((2*x[1]-1.)*temp - (4*F*(x[1]-CCircle[1])*(prody/den)));
    }
    else if(dim==3) {
        dsol(2,0) = B*prodx*prody*((2*x[2]-1.)*temp - (4*F*(x[2]-CCircle[2])*(prodz/den)));
    }
}

bool GradientAndLaplacianOnCorners(TPZInterpolatedElement *el,REAL &Grad,REAL &Laplacian) {
	Grad = 0.0;
	Laplacian = 0.0;
	if(!el) return false;
	int nstates = el->Material()->NStateVariables();
    int dim = el->Dimension(), i, idsol;
    TPZManVector<STATE,3> sol(nstates,(STATE)0.0);
    TPZFMatrix<STATE> dsol(dim,nstates);
	dsol.Zero();
	TPZVec<REAL> qsi(3,0.0);
	TPZManVector<REAL,3> x(3,0.0);
	TPZManVector<STATE,9> deriv2(9,0.0);

	REAL GradTemp, LaplacianTemp;
	TPZManVector<STATE,3> GradUSquare(3,0.0);

	int ncorners = el->NCornerConnects();
	// Computing on all corners of the element
	for(i=0;i<ncorners+1;i++) {
		// After computing over corners then must to be computed on the center of the element
		if(i==ncorners)
			i=el->NConnects()-1;
		// Computing gradient and gradient norm maxime
	    GradTemp = 0.0;
		LaplacianTemp = 0.0;
		el->Reference()->CenterPoint(i, qsi);
		el->Reference()->X(qsi,x);
		ExactSolutionSphere(x,sol,dsol,deriv2);
		for(idsol=0;idsol<dim;idsol++) {
			GradUSquare[idsol] = dsol.GetVal(idsol,0)*dsol.GetVal(idsol,0);
			GradTemp += GradUSquare[idsol];
		}
		GradTemp = sqrt(GradTemp);
		Grad = (Grad > GradTemp) ? Grad : GradTemp;

		// Computing the curvature on the line with gradient vector as vector direction
		if(dim==1)
			LaplacianTemp = fabs(deriv2[0]);
		else if(dim==2) {
			LaplacianTemp = fabs(deriv2[0]*GradUSquare[0]+deriv2[3]*GradUSquare[1]+(dsol.GetVal(0,0)*dsol.GetVal(1,0)*(deriv2[1]+deriv2[4])));
		}
		else {
			LaplacianTemp = deriv2[0]*GradUSquare[0]+deriv2[3]*GradUSquare[1]+deriv2[6]*GradUSquare[2];
			GradTemp = LaplacianTemp + dsol(0,0)*dsol(1,0)*(deriv2[1]+deriv2[4]);
			GradTemp += dsol(0,0)*dsol(2,0)*(deriv2[2]+deriv2[7]);
			GradTemp += dsol(1,0)*dsol(2,0)*(deriv2[5]+deriv2[8]);
			LaplacianTemp = fabs(GradTemp);
		}
		Laplacian = (Laplacian > LaplacianTemp) ? Laplacian : LaplacianTemp;
	}
	return true;
}

STATE GradientNorm(TPZInterpolatedElement *el) {
	STATE Grad = 0.0;
	if(!el) return false;
	int nstates = el->Material()->NStateVariables();
    int dim = el->Dimension(), i, idsol;
    TPZManVector<STATE,3> sol(nstates,(STATE)0.0);
    TPZFMatrix<STATE> dsol(dim,nstates);
	dsol.Zero();
	TPZVec<REAL> qsi(3,0.0);
	TPZManVector<REAL,3> x(3,0.0);
	TPZManVector<STATE,9> deriv2(9,0.0);

	REAL GradTemp;
	TPZManVector<STATE,3> GradUSquare(3,0.0);

	int ncorners = el->NCornerConnects();
	// Computing on all corners of the element
	for(i=0;i<ncorners+1;i++) {
		// After computing over corners then must to be computed on the center of the element
		if(i==ncorners)
			i=el->NConnects()-1;
		// Computing gradient and gradient norm maxime
	    GradTemp = 0.0;
		el->Reference()->CenterPoint(i, qsi);
		el->Reference()->X(qsi,x);
		ExactSolutionSphere(x,sol,dsol,deriv2);
		for(idsol=0;idsol<dim;idsol++) {
			GradUSquare[idsol] = dsol.GetVal(idsol,0)*dsol.GetVal(idsol,0);
			GradTemp += GradUSquare[idsol];
		}
		GradTemp = sqrt(GradTemp);
		Grad = (Grad > GradTemp) ? Grad : GradTemp;

	}
	return Grad;
}
STATE Laplacian(TPZInterpolatedElement *el) {
	REAL Laplacian;
	if(!el) return false;
	int nstates = el->Material()->NStateVariables();
    int dim = el->Dimension(), i, idsol;
    TPZManVector<STATE,3> sol(nstates,(STATE)0.0);
    TPZFMatrix<STATE> dsol(dim,nstates);
	dsol.Zero();
	TPZVec<REAL> qsi(3,0.0);
	TPZManVector<REAL,3> x(3,0.0);
	TPZManVector<STATE,9> deriv2(9,0.0);

	REAL LaplacianTemp, Temp;
	TPZManVector<STATE,3> GradUSquare(3,0.0);

	int ncorners = el->NCornerConnects();
	// Computing on all corners of the element
	for(i=0;i<ncorners+1;i++) {
		// After computing over corners then must to be computed on the center of the element
		if(i==ncorners)
			i=el->NConnects()-1;
		// Computing gradient and gradient norm maxime
		LaplacianTemp = 0.0;
		el->Reference()->CenterPoint(i, qsi);
		el->Reference()->X(qsi,x);
		ExactSolutionSphere(x,sol,dsol,deriv2);
		for(idsol=0;idsol<dim;idsol++) {
			GradUSquare[idsol] = dsol.GetVal(idsol,0)*dsol.GetVal(idsol,0);
		}

		// Computing the curvature on the line with gradient vector as vector direction
		if(dim==1)
			LaplacianTemp = fabs(deriv2[0]);
		else if(dim==2) {
			LaplacianTemp = fabs(deriv2[0]*GradUSquare[0]+deriv2[3]*GradUSquare[1]+(dsol.GetVal(0,0)*dsol.GetVal(1,0)*(deriv2[1]+deriv2[4])));
		}
		else {
			LaplacianTemp = deriv2[0]*GradUSquare[0]+deriv2[3]*GradUSquare[1]+deriv2[6]*GradUSquare[2];
			Temp = LaplacianTemp + dsol(0,0)*dsol(1,0)*(deriv2[1]+deriv2[4]);
			Temp += dsol(0,0)*dsol(2,0)*(deriv2[2]+deriv2[7]);
			Temp += dsol(1,0)*dsol(2,0)*(deriv2[5]+deriv2[8]);
			LaplacianTemp = fabs(Temp);
		}
		Laplacian = (Laplacian > LaplacianTemp) ? Laplacian : LaplacianTemp;
	}
	return Laplacian;
}

/** We are considering - f, because is as TPZMatPoisson3d was implemented in Contribute method */
void RightTermCircle(const TPZVec<REAL> &x, TPZVec<STATE> &force, TPZFMatrix<STATE> &dforce) {
	int dim = dforce.Rows();
    TPZManVector<STATE,3> sol(1,(STATE)0.0);
    TPZFMatrix<STATE> dsol(dim,1);
	dsol.Zero();
	TPZVec<REAL> qsi(3,0.0);
	TPZManVector<STATE,9> deriv2(9,0.0);

	// Computing gradient and gradient norm maxime
	ExactSolutionSphere(x,sol,dsol,deriv2);
	force[0] = 0.0;
	for(int i=0;i<dim;i++)
		force[0] += deriv2[3*i];
	force[0] *= (ValueK);
/*	REAL B = ValueK/M_PI;
	REAL F = 2*sqrt(ValueK);
	REAL Coeff;
	if(dim==1)
		 Coeff = -2.;
	else if(dim==2)
		Coeff = 8.;
	else
		Coeff = -32.;
	B *= (2.*Coeff);

	REAL Prod, prodx, prody, prodz;
    REAL Soma, arc, den;
    prodx = x[0]*(x[0] - 1.);
	if(dim==1) {
		arc = F*(RCircle*RCircle - (x[0]-CCircle[0])*(x[0]-CCircle[0]));
        prody = prodz = 1.;
        Soma = 1.;
	}
	else if(dim == 2) {
		arc = F*(RCircle*RCircle - ((x[0]-CCircle[0])*(x[0]-CCircle[0]) + (x[1]-CCircle[1])*(x[1]-CCircle[1])));
        prody = x[1]*(x[1] - 1.);
        prodz = 1.;
        Soma = prodx + prody;
	}
	else {
		arc = F*(RCircle*RCircle - ((x[0]-CCircle[0])*(x[0]-CCircle[0]) + (x[1]-CCircle[1])*(x[1]-CCircle[1]) + (x[2]-CCircle[2])*(x[2]-CCircle[2])));
		prody = x[1]*(x[1] - 1.);
		prodz = x[2]*(x[2] - 1.);
        Soma = prodx*prody + prodx*prodz + prody*prodz;
	}
    Prod = prodx*prody*prodz;
	den = 1. + arc*arc;
	force[0] = Soma*(M_PI + 2.*atan(arc));
	force[0] += (-2.)*(F/den)*(5*dim*Prod+Soma);
	force[0] += (F*Prod*arc*(8.*arc-0.5*F)/(den*den));
	force[0] *= B;
	*/
}
