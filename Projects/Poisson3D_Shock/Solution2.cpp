#include "pzreal.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzmaterial.h"
#include "pzintel.h"

#include "problem.h"



REAL VarTimesOneMinusVar(int var,int dim,const TPZVec<REAL> &x) {
	if(var < dim)
		return (x[var]*(1. - x[var]));
	return 1.;
}

// Computes ||x-x0||^2
REAL SquareNorm(int dim,const TPZManVector<REAL> &x,TPZVec<REAL> &x0) {
	REAL norm = 0.0;
	int i;
	for(i=0;i<dim;i++)
		norm += (x[i] - x0[i])*(x[i] - x0[i]);
	return norm;
}
REAL ArgumentArc(int dim,REAL Radius,const TPZVec<REAL> &x,TPZVec<REAL> &x0) {
	return (F*(Radius*Radius - SquareNorm(dim,x,x0)));
}
REAL PolinomicValue(int var,int dim,const TPZVec<REAL> &x,TPZVec<REAL> &x0) {
	if(var<dim)
		return (5*x[var]*x[var]-4.*x[var]*x0[var]-3.*x[var]+2.*x0[var]);
	return 1.;
}
void ExactSolutionArcTangent(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol, TPZVec<STATE> &ddsol) {
	ExactSolutionArcTangent(x,sol,dsol);
	REAL B = 5.;
	if(ModelDimension==1)
		B *= 0.25;
	else if(ModelDimension==3)
		B *= 4;
    REAL arc = ArgumentArc(ModelDimension,RCircle,x,CCircle);
    REAL temp, den;
    temp = 0.5*M_PI + atan(arc);
	den = 1. + arc*arc;
	REAL prodx = VarTimesOneMinusVar(0,ModelDimension,x);
	REAL prody = VarTimesOneMinusVar(1,ModelDimension,x);
	REAL prodz = VarTimesOneMinusVar(2,ModelDimension,x);
	//Derivatives of the second order
	// ddsol = {d2u/dx2, d2u/dxdy, d2u/dxdz, d2u/dy2, d2u/dydx, d2u/dydz, d2u/dz2, d2u/dzdx, d2u/dzdy}
	for(int i=0;i<9;i++)
		ddsol[i] = 0.0;
	REAL poli = PolinomicValue(0,ModelDimension,x,CCircle);
	REAL sum2 = F*(poli/den);
	REAL sum3 = 4.*F*F*prodx*(x[0]-CCircle[0])*(x[0]-CCircle[0])*(arc/(den*den));
	ddsol[0] = -2.*B*prody*prodz*(temp-sum2+sum3);
	if(ModelDimension > 1) {
		poli = PolinomicValue(1,ModelDimension,x,CCircle);
		sum2 = F*(poli/(den));
		sum3 = 4.*F*F*prody*(x[1]-CCircle[1])*(x[1]-CCircle[1])*(arc/(den*den));
		ddsol[3] = -2.*B*prodx*prodz*(temp-sum2+sum3);
		poli = (1.-2*x[0])*(1.-2*x[1])*temp;
		sum2 = prodx*(1.-2*x[1])*(x[0]-CCircle[0])+prody*(1.-2*x[0])*(x[1]-CCircle[1]);
		sum3 = prodx*prody*(x[0]-CCircle[0])*(x[1]-CCircle[1])*arc;
		ddsol[4] = B*prodz*(poli-((2.*F*sum2)/den)-((8.*F*F*sum3)/(den*den)));
		ddsol[1] = ddsol[4];
	}
	if(ModelDimension==3) {  // Revisar
		poli = PolinomicValue(2,ModelDimension,x,CCircle);
		sum2 = F*(poli/(den));
		sum3 = 4.*F*F*prodz*(x[2]-CCircle[2])*(x[2]-CCircle[2])*(arc/(den*den));
		ddsol[6] = -2.*B*prodx*prody*(temp-sum2+sum3);
		poli = (1.-2*x[0])*(1.-2*x[2])*temp;
		sum2 = prodx*(1.-2*x[2])*(x[0]-CCircle[0])+prodz*(1.-2*x[0])*(x[2]-CCircle[2]);
		sum3 = prodx*prodz*(x[0]-CCircle[0])*(x[2]-CCircle[2])*arc;
		ddsol[7] = B*prody*(poli-((2.*F*sum2)/den)-((8.*F*F*sum3)/(den*den)));
		ddsol[2] = ddsol[7];
		poli = (1.-2*x[1])*(1.-2*x[2])*temp;
		sum2 = prody*(1.-2*x[2])*(x[1]-CCircle[1])+prodz*(1.-2*x[1])*(x[2]-CCircle[2]);
		sum3 = prody*prodz*(x[1]-CCircle[1])*(x[2]-CCircle[2])*arc;
		ddsol[8] = B*prodx*(poli-((2.*F*sum2)/den)-((8.*F*F*sum3)/(den*den)));
		ddsol[5] = ddsol[8];
	}
}
void ExactSolutionArcTangent(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) {
	dsol.Zero();
	REAL B = 5.;
	if(ModelDimension==1)
		B *= 0.25;
	else if(ModelDimension==3)
		B *= 4;
	// Argument value (arc) to compute ArcTangent( arc )
    REAL arc = ArgumentArc(ModelDimension,RCircle,x,CCircle);
    REAL Prod, temp, den;
	REAL prodx = VarTimesOneMinusVar(0,ModelDimension,x);
	REAL prody = VarTimesOneMinusVar(1,ModelDimension,x);
	REAL prodz = VarTimesOneMinusVar(2,ModelDimension,x);

    Prod = prodx*prody*prodz;
    temp = 0.5*M_PI + atan(arc);
    sol[0] = B*Prod*temp;
	den = 1. + arc*arc;
	
	// Derivatives of the first order
    dsol(0,0) = B*prody*prodz*((1.-2*x[0])*temp - (2*F*(x[0]-CCircle[0])*(prodx/den)));
    if(ModelDimension>1) {
        dsol(1,0) = B*prodx*prodz*((1.-2*x[1])*temp - (2*F*(x[1]-CCircle[1])*(prody/den)));
    }
    else if(ModelDimension>2) {
        dsol(2,0) = B*prodx*prody*((1.-2*x[2])*temp - (2*F*(x[2]-CCircle[2])*(prodz/den)));
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
		ExactSolutionArcTangent(x,sol,dsol,deriv2);
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
bool GradientAndLaplacian(TPZInterpolatedElement *el,REAL &Grad,REAL &Laplacian) {
	Grad = 0.0;
	Laplacian = 0.0;
	if(!el) return false;
	int nstates = el->Material()->NStateVariables();
    int dim = el->Dimension(), idsol;
    TPZManVector<STATE,3> sol(nstates,(STATE)0.0);
    TPZFMatrix<STATE> dsol(dim,nstates);
	dsol.Zero();
	TPZVec<REAL> qsi(3,0.0);
	TPZManVector<REAL,3> x(3,0.0);
	TPZManVector<STATE,9> deriv2(9,0.0);

	TPZManVector<STATE,3> GradUSquare(3,0.0);

	// Computing on the center of the element
	// Computing gradient and gradient norm maxime
	el->Reference()->CenterPoint(el->NConnects()-1, qsi);
	el->Reference()->X(qsi,x);
	ExactSolutionArcTangent(x,sol,dsol,deriv2);
	for(idsol=0;idsol<dim;idsol++) {
		GradUSquare[idsol] = dsol.GetVal(idsol,0)*dsol.GetVal(idsol,0);
		Grad += GradUSquare[idsol];
	}
	Grad = sqrt(Grad);

	// Computing the curvature on the line with gradient vector as vector direction
	if(dim==1)
		Laplacian = fabs(deriv2[0]);
	else if(dim==2) {
		Laplacian = fabs(deriv2[0]*GradUSquare[0]+deriv2[3]*GradUSquare[1]+(dsol.GetVal(0,0)*dsol.GetVal(1,0)*(deriv2[1]+deriv2[4])));
	}
	else {
		Laplacian = deriv2[0]*GradUSquare[0]+deriv2[3]*GradUSquare[1]+deriv2[6]*GradUSquare[2];
		REAL Temp = Laplacian + dsol(0,0)*dsol(1,0)*(deriv2[1]+deriv2[4]);
		Temp += dsol(0,0)*dsol(2,0)*(deriv2[2]+deriv2[7]);
		Temp += dsol(1,0)*dsol(2,0)*(deriv2[5]+deriv2[8]);
		Laplacian = fabs(Temp);
	}
	return true;
}
bool LaplacianValue(TPZInterpolatedElement *el,REAL &Laplacian) {
	Laplacian = 0.0;
	if(!el) return false;
	int nstates = el->Material()->NStateVariables();
    int dim = el->Dimension();
    TPZManVector<STATE,3> sol(nstates,(STATE)0.0);
    TPZFMatrix<STATE> dsol(dim,nstates);
	dsol.Zero();
	TPZVec<REAL> qsi(3,0.0);
	TPZManVector<REAL,3> x(3,0.0);
	TPZManVector<STATE,9> deriv2(9,0.0);
    TPZVec<REAL> force(1,0.0);
    
	// Computing on the center of the element
	// Computing gradient and gradient norm maxime
	el->Reference()->CenterPoint(el->NConnects()-1, qsi);
	el->Reference()->X(qsi,x);
    
    RightTermArcTangent(x,force,dsol);
    Laplacian = fabs(force[0]/ValueK);
	return true;
}

void ComputingMaxLaplacian(TPZCompMesh *cmesh,REAL &MaxLaplacian,REAL &MinLaplacian) {
	MaxLaplacian = 0.0;
    MinLaplacian = 1.e+5;
	REAL Laplace;
	long nels = cmesh->NElements();
	TPZInterpolatedElement *el;
    
	for(long i=0L;i<nels;i++) {
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el || el->Dimension()!=cmesh->Dimension()) continue;
		// If error is small and laplacian value is very little then the order will be minimized
		if(!LaplacianValue(el,Laplace))
			DebugStop();
        
        MinLaplacian = (Laplace < MinLaplacian) ? Laplace : MinLaplacian;
		MaxLaplacian = (Laplace > MaxLaplacian) ? Laplace : MaxLaplacian;
	}
}

void ComputingMaxGradientAndLaplacian(TPZCompMesh *cmesh,REAL &MaxGrad,REAL &MaxLaplacian) {
	MaxGrad = 0.0;
	MaxLaplacian = 0.0;
	REAL Grad, Laplacian;
	long nels = cmesh->NElements();
	TPZInterpolatedElement *el;

	for(long i=0L;i<nels;i++) {
		el = dynamic_cast<TPZInterpolatedElement* >(cmesh->ElementVec()[i]);
		if(!el || el->Dimension()!=cmesh->Dimension()) continue;
		// If error is small and laplacian value is very little then the order will be minimized
		if(!GradientAndLaplacian(el,Grad,Laplacian))
			DebugStop();

		MaxGrad = (Grad > MaxGrad) ? Grad : MaxGrad;
		MaxLaplacian = (Laplacian > MaxLaplacian) ? Laplacian : MaxLaplacian;
	}
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
		ExactSolutionArcTangent(x,sol,dsol,deriv2);
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
	REAL Laplacian = 0.0;
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
		ExactSolutionArcTangent(x,sol,dsol,deriv2);
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
void RightTermArcTangent(const TPZVec<REAL> &x, TPZVec<STATE> &force, TPZFMatrix<STATE> &dforce) {
    TPZManVector<STATE,3> sol(1,(STATE)0.0);
    TPZFMatrix<STATE> dsol(ModelDimension,1);
	dsol.Zero();
	TPZVec<REAL> qsi(3,0.0);
	TPZManVector<STATE,9> deriv2(9,0.0);

	// Computing gradient and gradient norm maxime
	ExactSolutionArcTangent(x,sol,dsol,deriv2);
	force[0] = 0.0;
	for(int i=0;i<ModelDimension;i++)
		force[0] += deriv2[3*i];
	force[0] *= (ValueK);
}
