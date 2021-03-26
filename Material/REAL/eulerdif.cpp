/**
 * \file
 * @brief Contains implementations of the methods related with numerical diffusivity coefficient for SUPG.
 */
#include "eulerdif.h" 
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzreal.h"
using namespace std;

STATE TEulerDiffusivity::fGamma = 1.4;


STATE TEulerDiffusivity::Pressure(TPZVec<STATE> &U) {
	if(U[0] > -1.e-6 && U[0] < 1.e-6) return (fGamma-1.)*U[3];
	STATE velocity = (U[1]*U[1]+U[2]*U[2])/U[0];
	return (fGamma-1.)*(U[3]-(.5*velocity));
}

void TEulerDiffusivity::Flux(TPZVec<STATE> &U,TPZVec<STATE> &flux) {
	//Variavel : u = (densidade, densidade*velocidade, energia)
	//Funcao : f(u)=(densid*veloc, densid*(velocid)^2+pressao, veloc*(energia+pressao))
	//Funcao : f(u)=( u2, ((u2*u2)/u1)+pressao, (u2/u1)*(u3+pressao))
	
	STATE pressao = Pressure(U);
	if(U[0] > -1.e-6 && U[0] < 1.e-6) {
		flux[0] = flux[2] = flux[3] = flux[4] = flux[5] = flux[7] = 0.;
		cout << "EulerLaw2D::Flux find nule density.\n";
		flux[1] = flux[6] = pressao;  // Jorge 17/08
		return;
	}
	
	flux[0] = U[1];
	flux[1] = ((U[1]/U[0])*U[1]) + pressao;
	flux[2] = U[1]*(U[2]/U[0]);
	flux[3] = (U[1]/U[0])*(U[3]+pressao);
	flux[4] = U[2];
	flux[5] = U[2]*(U[1]/U[0]);
	flux[6] = ((U[2]/U[0])*U[2]) + pressao;
	flux[7] = (U[2]/U[0])*(U[3] + pressao);
}

void TEulerDiffusivity::JacobFlux(TPZVec<STATE> &U,TPZFMatrix<STATE> &Ajacob,TPZFMatrix<STATE> &Bjacob) {
	if(U[0] > -1.e-6 && U[0] < 1.e-6) {
		cout << "\nERRO : Densidade NULA. Falla jacobiano.\n";
		exit(1);
	}
	
	STATE velx, vely, ener, aux, gammam1 = fGamma-1.;
	velx = U[1]/U[0];
	vely = U[2]/U[0];
	ener = U[3]/U[0];
	
	Ajacob(0,0)=Ajacob(0,2)=Ajacob(0,3)=Ajacob(2,3)=0.;
	Bjacob(0,0)=Bjacob(0,1)=Bjacob(0,3)=Bjacob(1,3)=0.;
	Ajacob(0,1)=Bjacob(0,2)=1.;
	Ajacob(1,3)=Bjacob(2,3)=gammam1;
	Ajacob(1,1)=(3.-fGamma)*velx;
	Bjacob(2,2)=(3.-fGamma)*vely;
	Ajacob(1,2)=-1.*gammam1*vely;
	Bjacob(2,1)=-1.*gammam1*velx;
	Ajacob(3,2)=Bjacob(3,1)=Ajacob(1,2)*velx;
	Ajacob(2,0)=Bjacob(1,0)=-1.*velx*vely;
	Ajacob(2,1)=Bjacob(1,1)=vely;
	Ajacob(2,2)=Bjacob(1,2)=velx;
	Ajacob(3,3)=fGamma*velx;
	Bjacob(3,3)=fGamma*vely;
	aux = gammam1*(velx*velx+vely*vely);
	Ajacob(1,0) = 0.5*aux - (velx*velx);
	Ajacob(3,1) = (fGamma*ener) - (.5*gammam1*(3*velx*velx+vely*vely));
	Bjacob(3,2) = (fGamma*ener) - (.5*gammam1*(velx*velx+3*vely*vely));
	Bjacob(2,0) = 0.5*aux - (vely*vely);
	Ajacob(3,0) = (aux-(fGamma*ener))*velx;
	Bjacob(3,0) = (aux-(fGamma*ener))*vely;
}

void TEulerDiffusivity::ValJacobFlux(TPZVec<STATE> &U,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) {
	valjacob.Zero();
	STATE rho = U[0];
	STATE vx, vy, /*lambda,*/ cval, cval2, u, kt, st, aa, bb;
	STATE gM1 = fGamma-1.;
	vx = U[1]/rho;
	vy = U[2]/rho;
	cval2 = (fGamma*Pressure(U))/rho;
	if(cval2<0) {
		cout << "TEulerLaw2D::ValJacobFlux. velocidade do som negativa.\n";
		return;
	}
	cval = sqrt(cval2);
	//  d = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
	// if(IsZero(d)) return 2;
	
	u = sqrt(vx*vx + vy*vy);
	STATE Mach = u/cval;
	STATE Mach2 = u*u/cval2;
	if(u>-1.e-6 && u< 1.e-6) {
		kt=1.;
		st=0.;
	} else {
		kt = vx/u;
		st = vy/u;
	}
	aa = normal[0]*kt + normal[1]*st;
	bb = normal[1]*kt - normal[0]*st;
	/*	
	 cout << "normal = ("<<normal[0]<<","<<normal[1]<<")" << endl;
	 cout << "aa = " << aa << "  bb = " << bb << endl;
	 cout << "kt = " << kt << "  st = " << st << endl;
	 cout << "cval2 = " << cval2 << "  Mach = " << Mach << endl;
	 */
	STATE aabbnorm2 = aa*aa+bb*bb;
	STATE aabbnorm = sqrt(aabbnorm2);
	//Para o computo da matrix direita
	TPZFMatrix<STATE> right(4,4);
	TPZVec<STATE> diag(4);
	STATE aascal = aa/aabbnorm;
	STATE bbscal = bb/aabbnorm;
	right(0,0) = -(bbscal*cval*Mach);
	right(0,1) = bbscal*kt + aascal*st;
	right(0,2) = -(aascal*kt) + bbscal*st;
	right(0,3) = 0.;
	
	right(1,0) = 1 - gM1*Mach2/2.;
	right(1,1) = (gM1*kt*Mach)/cval;
	right(1,2) = (gM1*Mach*st)/cval;
	right(1,3) = -gM1/cval2;
	
	right(2,0) = (cval2*Mach*(2*aascal + gM1*Mach))/4.;
	right(2,1) = -(cval*(aascal*kt + gM1*kt*Mach - bbscal*st))/2.;
	right(2,2) = -(cval*(bbscal*kt + (aascal + gM1*Mach)*st))/2.;
	right(2,3) = gM1/2.;
	
	right(3,0) = (cval2*Mach*(-2*aascal + gM1*Mach))/4.;
	right(3,1) = (cval*(aascal*kt - gM1*kt*Mach - bbscal*st))/2.;
	right(3,2) = (cval*(bbscal*kt + (aascal - gM1*Mach)*st))/2.;
	right(3,3) = gM1/2.;
	
	//matriz esquerda RtQX*val(autovalores)
	TPZFMatrix<STATE> left(4,4);
	left(0,0) = 0.;
	left(0,1) = 1.;
	left(0,2) = 1/cval2;
	left(0,3) = 1/cval2;
	
	left(1,0) = bbscal*kt + aascal*st;
	left(1,1) = cval*kt*Mach;
	left(1,2) = (-(aascal*kt) + kt*Mach + bbscal*st)/cval;
	left(1,3) = (aascal*kt + kt*Mach - bbscal*st)/cval;
	
	left(2,0) = -(aascal*kt) + bbscal*st;
	left(2,1) = cval*Mach*st;
	left(2,2) = -((bbscal*kt + aascal*st - Mach*st)/cval);
	left(2,3) = (bbscal*kt + (aascal + Mach)*st)/cval;
	
	
	left(3,0) = bbscal*cval*Mach;
	left(3,1) = (cval2*Mach2)/2.;
	left(3,2) = 1/gM1 + (Mach*(-2*aascal + Mach))/2.;
	left(3,3) = 1/gM1 + (Mach*(2*aascal + Mach))/2.;
	
	//right.Print("right eigenvalues");
	//left.Print("left eigenvalues");
	diag[0] = fabs(aa*cval*Mach);
	diag[1] = fabs(aa*cval*Mach);
	diag[2] = fabs(-(aabbnorm*cval) + aa*cval*Mach);
	diag[3] = fabs(cval*(aabbnorm + aa*Mach));
	
	int i,j,k;
	for(i=0;i<4;i++) {
		for(j=0;j<4;j++) {
			valjacob(i,j)=0.;
			for(k=0;k<4;k++) {
				valjacob(i,j) += left(i,k)*diag[k]*right(k,j);
			}
		}
	}
	//valjacob.Print("valjacob");
	return;
}

//Dada uma matrix da forma {{0,a,b,x},{c,d,e,f},{g,h,i,j},{k,l,m,n}} retorna a inversa dela
void TEulerDiffusivity::InverseJacob(TPZFMatrix<STATE> &jac) {
	STATE a=jac.GetVal(0,0), b=jac.GetVal(0,1), c=jac.GetVal(0,2), d=jac.GetVal(0,3);
	STATE e=jac.GetVal(1,0), f=jac.GetVal(1,1), g=jac.GetVal(1,2), h=jac.GetVal(1,3);
	STATE i=jac.GetVal(2,0), j=jac.GetVal(2,1), k=jac.GetVal(2,2), l=jac.GetVal(2,3);
	STATE m=jac.GetVal(3,0), n=jac.GetVal(3,1), r=jac.GetVal(3,2), s=jac.GetVal(3,3);
	STATE BHDF = b*h-d*f, BKCJ = b*k-c*j, BSDN = b*s-d*n, CFBG = c*f-b*g;
	STATE CLDK = c*l-d*k, CNBR = c*n-b*r, DGCH = d*g-c*h, DJBL = d*j-b*l;
	STATE DRCS = d*r-c*s, FRGN = f*r-g*n, GJFK = g*j-f*k, GSHR = g*s-h*r;
	STATE FSHN = f*s-h*n, LNJS = l*n-j*s, HJFL = h*j-f*l, HKGL = h*k-g*l;
	STATE KNJR = k*n-j*r, LRKS = l*r-k*s;
	STATE det = DJBL*(g*m-e*r)+BHDF*(k*m-i*r)+BSDN*(g*i-e*k);
	det += (FSHN*(a*k-c*i)+LNJS*(a*g-c*e)+HJFL*(a*r-c*m));
	if(det > -1.e-6 && det < 1.e-6) {
		cout << "TEulerLaw2D::InverseJacob. Determinante zero, matriz singular.\n";
		exit(1);
	}
	
	det = 1./det;
	jac(0,0) = det*(k*FSHN+g*LNJS+r*HJFL);
	jac(1,0) = det*(m*HKGL+i*GSHR+e*LRKS);
	jac(2,0) = (-det)*(m*HJFL+i*FSHN+e*LNJS);
	jac(3,0) = det*(m*GJFK+i*FRGN+e*KNJR);
	jac(0,1) = det*(d*KNJR+b*LRKS-c*LNJS);
	jac(1,1) = det*(m*CLDK+i*DRCS-a*LRKS);
	jac(2,1) = det*(m*DJBL+i*BSDN+a*LNJS);
	jac(3,1) = det*(m*BKCJ+i*CNBR-a*KNJR);
	jac(0,2) = det*(d*FRGN-c*FSHN+b*GSHR);
	jac(1,2) = det*(m*DGCH-e*DRCS-a*GSHR);
	jac(2,2) = det*(m*BHDF-e*BSDN+a*FSHN);
	jac(3,2) = det*(m*CFBG-e*CNBR-a*FRGN);
	jac(0,3) = det*(d*GJFK-c*HJFL+b*HKGL);
	jac(1,3) = (-det)*(i*DGCH+e*CLDK+a*HKGL);
	jac(2,3) = det*(a*HJFL-i*BHDF-e*DJBL);
	jac(3,3) = (-det)*(i*CFBG+e*BKCJ+a*GJFK);
}

void TEulerDiffusivity::InvJacob2d(TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &jacinv) {
	
	REAL tmp[2][2];
	// d(ksi) / dx
	tmp[0][0] = jacinv(0,0)*axes(0,0)+jacinv(0,1)*axes(1,0);
	// d(ksi) / dy
	tmp[0][1] = jacinv(0,0)*axes(0,1)+jacinv(0,1)*axes(1,1);
	// d(eta) / dx
	tmp[1][0] = jacinv(1,0)*axes(0,0)+jacinv(1,1)*axes(1,0);
	// d(eta) / dy
	tmp[1][1] = jacinv(1,0)*axes(0,1)+jacinv(1,1)*axes(1,1);
	
	jacinv(0,0) = tmp[0][0];
	jacinv(0,1) = tmp[0][1];
	jacinv(1,0) = tmp[1][0];
	jacinv(1,1) = tmp[1][1];
}
	
void TEulerDiffusivity::MatrixDiff(TPZVec<STATE> &sol,TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &jacinv,
								   TPZFMatrix<STATE> &ATauA,TPZFMatrix<STATE> &ATauB,TPZFMatrix<STATE> &BTauA,
								   TPZFMatrix<STATE> &BTauB) {
	
	//Computando o jacobiano da transformacao do elemento mestre ao elemento triangular
	InvJacob2d(axes,jacinv);
	
	//int dsolr = dsol.Rows(), dsolc = dsol.Cols();
	// Vetor de coeficientes dos jacobianos (alfa,beta) onde alfa*B+beta*B
	TPZVec<REAL> Beta(2);
	
	TPZFMatrix<STATE> Ajacob(4,4);
	TPZFMatrix<STATE> Bjacob(4,4);
	TPZFMatrix<STATE> valjacob(4,4);
	TPZFMatrix<STATE> invvaljacob(4,4,0.);
	
	JacobFlux(sol,Ajacob,Bjacob);
	
	// Calculando Tau
	// Matrix |	jacinv(0,0)*A + jacinv(0,1)*B|
	Beta[0] = jacinv(0,0); Beta[1] = jacinv(0,1);
	ValJacobFlux(sol,valjacob,Beta);
	//valjacob.Print("Jacob ValAbs dksi");
	
	invvaljacob = valjacob;
	Beta[0] = jacinv(1,0); Beta[1] = jacinv(1,1);
	ValJacobFlux(sol,valjacob,Beta);
	//valjacob.Print("Jacob ValAbs deta");
	
	invvaljacob += valjacob;
	InverseJacob(invvaljacob);
	//invvaljacob.Print("InvValJacob");
	
	ATauA.Zero();
	ATauB.Zero();
	BTauA.Zero();
	BTauB.Zero();
	int i,j,k,l;
	for(i=0;i<4;i++) {
		for(j=0;j<4;j++) {
			for(k=0;k<4;k++) {
				for(l=0;l<4;l++) {
					ATauA(i,j) += Ajacob(i,k)*invvaljacob(k,l)*Ajacob(l,j);
					ATauB(i,j) += Ajacob(i,k)*invvaljacob(k,l)*Bjacob(l,j);
					BTauA(i,j) += Bjacob(i,k)*invvaljacob(k,l)*Ajacob(l,j);
					BTauB(i,j) += Bjacob(i,k)*invvaljacob(k,l)*Bjacob(l,j);
				}
			}
		}
	}
}

void TEulerDiffusivity::JacobFlux(TPZVec<STATE> &U,TPZFMatrix<STATE> &jacob,TPZVec<REAL>
								  &normal) {
	
	STATE velx, vely, ener, aux, gammam1 = fGamma-1.;
	velx = U[1]/U[0];
	vely = U[2]/U[0];
	ener = U[3]/U[0];
	STATE alfa = normal[0], beta= normal[1];
	aux = gammam1*(velx*velx+vely*vely);
	
	jacob(0,0)=jacob(0,3)=0.;
	jacob(0,1)=alfa; jacob(0,2)=beta;
	jacob(1,3)= alfa*gammam1;
	jacob(2,3)= beta*gammam1;
	jacob(1,1)= alfa*((3.-fGamma)*velx)+beta*vely;
	jacob(1,2)=-alfa*gammam1*vely+beta*velx;
	jacob(3,2)=-alfa*gammam1*velx*vely+beta*((fGamma*ener) - .5*gammam1*(3*vely*vely+velx*velx));
	jacob(2,0)=-alfa*velx*vely+beta*(0.5*aux - (vely*vely));
	jacob(2,1)=alfa*vely-beta*gammam1*velx;
	jacob(2,2)=alfa*velx+beta*(3.-fGamma)*vely;
	jacob(3,3)=fGamma*(alfa*velx+beta*vely);
	jacob(1,0) = alfa*(0.5*aux - (velx*velx))-beta*(velx*vely);
	jacob(3,1) = alfa*((fGamma*ener) - .5*gammam1*(3*velx*velx+vely*vely))-beta*gammam1*vely*velx;
	jacob(3,0) = (alfa*velx+beta*vely)*(aux - fGamma*ener);
}

int TEulerDiffusivity::main() {
	
	TEulerDiffusivity eul;
	TPZVec<STATE> U(4);
	U[0] = 5.4;
	U[1] = -.5;
	U[2] = -.25;
	U[3] = 50.8;
	TPZVec<STATE> flux(8);
	eul.Flux(U,flux);
	
	TPZVec<REAL> normal(2);
	normal[0] = 0.5;
	normal[1] = sqrt(.75);
	TPZFMatrix<STATE> A(4,4);
	TPZFMatrix<STATE> B(4,4);
	TPZFMatrix<STATE> jacob(4,4);
	eul.JacobFlux(U,A,B);
	eul.JacobFlux(U,jacob,normal);
	A *= normal[0];
	B *= normal[1];
	int i,j;
	for(j=0;j<4;j++)
		for(i=0;i<4;i++)
		    A(i,j) += B(i,j);
	
	eul.ValJacobFlux(U,jacob,normal);
	jacob.Print("ValJacob");
	TPZFMatrix<REAL> axes(3,3,0.),jacinv(2,2,0.);
    TPZFNMatrix<16,STATE> ATA(4,4), ATB(4,4), BTA(4,4),BTB(4,4);
	axes(0,0)=1.5; axes(0,1)=2.3; axes(1,0)=.1; axes(1,1)=-1.9;
	jacinv(0,0) = 1.3;jacinv(0,1) = -2.4;jacinv(1,0)=3.5;jacinv(1,1)=5.2;
	eul.MatrixDiff(U,axes,jacinv,ATA,ATB,BTA,BTB);
	ATA.Print("ATA ");
	ATB.Print("ATB ");
	BTA.Print("BTA ");
	BTB.Print("BTB ");
	return 0;
}
