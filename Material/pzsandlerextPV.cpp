//
//  pzsandlerextPV.cpp
//  PZ
//
//  Created by Diogo Cecilio on 9/3/13.
//
//

#include "pzsandlerextPV.h"



TPZSandlerExtended::TPZSandlerExtended(TPZSandlerExtended & copy)
{
    fA=copy.fA;
    fB=copy.fB;
    fC=copy.fC;
    fD=copy.fD;
    fK=copy.fK;
    fG=copy.fG;
    fW=copy.fW;
    fR=copy.fR;
    fPhi=copy.fPhi;
    fN=copy.fN;
    fPsi=copy.fPsi;
}

TPZSandlerExtended::TPZSandlerExtended(REAL A, REAL B,REAL C, REAL D,REAL K,REAL G,REAL W,REAL R,REAL Phi,REAL N,REAL Psi):
fA(A),fB(B),fC(C),fD(D),fW(W),fK(K),fG(G),fR(R),fPhi(Phi),fN(N),fPsi(Psi)
{
    
}

TPZSandlerExtended::~TPZSandlerExtended()
{
    
}

REAL TPZSandlerExtended::F(REAL x,REAL phi)
{
    return (fA - fC*exp(x *fB) -fPhi*x);
}

REAL TPZSandlerExtended::X(REAL k)
{
    return (k - fR * F(k, fPhi));
}

REAL TPZSandlerExtended::EpsEqX(REAL X)
{
   return (fW*( exp(fD*X) - 1 ));
}

REAL TPZSandlerExtended::EpsEqk(REAL k)
{
    return EpsEqX(X(k));
}

void TPZSandlerExtended::Firstk(REAL &k)
{
    REAL f,df,kn1,kn,resnorm,diff;
    int counter =1;
    resnorm=1;
    kn=0.;//chute inicial
    while (resnorm>1.e-12 && counter<30) {
        
        f=EpsEqk(kn);
        df =fD*exp(fD*(kn - (fA - fC*exp(fB*kn))*fR))*(1 + fB*fC*exp(fB*kn)*fR)*fW;
        kn1=kn-f/df;
        diff=kn1-kn;
        resnorm=sqrt(diff*diff);
        kn=kn1;
        counter++;
 
    }
    k=kn1;
}


REAL TPZSandlerExtended::ResL(TPZTensor<REAL>::TPZDecomposed pt, REAL theta,REAL beta,REAL k,REAL kprev )
{
 
  REAL I1tr=(pt.fEigenvalues[0])+(pt.fEigenvalues[1])+(pt.fEigenvalues[2]);
  REAL I1 =fR*F(k,fPhi)* cos(theta) +k;
  REAL delepsp = EpsEqk(k) - EpsEqk(kprev);
  return  (3*fK*delepsp - (I1tr - I1));
}

TPZManVector<REAL> TPZSandlerExtended::FromHWCylToPrincipal(TPZManVector<REAL> &HWCylCoords)
{
    TPZManVector<REAL,3> SigCoords(3,0.);
    SigCoords[0]=(1/sqrt(3))*HWCylCoords[0]+sqrt(2./3.)*HWCylCoords[1]*cos(HWCylCoords[2]);
    SigCoords[1]=(1/sqrt(3))*HWCylCoords[0]+sqrt(2./3.)*HWCylCoords[1]*cos(HWCylCoords[2]-(2*M_PI/3));
    SigCoords[2]=(1/sqrt(3))*HWCylCoords[0]+sqrt(2./3.)*HWCylCoords[1]*cos(HWCylCoords[2]+(2*M_PI/3));
    return SigCoords;
}

TPZManVector<REAL>  TPZSandlerExtended::FromHWCylToCart(TPZManVector<REAL> &HWCylCoords)
{
    TPZManVector<REAL,3> cart(3,0.);
    cart[0]=HWCylCoords[0];
    cart[1]=HWCylCoords[1]*cos(HWCylCoords[2]);
    cart[2]=HWCylCoords[1]*sin(HWCylCoords[2]);
    return cart;
}

TPZManVector<REAL> TPZSandlerExtended::FromPrincipalToHWCart(TPZTensor<REAL>::TPZDecomposed &PrincipalCoords)
{
    TPZManVector<REAL,3> cartsol(3,0.);
    TPZFMatrix<REAL> Rot(3,3,0.),temp(3,1,0.),cart(3,1,0.);
    temp(0,0)=PrincipalCoords.fEigenvalues[0];
    temp(1,0)=PrincipalCoords.fEigenvalues[1];
    temp(2,0)=PrincipalCoords.fEigenvalues[2];
    GetRotMatrix(Rot);
    Rot.Multiply(temp,cart);
    cartsol[0]=cart(0,0);
    cartsol[1]=cart(1,0);
    cartsol[2]=cart(2,0);
    return cartsol;
}


TPZManVector<REAL> TPZSandlerExtended::F1Cyl(REAL xi,REAL beta)
{
    REAL gamma = 0.5*(1 + (1 - sin(3*beta))/fPsi + sin(3*beta));
    REAL I1 = xi*sqrt(3);
    REAL sqrtj2 = ( F(I1, fPhi) - fN )/gamma;
    REAL rho = sqrt(2)*sqrtj2;
    TPZManVector<REAL,3> solcyl(3,0.);
    solcyl[0]=xi,
    solcyl[1]=rho;
    solcyl[2]=beta;
    return solcyl;
    
    
}

TPZManVector<REAL> TPZSandlerExtended::F2Cyl(REAL theta,REAL beta,REAL k)
{
    REAL gamma = 0.5*(1 + (1 - sin(3*beta))/fPsi + sin(3*beta));
    REAL I1 = k + fR*F(k, fPhi)*cos(theta);
    REAL sqrtj2 = (F(k, fPhi)- fN)*sin(theta)/gamma;
    REAL rho = sqrt(2)*sqrtj2;
    REAL xi=I1/sqrt(3);
    TPZManVector<REAL,3> solcyl(3,0.);
    solcyl[0]=xi,
    solcyl[1]=rho;
    solcyl[2]=beta;
    //TPZManVector<REAL> principal = FromHWCylToPrincipal(solcyl);
    return solcyl;
    
}

void TPZSandlerExtended::GetRotMatrix(TPZFMatrix<REAL> &Rot)
{
    Rot.Resize(3,3);
    Rot(0,0)=1./sqrt(3.);    Rot(0,1)=1./sqrt(3.);    Rot(0,2)=1./sqrt(3.);
    Rot(1,0)=sqrt(2./3.);    Rot(1,1)=-1./sqrt(6.);   Rot(1,2)=-1./sqrt(6.);
    Rot(2,0)=0;              Rot(2,1)=1./sqrt(2.);    Rot(2,2)=-1./sqrt(2.);
}

REAL TPZSandlerExtended::DistF1(TPZTensor<REAL>::TPZDecomposed pt,REAL xi,REAL beta)
{
    TPZFMatrix<REAL> Rot;
    GetRotMatrix(Rot);
    TPZManVector<REAL> cyl =F1Cyl(xi,beta);
    TPZManVector<REAL,3> cart = FromHWCylToCart(cyl);
    TPZManVector<REAL,3> carttrial=FromPrincipalToHWCart(pt);
    return ((1/(3*fK))*(carttrial[0]-cart[0])*(carttrial[0]-cart[0]))
    +(1/(2*fG))*((carttrial[1]-cart[1])*(carttrial[1]-cart[1])+(carttrial[2]-cart[2])*(carttrial[2]-cart[2]));
    
}

REAL TPZSandlerExtended::DistF2(TPZTensor<REAL>::TPZDecomposed pt,REAL theta,REAL beta,REAL k)
{
    TPZFMatrix<REAL> Rot;
    GetRotMatrix(Rot);
    TPZManVector<REAL> cyl =F2Cyl(theta,beta,k);
    TPZManVector<REAL,3> cart = FromHWCylToCart(cyl);
    TPZManVector<REAL,3> carttrial=FromPrincipalToHWCart(pt);
    return ((1/(3*fK))*(carttrial[0]-cart[0])*(carttrial[0]-cart[0]))
    +(1/(2*fG))*((carttrial[1]-cart[1])*(carttrial[1]-cart[1])+(carttrial[2]-cart[2])*(carttrial[2]-cart[2]));
}


TPZFMatrix<REAL> TPZSandlerExtended::DDistFunc1(TPZTensor<REAL>::TPZDecomposed pt,REAL xi,REAL beta)
{
    REAL sig1,sig2,sig3,sinbeta,cosbeta,sin3beta,cos3beta;
    sinbeta=sin(beta);
    cosbeta=cos(beta);
    cos3beta=cos(3*beta);
    sin3beta=sin(3*beta);
    sig1 = pt.fEigenvalues[0];
    sig2 = pt.fEigenvalues[1];
    sig3 = pt.fEigenvalues[2];
    REAL ddistf1dxi = (-2*(sig1/sqrt (3) + sig2/sqrt (3) + sig3/sqrt (3) - xi))/(3.*fK) +((-4*sqrt (2)*(-(sqrt (3)*fB*fC*exp (sqrt (3)*fB*xi)) -
    sqrt (3)*fPhi)*cosbeta*(sqrt (0.6666666666666666)*sig1 - sig2/sqrt (6) -sig3/sqrt (6) -(2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -
    sqrt (3)*xi*fPhi)*cosbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta) -
     (4*sqrt (2)*(-(sqrt (3)*fB*fC*exp (sqrt (3)*fB*xi)) -sqrt (3)*fPhi)*sinbeta*(sig2/sqrt (2) - sig3/sqrt (2) -(2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -sqrt (3)*xi*fPhi)*sinbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta))/(2.*fG);
    
    REAL ddistf1dbeta=(2*((2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -sqrt (3)*xi*fPhi)*(3*cos3beta - (3*cos3beta)/fPsi)*
    sinbeta)/pow (1 + (1 - sin3beta)/fPsi + sin3beta,2) - (2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -
    sqrt (3)*xi*fPhi)*cosbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta))*(sig2/sqrt (2) -sig3/sqrt (2) - (2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -
    sqrt (3)*xi*fPhi)*sinbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta)) +
    2*(sqrt (0.6666666666666666)*sig1 - sig2/sqrt (6) -sig3/sqrt (6) - (2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -
    sqrt (3)*xi*fPhi)*cosbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta))*((2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -
    sqrt (3)*xi*fPhi)*cosbeta*(3*cos3beta - (3*cos3beta)/fPsi))/pow (1 + (1 - sin3beta)/fPsi + sin3beta,2) +
    (2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -sqrt (3)*xi*fPhi)*sinbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta)))/(2.*fG);
    TPZFMatrix<REAL> resf1(2,1,0.);
    resf1(0,0)=ddistf1dxi;
    resf1(1,0)=ddistf1dbeta;
    return resf1;
    
}

TPZFMatrix<REAL> TPZSandlerExtended::DDistFunc2(TPZTensor<REAL>::TPZDecomposed pt,REAL theta,REAL beta,REAL k,REAL kprev)
{
    REAL costheta, sintheta,c1,c2,c3,c4,c5,c6,c7,sig1,sig2,sig3,sinbeta,cosbeta,sin3beta,cos3beta;
    sinbeta=sin(beta);
    cosbeta=cos(beta);
    cos3beta=cos(3*beta);
    sin3beta=sin(3*beta);
    sig1 = pt.fEigenvalues[0];
    sig2 = pt.fEigenvalues[1];
    sig3 = pt.fEigenvalues[2];
    c7=1. + (1. - sin3beta)/fPsi +sin3beta;
    costheta=cos(theta);
    sintheta =sin(theta);
    c1=(sig1+sig2+sig3)/sqrt(3.);
    c2=fR*F(k,fPhi)/sqrt(3.);
    c3= sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c4=2*sqrt(2)*(F(k,fPhi)-fN)*cos(beta)/c7;
    c5=(1./sqrt(2.))*(sig2-sig3);
    c6=2.*sqrt(2.)*(F(k,fPhi)-fN)*sin(beta)/c7;
    REAL ddistf2dtheta=(2.*c2*(c1 - k/sqrt(3.) - c2*costheta)*sintheta)/(3.*fK) +
    (-2*c4*costheta*(c3 - c4*sintheta) - 2.*c6*costheta*(c5 - c6*sintheta))/(2.*fG);

    
    c1=pow(sig1/sqrt(3.) + sig2/sqrt(3.) + sig3/sqrt(3.) - (k + fR*(fA - fC*exp(fB*k) - k*fPhi)*costheta)/sqrt(3.),2.)/(3.*fK);
    c2=c3;
    c3=2.*sqrt(2.)*(F(k,fPhi)- fN)*sintheta;
    c4=c5;
    c5=c3;
    c6=1./fPsi;
    REAL ddistf2dbeta =(2*(c2 - (c3*cosbeta)/(1 + c6*(1 - sin3beta) +sin3beta))*
                   ((c3*cosbeta*(3*cos3beta - 3*c6*cos3beta))/pow (1 + c6*(1 - sin3beta) + sin3beta,2) + (c3*sinbeta)/(1 + c6*(1 - sin3beta) + sin3beta)) +2*((c5*(3*cos3beta - 3*c6*cos3beta)*sinbeta)/pow (1 + c6*(1 - sin3beta) + sin3beta,2) - (c5*cosbeta)/(1 + c6*(1 - sin3beta) +sin3beta))*(c4 - (c5*sinbeta)/(1 + c6*(1 - sin3beta) +sin3beta)))/(2.*fG);
    

    REAL resL = ResL(pt, theta, beta,k,kprev);
    
    TPZFMatrix<REAL> resf2(3,1,0.);
    resf2(0,0)=ddistf2dtheta;
    resf2(1,0)=ddistf2dbeta;
    resf2(2,0)=resL;
    
    return  resf2 ;
}

TPZFMatrix<REAL> TPZSandlerExtended::D2DistFunc1(TPZTensor<REAL>::TPZDecomposed pt,REAL xi,REAL beta)
{

    REAL d2distf1dxixi,d2distf1dxibeta;
    REAL d2distf1dbetaxi,d2distf1dbetabeta;
    REAL sinbeta,cos3beta,sintheta,costheta,cosbeta,cos2beta,cos7beta,cos6beta,sin5beta,sin7beta,cos4beta,sin2theta,sin3beta,sig1,sig2,sig3;
    cosbeta=cos(beta);
    sin3beta=sin(3*beta);
    sin5beta=sin(5*beta);
    sinbeta=sin(beta);
    cos2beta=cos(2*beta);
    cos3beta=cos(3*beta);
    cos7beta=cos(7*beta);
    cos4beta=cos(4*beta);
    cos6beta=cos(6*beta);
    cos4beta=cos(4*beta);
    sin7beta=sin(7*beta);
    sig1 = pt.fEigenvalues[0];
    sig2 = pt.fEigenvalues[1];
    sig3 = pt.fEigenvalues[2];
    d2distf1dxixi=2/(3.*fK) + ((16*pow(cosbeta,2)*pow(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi,2))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) +
                 (16*pow(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi,2)*pow(sinbeta,2))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) +
                 (12*sqrt(2)*cosbeta*exp(sqrt(3)*fB*xi)*pow(fB,2)*fC*(sqrt(0.6666666666666666)*sig1 - sig2/sqrt(6) - sig3/sqrt(6) -
                (2*sqrt(2)*cosbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                 (12*sqrt(2)*exp(sqrt(3)*fB*xi)*pow(fB,2)*fC*sinbeta*(sig2/sqrt(2) - sig3/sqrt(2) - (2*sqrt(2)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta))/(2.*fG);
    
    d2distf1dxibeta=((4*sqrt(2)*cosbeta*(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi)*(3*cos3beta - (3*cos3beta)/fPsi)*
      (sqrt(0.6666666666666666)*sig1 - sig2/sqrt(6) - sig3/sqrt(6) - (2*sqrt(2)*cosbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/
     pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) + (4*sqrt(2)*(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi)*sinbeta*
                                                    (sqrt(0.6666666666666666)*sig1 - sig2/sqrt(6) - sig3/sqrt(6) - (2*sqrt(2)*cosbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/
     (1 + (1 - sin3beta)/fPsi + sin3beta) - (4*sqrt(2)*(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi)*sinbeta*
                                             ((-2*sqrt(2)*cosbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                                              (2*sqrt(2)*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2)))/(1 + (1 - sin3beta)/fPsi + sin3beta) -
     (4*sqrt(2)*cosbeta*(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi)*(sig2/sqrt(2) - sig3/sqrt(2) -
                                                                                  (2*sqrt(2)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
     (4*sqrt(2)*(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi)*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta*
      (sig2/sqrt(2) - sig3/sqrt(2) - (2*sqrt(2)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) -
     (4*sqrt(2)*cosbeta*(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi)*((2*sqrt(2)*cosbeta*(3*cos3beta - (3*cos3beta)/fPsi)*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/
                                                                                  pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) + (2*sqrt(2)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta))/
    (2.*fG);
    
    d2distf1dbetaxi=(2*((2*sqrt(2)*cosbeta*(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi)*(3*cos3beta - (3*cos3beta)/fPsi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) +
                        (2*sqrt(2)*(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi)*sinbeta)/(1 + (1 - sin3beta)/fPsi + sin3beta))*
                     (sqrt(0.6666666666666666)*sig1 - sig2/sqrt(6) - sig3/sqrt(6) - (2*sqrt(2)*cosbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)) -
                     (4*sqrt(2)*(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi)*sinbeta*((-2*sqrt(2)*cosbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                                                                                                  (2*sqrt(2)*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2)))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                     2*((-2*sqrt(2)*cosbeta*(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                        (2*sqrt(2)*(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi)*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta)/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2))*
                     (sig2/sqrt(2) - sig3/sqrt(2) - (2*sqrt(2)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)) -
                     (4*sqrt(2)*cosbeta*(-(sqrt(3)*exp(sqrt(3)*fB*xi)*fB*fC) - sqrt(3)*fPhi)*((2*sqrt(2)*cosbeta*(3*cos3beta - (3*cos3beta)/fPsi)*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/
                                                                                                  pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) + (2*sqrt(2)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta))/
    (2.*fG);
    
    
    d2distf1dbetabeta=(2*(sqrt(0.6666666666666666)*sig1 - sig2/sqrt(6) - sig3/sqrt(6) - (2*sqrt(2)*cosbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta))*
                       ((-4*sqrt(2)*cosbeta*pow(3*cos3beta - (3*cos3beta)/fPsi,2)*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,3) +
                        (2*sqrt(2)*cosbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                        (2*sqrt(2)*cosbeta*(-9*sin3beta + (9*sin3beta)/fPsi)*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) -
                        (4*sqrt(2)*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2)) +
                       2*pow((-2*sqrt(2)*cosbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                               (2*sqrt(2)*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2),2) +
                       2*pow((2*sqrt(2)*cosbeta*(3*cos3beta - (3*cos3beta)/fPsi)*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) +
                               (2*sqrt(2)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta),2) +
                       2*(sig2/sqrt(2) - sig3/sqrt(2) - (2*sqrt(2)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta))*
                       ((4*sqrt(2)*cosbeta*(3*cos3beta - (3*cos3beta)/fPsi)*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) -
                        (4*sqrt(2)*pow(3*cos3beta - (3*cos3beta)/fPsi,2)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,3) +
                        (2*sqrt(2)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                        (2*sqrt(2)*(-9*sin3beta + (9*sin3beta)/fPsi)*sinbeta*(fA - exp(sqrt(3)*fB*xi)*fC - fN - sqrt(3)*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2)))/(2.*fG);
    
    TPZFMatrix<REAL> JAC(2,2,0);
    JAC(0,0)=d2distf1dxixi;JAC(0,1)=d2distf1dxibeta;
    JAC(1,0)=d2distf1dbetaxi;JAC(1,1)=d2distf1dbetabeta;
    return JAC;
    
}

TPZFMatrix<REAL> TPZSandlerExtended::D2DistFunc2(TPZTensor<REAL>::TPZDecomposed pt,REAL theta,REAL beta,REAL k)
{
    REAL d2distf2dthetatheta,d2distf2dthetabeta,d2distf2dthetak;
    REAL d2distf2dbetatheta,d2distf2dbetabeta,d2distf2dbetak;
    REAL dresktheta,dreskbeta,dreskk;
    REAL sinbeta,cos3beta,sintheta,costheta,cosbeta,cos2beta,cos7beta,cos6beta,sin5beta,sin7beta,cos4beta,sin2theta,sin3beta,cos2theta,c1,c2,c3,c4,c5,c6,c7,c8,c9,gamma,sig1,sig2,sig3;
    cos2theta=cos(2*theta);
    sintheta=sin(theta);
    costheta=cos(theta);
    cosbeta=cos(beta);
    sin3beta=sin(3*beta);
    sin5beta=sin(5*beta);
    sinbeta=sin(beta);
    cos2beta=cos(2*beta);
    cos3beta=cos(3*beta);
    cos7beta=cos(7*beta);
    cos4beta=cos(4*beta);
    cos6beta=cos(6*beta);
    cos4beta=cos(4*beta);
    sin7beta=sin(7*beta);
    sig1 = pt.fEigenvalues[0];
    sig2 = pt.fEigenvalues[1];
    sig3 = pt.fEigenvalues[2];
    
    
    
    c1=(sig1+sig2+sig3)/sqrt(3.);
    c2=(fR*(fA - fC*exp(fB*k) - k*fPhi))/sqrt(3);
    c3= sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c4=(2*sqrt(2)*(fA - fC*exp(fB*k) - fN - k*fPhi)*cosbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta);
    c5=(1./sqrt(2.))*(sig2-sig3);
    c6=(2*sqrt(2)*(fA - fC*exp(fB*k) - fN - k*fPhi)*sinbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta);
    d2distf2dthetatheta=(2*c2*fG*(3*c1 - sqrt(3)*k)*costheta + (-6*pow(c2,2)*fG + 9*(pow(c4,2) + pow(c6,2))*fK)*cos2theta + 9*(c3*c4 + c5*c6)*fK*sintheta)/(9.*fG*fK);
    
    
    
    c1=(2*fR*(fA - exp(fB*k)*fC - fPhi*k)*(-((k + costheta*fR*(fA - exp(fB*k)*fC - fPhi*k))/sqrt(3)) + sig1/sqrt(3) + sig2/sqrt(3) + sig3/sqrt(3))*sintheta)/(3.*sqrt(3)*fK);
    c2=4*sqrt(2)*costheta*(fA - exp(fB*k)*fC - fN - fPhi*k);
    c3= sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c4=2*sqrt(2)*(fA - exp(fB*k)*fC - fN - fPhi*k)*sintheta;
    c5=1/fPsi;
    c6=(1./sqrt(2.))*(sig2-sig3);
    d2distf2dthetabeta=((c2*(3*cos3beta - 3*c5*cos3beta)*cosbeta*(c3 - (c4*cosbeta)/(1 + c5*(1 - sin3beta) + sin3beta)))/pow(1 + c5*(1 - sin3beta) + sin3beta,2) +
                        (c2*(c3 - (c4*cosbeta)/(1 + c5*(1 - sin3beta) + sin3beta))*sinbeta)/(1 + c5*(1 - sin3beta) + sin3beta) -
                        (c2*sinbeta*(-((c4*cosbeta)/(1 + c5*(1 - sin3beta) + sin3beta)) + (c4*(3*cos3beta - 3*c5*cos3beta)*sinbeta)/pow(1 + c5*(1 - sin3beta) + sin3beta,2)))/
                        (1 + c5*(1 - sin3beta) + sin3beta) - (c2*cosbeta*(c6 - (c4*sinbeta)/(1 + c5*(1 - sin3beta) + sin3beta)))/(1 + c5*(1 - sin3beta) + sin3beta) +
                        (c2*(3*cos3beta - 3*c5*cos3beta)*sinbeta*(c6 - (c4*sinbeta)/(1 + c5*(1 - sin3beta) + sin3beta)))/pow(1 + c5*(1 - sin3beta) + sin3beta,2) -
                        (c2*cosbeta*((c4*(3*cos3beta - 3*c5*cos3beta)*cosbeta)/pow(1 + c5*(1 - sin3beta) + sin3beta,2) + (c4*sinbeta)/(1 + c5*(1 - sin3beta) + sin3beta)))/
                        (1 + c5*(1 - sin3beta) + sin3beta))/(2.*fG);
    
    
    c1=(sig1+sig2+sig3)/sqrt(3.);
    c2=fR*costheta;
    c3=2*fR*sintheta;
    c4=2*sqrt(2)*cosbeta*costheta;
    c5=sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c6=sqrt(2)*cosbeta*sintheta;
    c7=2*sqrt(2)*costheta*sinbeta;
    c8=(1./sqrt(2.))*(sig2-sig3);
    c9=sqrt(2)*sinbeta*costheta;
    gamma=0.5*(1 + (1 - sin3beta)/fPsi +  sin3beta);
    d2distf2dthetak=-(c3*(1 + c2*(-(exp(fB*k)*fB*fC) - fPhi))*(fA - exp(fB*k)*fC - fPhi*k))/(9.*fK) + (c3*(-(exp(fB*k)*fB*fC) - fPhi)*(c1 - (k + c2*(-(fC*exp(fB*k)) + fA - fPhi*k))/sqrt(3)))/(3.*sqrt(3)*fK) +
    ((c4*c6*(-(exp(fB*k)*fB*fC) - fPhi)*(fA - exp(fB*k)*fC - fN - fPhi*k))/pow(gamma,2) + (c7*c9*(-(exp(fB*k)*fB*fC) - fPhi)*(fA - exp(fB*k)*fC - fN - fPhi*k))/pow(gamma,2) -
     (c4*(-(exp(fB*k)*fB*fC) - fPhi)*(c5 - (c6*(fA - exp(fB*k)*fC - fN - fPhi*k))/gamma))/gamma - (c7*(-(exp(fB*k)*fB*fC) - fPhi)*(c8 - (c9*(fA - exp(fB*k)*fC - fN - fPhi*k))/gamma))/gamma)/(2.*fG);
    
    
    
    //////////////
    
    c1=(2*sqrt(2)*(3*cos3beta - (3*cos3beta)/fPsi)*(fA - exp(fB*k)*fC - fN - fPhi*k)*sinbeta)/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2);
    c2=(2*sqrt(2)*cosbeta*(fA - exp(fB*k)*fC - fN - fPhi*k))/(1 + (1 - sin3beta)/fPsi + sin3beta);
    c3=(1./sqrt(2.))*(sig2-sig3);
    c4=(2*sqrt(2)*(fA - exp(fB*k)*fC - fN - fPhi*k)*sinbeta)/(1 + (1 - sin3beta)/fPsi + sin3beta);
    c5=sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c6=(2*sqrt(2)*(3*cos3beta - (3*cos3beta)/fPsi)*(fA - exp(fB*k)*fC - fN - fPhi*k)*cosbeta)/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2);
    c7=(2*sqrt(2)*sinbeta*(fA - exp(fB*k)*fC - fN - fPhi*k))/(1 + (1 - sin3beta)/fPsi + sin3beta);
    d2distf2dbetatheta=(2*(c6*costheta + c7*costheta)*(c5 - c2*sintheta) - 2*c4*costheta*(c1*sintheta - c2*sintheta) + 2*(c1*costheta - c2*costheta)*(c3 - c4*sintheta) - 2*c2*costheta*(c6*sintheta + c7*sintheta))/(2.*fG);
    
    c1=pow(-((k + costheta*fR*(fA - exp(fB*k)*fC - fPhi*k))/sqrt(3)) + sig1/sqrt(3) + sig2/sqrt(3) + sig3/sqrt(3),2)/(3.*fK);
    c2=sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c3=2*sqrt(2)*sintheta*(fA - exp(fB*k)*fC - fN - k*fPhi);
    c4=(1./sqrt(2.))*(sig2-sig3);
    c5=2*sqrt(2)*(fA -exp(fB*k)*fC - fN - fPhi*k)*sintheta;
    c6=1/fPsi;
    d2distf2dbetabeta=(2*pow(c5,2)*pow(cosbeta,2)*pow(1 + c6 + 2*(-1 + c6)*sin3beta - 6*(-1 + c6)*sinbeta,2) + 2*pow(c3*(-1 + c6)*(2*cos2beta + cos4beta) - c3*(1 + c6)*sinbeta,2) +2*c3*cosbeta*(-(c2*(1 + c6)) + c3*cosbeta + c2*(-1 + c6)*sin3beta)*(15 - 34*c6 + 15*pow(c6,2) - 6*pow(-1 + c6,2)*
    cos2beta + 6*pow(-1 + c6,2)*cos4beta + 2*pow(-1 + c6,2)*cos6beta -13*(-1 + pow(c6,2))*sin3beta + 12*(-1 + pow(c6,2))*sinbeta) + c5*(-(c4*(1 + c6)) + c4*(-1 + c6)*sin3beta + c5*sinbeta)*((1 - pow(c6,2))*cos2beta + 13*(-1 + pow(c6,2))*cos4beta + 2*(pow(-1 + c6,2)*(-4*sin5beta + sin7beta) + 4*(3 + c6*(-7 + 3*c6))*sinbeta)))/(2.*fG*pow(1 + c6 - (-1 + c6)*sin3beta,4));
    
    
    c1=(2*sqrt(2)*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta*sintheta)/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2);
    c2=(2*sqrt(2)*cosbeta*sintheta)/(1 + (1 - sin3beta)/fPsi + sin3beta);
    c3=(2*sqrt(2)*sinbeta*sintheta)/(1 + (1 - sin3beta)/fPsi + sin3beta);
    c4=(1./sqrt(2.))*(sig2-sig3);
    c5=c2;
    c6=sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c7=(2*sqrt(2)*(3*cos3beta - (3*cos3beta)/fPsi)*cosbeta*sintheta)/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2);
    c8=c3;
    d2distf2dbetak=(-2*c3*(-(exp(fB*k)*fB*fC) - fPhi)*(c1*(fA - exp(fB*k)*fC - fN - fPhi*k) - c2*(fA - exp(fB*k)*fC - fN - fPhi*k)) +
     2*(c1*(-(exp(fB*k)*fB*fC) - fPhi) - c2*(-(exp(fB*k)*fB*fC) - fPhi))*(c4 - c3*(fA - exp(fB*k)*fC - fN - fPhi*k)) +
     2*(c7*(-(exp(fB*k)*fB*fC) - fPhi) + c8*(-(exp(fB*k)*fB*fC) - fPhi))*(c6 - c5*(fA - exp(fB*k)*fC - fN - fPhi*k)) -
     2*c5*(-(exp(fB*k)*fB*fC) - fPhi)*(c7*(fA - exp(fB*k)*fC - fN - fPhi*k) + c8*(fA - exp(fB*k)*fC - fN - fPhi*k)))/(2.*fG);
    
    
    
    dreskk=1 + costheta*(-(exp(fB*k)*fB*fC) - fPhi)*fR +
    3*exp(fD*(-((fA - exp(fB*k)*fC)*fR) + k))*fD*fK*
    (1 + exp(fB*k)*fB*fC*fR)*fW;
    
    
    dresktheta=-fR*(fA - fC * exp(fB*k) - k*fPhi)* sintheta;
    
    TPZFMatrix<REAL> JAC(3,3,0);
    JAC(0,0)= d2distf2dthetatheta;JAC(0,1)=d2distf2dthetabeta;JAC(0,2)=d2distf2dthetak;
    JAC(1,0)=d2distf2dbetatheta;JAC(1,1)=d2distf2dbetabeta;JAC(1,2)=d2distf2dbetak;
    JAC(2,0)=dresktheta;JAC(2,1)=0;JAC(2,2)=dreskk;
    return JAC;
}



void TPZSandlerExtended::YieldFunction(TPZTensor<REAL>::TPZDecomposed &sigma, TPZVec<REAL> &yield,REAL &kprev)
{

    REAL I1,J2,J3,beta,temp,gamma;
    sigma.ComputeI1(I1);
    sigma.ComputeJ2(J2);
    sigma.ComputeJ3(J3);
    if (J2<1.e-6) {
        J2=1.e-6;
    }
    temp=J3/(sqrt(J2*J2*J2));
    beta=(1/3)*acos((3*sqrt(3)/2)*temp);
    gamma = 0.5*(1 + (1 - sin(3*beta))/fPsi + sin(3*beta));
    yield[0]=sqrt(J2)-(pow(((I1-kprev)/(fR*F(kprev,fPhi)))+(gamma*sqrt(J2)/F(kprev,fPhi)),2)-1);
    yield[1]=sqrt(J2)-F(I1,fPhi);
    
    
}

void TPZSandlerExtended::ProjectF1(TPZTensor<REAL>::TPZDecomposed &sigmatrial, TPZTensor<REAL>::TPZDecomposed &sigproj)
{
    REAL xi=-M_PI;
    REAL deltaxi=M_PI/100;
    REAL beta = 0;
    REAL deltabeta=2*M_PI/200;
    //TPZFMatrix<REAL> store(20,20,0.);
    TPZManVector< std::pair<REAL,REAL>,200 > pairreal(200);
    TPZManVector<REAL,200> dist(200,0.);
    int count=0;
    for(int i=0;i<200;i++)
    {
        dist[count]=sqrt(DistF1(sigmatrial, xi,beta));
        pairreal[count]=std::make_pair(xi,beta);
        xi+=deltaxi;
        beta+=deltabeta;
        count++;
    }
    REAL temp;
    pair<REAL,REAL> temp2;
    for (int i=0; i<199;i++)
    {
        
        for (int j=0; j<199; j++)
        {
            
            if(dist[j]>dist[j+1])
            {
                temp=dist[j];
                temp2=pairreal[j];
                
                pairreal[j]=pairreal[j+1];
                pairreal[j+1]=temp2;
                
                dist[j]=dist[j+1];
                dist[j+1]=temp;
                
            }
            
        }
        
    }
    REAL xiguess,betaguess,resnorm;
    xiguess=pairreal[0].first;
    betaguess=pairreal[0].second;
    resnorm=1;
    int counter=1;
    TPZFMatrix<REAL> xn1(2,1,0.),xn(2,1,0.),jac,invjac,sol(2,1,0.),fxn(2,1,0.),diff(2,1,0.);
    xn(0,0)=xiguess;
    xn(1,0)=betaguess;
    while (resnorm > 10e-12 && counter < 30)
    {

        jac=D2DistFunc1(sigmatrial, xn[0],xn[1]);
        fxn=DDistFunc1(sigmatrial, xn[0],xn[1]);
        jac.Inverse(invjac);
        invjac.Multiply(fxn,sol);
        xn1=xn-sol;
        diff=xn1-xn;
        resnorm=Norm(diff);
        xn=xn1;
        counter++;
        
    }

    REAL sin3beta,gamma,I1,sqrtj2;
    sin3beta = sin(3*xn1[1]);
    gamma=0.5*(1 + (1 - sin3beta)/fPsi +  sin3beta);
    I1 = xn1[0]*sqrt(3);
    sqrtj2 = (F(I1,fPhi) -fN)/gamma;
    TPZVec<REAL> i1sqrtj2(2,0.);
    i1sqrtj2[0]=I1;
    i1sqrtj2[1]=sqrtj2;
    TPZTensor<REAL> sigtrial,sigcor;
    sigtrial.XX()=sigmatrial.fEigenvalues[0];
    sigtrial.YY()=sigmatrial.fEigenvalues[1];
    sigtrial.ZZ()=sigmatrial.fEigenvalues[2];
    sigtrial.Adjust(i1sqrtj2,sigcor);
    TPZTensor<REAL>::TPZDecomposed tempp(sigcor);
    sigproj=tempp;
    
}

void TPZSandlerExtended::ProjectF2(TPZTensor<REAL>::TPZDecomposed &sigmatrial, TPZTensor<REAL>::TPZDecomposed &sigproj,REAL &kprev)
{
    REAL kn,theta,beta,deltabeta,deltatheta;
    Firstk(kn);
    
    theta=0.;
    beta=0;
    int nsearch=20;
    deltabeta=(2*M_PI)/20;
    deltatheta=(M_PI/2)/20;
    TPZManVector< std::pair<REAL,REAL>,20 > pairreal(20);
    TPZManVector<REAL,20> dist(20,0.);
    int count=0;
    for(int i=0;i<20;i++)
    {
        dist[count]=DistF2(sigmatrial, theta,beta,kn);
        pairreal[count]=std::make_pair(theta,beta);
        theta+=deltatheta;
        beta+=deltabeta;
        count++;
    }
    REAL temp;
    pair<REAL,REAL> temp2;
    for (int i=0; i<19;i++)
    {
        
        for (int j=0; j<19; j++)
        {
            
            if(dist[j]>dist[j+1])
            {
                temp=dist[j];
                temp2=pairreal[j];
                
                pairreal[j]=pairreal[j+1];
                pairreal[j+1]=temp2;
                
                dist[j]=dist[j+1];
                dist[j+1]=temp;
                
            }
            
        }
        
    }
    REAL thetaguess,betaguess,resnorm;
    thetaguess=pairreal[0].first;
    betaguess=pairreal[0].second;
    resnorm=1;
    int counter=1;
    TPZFMatrix<REAL> xn1(3,1,0.),xn(3,1,0.),jac,invjac,sol(3,1,0.),fxn(3,1,0.),diff(3,1,0.);
    xn(0,0)=thetaguess;
    xn(1,0)=betaguess;
    xn(2,0)=kn;
    while (resnorm > 10e-12 && counter < 30)
    {
        jac=D2DistFunc2(sigmatrial, xn[0],xn[1],xn[2]);
        fxn=DDistFunc2(sigmatrial, xn[0],xn[1],xn[2],kn);
        jac.Inverse(invjac);
        invjac.Multiply(fxn,sol);
        xn1=xn-sol;
        diff=xn1-xn;
        resnorm=Norm(diff);
        xn=xn1;
        counter++;
        cout<< "\n resnorm = "<<resnorm <<endl;
        cout<< "\n xn1 = "<<xn1 <<endl;
        cout<< "\n deltak = "<<xn1[2]-kn <<endl;
    }
    
    REAL sin3beta,gamma,I1,sqrtj2,thetasol,betasol,ksol;
    
    thetasol=xn1[0];
    betasol=xn1[1];
    ksol=xn1[2];
    sin3beta = sin(3*betasol);
    gamma=0.5*(1 + (1 - sin3beta)/fPsi +  sin3beta);
    I1 = ksol*fR*F(ksol,fPhi)*cos(thetasol);
    sqrtj2 = (F(ksol,fPhi) -fN)*sin(thetasol)/gamma;
    TPZVec<REAL> i1sqrtj2(2,0.);
    i1sqrtj2[0]=I1;
    i1sqrtj2[1]=sqrtj2;
    TPZTensor<REAL> sigtrial,sigcor;
    sigtrial.XX()=sigmatrial.fEigenvalues[0];
    sigtrial.YY()=sigmatrial.fEigenvalues[1];
    sigtrial.ZZ()=sigmatrial.fEigenvalues[2];
    sigtrial.Adjust(i1sqrtj2,sigcor);
    TPZTensor<REAL>::TPZDecomposed tempp(sigcor);
    sigproj=tempp;
}

void TPZSandlerExtended::ProjectRing(TPZTensor<REAL>::TPZDecomposed &sigmatrial, TPZTensor<REAL>::TPZDecomposed &sigproj,REAL &kprev)
{
    REAL kn,theta,beta,deltabeta,deltatheta;
    Firstk(kn);
    
    theta=0.;
    beta=0;
    int nsearch=20;
    deltabeta=(2*M_PI)/20;
    TPZManVector< std::pair<REAL,REAL>,20 > pairreal(20);
    TPZManVector<REAL,20> dist(20,0.);
    int count=0;
    for(int i=0;i<20;i++)
    {
        dist[count]=DistF2(sigmatrial, theta,beta,kn);
        pairreal[count]=std::make_pair(theta,beta);
        theta=M_PI/2;
        beta+=deltabeta;
        count++;
    }
    REAL temp;
    pair<REAL,REAL> temp2;
    for (int i=0; i<19;i++)
    {
        
        for (int j=0; j<19; j++)
        {
            
            if(dist[j]>dist[j+1])
            {
                temp=dist[j];
                temp2=pairreal[j];
                
                pairreal[j]=pairreal[j+1];
                pairreal[j+1]=temp2;
                
                dist[j]=dist[j+1];
                dist[j+1]=temp;
                
            }
            
        }
        
    }
    REAL thetaguess,betaguess,resnorm;
    thetaguess=pairreal[0].first;
    betaguess=pairreal[0].second;
    resnorm=1;
    int counter=1;
    TPZFMatrix<REAL> xn1(3,1,0.),xn(3,1,0.),jac,invjac,sol(3,1,0.),fxn(3,1,0.),diff(3,1,0.);
    xn(0,0)=M_PI/2;
    xn(1,0)=betaguess;
    xn(2,0)=kn;
    while (resnorm > 10e-12 && counter < 30)
    {
        jac=D2DistFunc2(sigmatrial,xn[0],xn[1],xn[2]);
        fxn=DDistFunc2(sigmatrial, xn[0],xn[1],xn[2],kn);
        jac.Inverse(invjac);
        invjac.Multiply(fxn,sol);
        
        xn1(0,0)=xn(0,0);
        xn1(1,0)=xn(1,0)-sol(1,0);
        xn1(2,0)=xn(2,0)-sol(2,0);
        
        diff=xn1-xn;
        resnorm=Norm(diff);
        xn=xn1;
        counter++;
        cout<< "\n resnorm = "<<resnorm <<endl;
        cout<< "\n xn1 = "<<xn1 <<endl;
        cout<< "\n deltak = "<<xn1[2]-kn <<endl;
    }
    
    REAL sin3beta,gamma,I1,sqrtj2,thetasol,betasol,ksol;
    
    thetasol=xn1[0];
    betasol=xn1[1];
    ksol=xn1[2];
    sin3beta = sin(3*betasol);
    gamma=0.5*(1 + (1 - sin3beta)/fPsi +  sin3beta);
    I1 = ksol*fR*F(ksol,fPhi)*cos(thetasol);
    sqrtj2 = (F(ksol,fPhi) -fN)*sin(thetasol)/gamma;
    TPZVec<REAL> i1sqrtj2(2,0.);
    i1sqrtj2[0]=I1;
    i1sqrtj2[1]=sqrtj2;
    TPZTensor<REAL> sigtrial,sigcor;
    sigtrial.XX()=sigmatrial.fEigenvalues[0];
    sigtrial.YY()=sigmatrial.fEigenvalues[1];
    sigtrial.ZZ()=sigmatrial.fEigenvalues[2];
    sigtrial.Adjust(i1sqrtj2,sigcor);
    TPZTensor<REAL>::TPZDecomposed tempp(sigcor);
    sigproj=tempp;
}

/**
 * Imposes the specified strain tensor and returns the correspondent stress state.
 *
 * @param[in] epsTotal Imposed total strain tensor
 * @param[out] sigma Resultant stress
 */
void TPZSandlerExtended::ApplyStrainComputeSigma(TPZPlasticState<REAL> &plasticstate,const TPZElasticResponse &ER, TPZTensor<REAL> &sigma)
{
    
}

