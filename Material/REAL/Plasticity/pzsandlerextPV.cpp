//
//  pzsandlerextPV.cpp
//  PZ
//
//  Created by Diogo Cecilio on 9/3/13.
//
//

#include "pzsandlerextPV.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.poroelastoplastic"));
#endif

TPZSandlerExtended::TPZSandlerExtended()
{
    
}

TPZSandlerExtended::TPZSandlerExtended(const TPZSandlerExtended & copy)
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
    fE=(9.*fK*fG)/(3.*fK+fG);
    fnu=((3.*fK)-(2.*fG))/(2*(3.*fK+fG));
    
}

TPZSandlerExtended::TPZSandlerExtended(STATE A, STATE B,STATE C, STATE D,STATE K,STATE G,STATE W,STATE R,STATE Phi,STATE N,STATE Psi):
fA(A),fB(B),fC(C),fD(D),fK(K),fG(G),fW(W),fR(R),fPhi(Phi),fN(N),fPsi(Psi)
{
    fE=(9.*fK*fG)/(3.*fK+fG);
    fnu=((3.*fK)-(2.*fG))/(2*(3.*fK+fG));
}

TPZSandlerExtended::~TPZSandlerExtended()
{
    
}

template <class T>
T TPZSandlerExtended::F(T x,STATE phi) const
{
    return (fA - fC*exp(x *fB) -fPhi*x);
}

template<class T>
T TPZSandlerExtended::X(T k) const
{
    return (k - fR * F(k, fPhi));
}

template<class T>
T TPZSandlerExtended::EpsEqX(T X) const
{
    return (fW*( exp(fD*X) - 1 ));
//   return fW* exp(fD*X);
}

template<class T>
T TPZSandlerExtended::EpsEqk(T k) const
{
    return EpsEqX(X(k));
}

void TPZSandlerExtended::Firstk(STATE &epsp,STATE &k) const
{
    STATE f,df,kn1,kn,resnorm,diff;
    int counter =1;
    resnorm=1;
    kn=epsp;//chute inicial
    while (resnorm>1.e-12 && counter<30) {
        
        f=EpsEqk(kn)-epsp;
        df =fD*exp(fD*(kn - (fA - fC*exp(fB*kn))*fR))*(1 + fB*fC*exp(fB*kn)*fR)*fW;
        //df=fD*exp(fD*(kn - fR*(fA - fC*exp(fB*kn) - kn*fPhi)))*fW*(1 - fR*(-(fB*fC*exp(fB*kn)) - fPhi));
        kn1=kn-f/df;
        diff=kn1-kn;
        resnorm=sqrt(diff*diff);
        kn=kn1;
        counter++;
        
    }
    k=kn1;
}

template<class T>
T TPZSandlerExtended::ResLF2(const TPZVec<T> &pt, T theta,T beta,T k,STATE kprev ) const
{
    
    T I1tr=(pt[0])+(pt[1])+(pt[2]);
    T I1 =fR*F(k,fPhi)* cos(theta) +k;
    T delepsp = EpsEqk(k) - EpsEqk(kprev);
    return  (3.*fK*delepsp - (I1tr - I1));
}

/// Compute the residual of the equation which defines the update of the damage variable
STATE TPZSandlerExtended::ResLF1(const TPZVec<STATE> &sigtrial, TPZVec<STATE> &sigproj,STATE k,STATE kprev ) const
{
    STATE I1tr=(sigtrial[0]+sigtrial[1]+sigtrial[2]);
    STATE I1 =(sigproj[0]+sigproj[1]+sigproj[2]);
    STATE delepsp = EpsEqk(k) - EpsEqk(kprev);
    return  (3.*fK*delepsp - (I1tr - I1));
    
}

/// Compute the derivative of the equation which determines the evolution of k
// the derivative are given in terms of k
STATE TPZSandlerExtended::DResLF1(const TPZVec<STATE> &sigtrial, const TPZVec<STATE> &sigproj, STATE k, STATE kprev) const
{
    STATE expfBk=exp(fB*k);
    STATE dreskk= 3.*exp(fD*(-((fA - expfBk*fC)*fR) + k))*fD*fK*
    (1. + expfBk*fB*fC*fR)*fW;
    return dreskk;
}


/// Compute the derivative of the equation which determines the evolution of k
void TPZSandlerExtended::DResLF2(const TPZVec<STATE> &pt, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &dresl) const
{
    STATE expfBk=exp(fB*k);
    STATE sintheta = sin(theta);
    STATE costheta = cos(theta);

    STATE dreskk=1. + costheta*(-(expfBk*fB*fC) - fPhi)*fR +
    3.*exp(fD*(-((fA - expfBk*fC)*fR) + k))*fD*fK*
    (1. + expfBk*fB*fC*fR)*fW;
    
    
    STATE dresktheta=-fR*(fA - fC * expfBk - k*fPhi)* sintheta;
    dresl[0]=dresktheta;
    dresl[1]=0;
    dresl[2]=dreskk;

    
}

void TPZSandlerExtended::FromHWCylToPrincipal(const TPZVec<STATE> &HWCylCoords, TPZVec<STATE> &HWCart)
{

    HWCart[0]=(1./sqrt(3.))*HWCylCoords[0]+sqrt(2./3.)*HWCylCoords[1]*cos(HWCylCoords[2]);
    HWCart[1]=(1./sqrt(3.))*HWCylCoords[0]+sqrt(2./3.)*HWCylCoords[1]*cos(HWCylCoords[2]-(2.*M_PI/3.));
    HWCart[2]=(1./sqrt(3.))*HWCylCoords[0]+sqrt(2./3.)*HWCylCoords[1]*cos(HWCylCoords[2]+(2.*M_PI/3.));
}

void  TPZSandlerExtended::FromHWCylToHWCart(const TPZVec<STATE> &HWCylCoords, TPZVec<STATE> &cart)
{
    cart[0]=HWCylCoords[0];
    cart[1]=HWCylCoords[1]*cos(HWCylCoords[2]);
    cart[2]=HWCylCoords[1]*sin(HWCylCoords[2]);
}

void TPZSandlerExtended::FromPrincipalToHWCart(const TPZVec<STATE> &PrincipalCoords, TPZVec<STATE> &HWCart)
{
    TPZFNMatrix<9,STATE> Rot(3,3,0.),temp(3,1,0.),cart(3,1,0.);
    temp(0,0)=PrincipalCoords[0];
    temp(1,0)=PrincipalCoords[1];
    temp(2,0)=PrincipalCoords[2];
    GetRotMatrix(Rot);
    Rot.Multiply(temp,cart);
    HWCart[0]=cart(0,0);
    HWCart[1]=cart(1,0);
    HWCart[2]=cart(2,0);
}

void TPZSandlerExtended::FromPrincipalToHWCyl(const TPZVec<STATE> &PrincipalCoords, TPZVec<STATE> &HWCyl)
{
    TPZFNMatrix<9,STATE> Rot(3,3,0.),temp(3,1,0.),cart(3,1,0.);
    temp(0,0)=PrincipalCoords[0];
    temp(1,0)=PrincipalCoords[1];
    temp(2,0)=PrincipalCoords[2];
    GetRotMatrix(Rot);
    Rot.Multiply(temp,cart);
    HWCyl[0]=cart(0,0);
    HWCyl[1]=sqrt(cart(1,0)*cart(1,0)+cart(2,0)*cart(2,0));
    HWCyl[2]=atan2(cart(2,0),cart(1,0));
//    HWCyl[2]=atan(cart(2,0)/cart(1,0));
}


void TPZSandlerExtended::F1Cyl(STATE xi, STATE beta, TPZVec<STATE> &f1cyl) const
{
    STATE sqrt2=sqrt(2);
    STATE sqrt3=sqrt(3);
    STATE gamma = 0.5*(1 + (1 - sin(3*beta))/fPsi + sin(3*beta));
    STATE I1 = xi*sqrt3;
    STATE sqrtj2 = ( F(I1, fPhi) - fN )/gamma;
    STATE rho = sqrt2*sqrtj2;
    f1cyl[0]=xi,
    f1cyl[1]=rho;
    f1cyl[2]=beta;
    
}

void TPZSandlerExtended::SurfaceParamF1(TPZVec<STATE> &sigproj, STATE &xi, STATE &beta) const
{
    TPZManVector<STATE> sigHWCyl(3);
    FromPrincipalToHWCyl(sigproj, sigHWCyl);
    xi = sigHWCyl[0];
    beta = sigHWCyl[2];
#ifdef DEBUG
    STATE dist = DistF1(sigproj, xi, beta);
    if (fabs(dist) > 1.e-8) {
        DebugStop();
    }
#endif
}



void TPZSandlerExtended::F2Cyl(STATE theta, STATE beta,STATE k, TPZVec<STATE> &f2cyl) const
{
    STATE sqrt2=sqrt(2);
    STATE sqrt3=sqrt(3);
    STATE gamma = 0.5*(1 + (1 - sin(3*beta))/fPsi + sin(3*beta));
    STATE var =fR*F(k, fPhi)*cos(theta);
    STATE I1 = k + var;
    STATE sqrtj2 = (F(k, fPhi)- fN)*sin(theta)/gamma;
    STATE rho = sqrt2*sqrtj2;
    STATE xi=I1/sqrt3;
    f2cyl[0]=xi,
    f2cyl[1]=rho;
    f2cyl[2]=beta;
    
}

void TPZSandlerExtended::SurfaceParamF2(TPZVec<STATE> &sigproj, STATE k, STATE &theta, STATE &beta) const
{
    TPZManVector<STATE> sigHWCyl(3);
    FromPrincipalToHWCyl(sigproj, sigHWCyl);
    STATE xi,rho;
    xi=sigHWCyl[0];
    rho=sigHWCyl[1];
    STATE I1 = sigHWCyl[0]*sqrt(3.);
    beta = sigHWCyl[2];
    STATE gamma = 0.5*(1 + (1 - sin(3*beta))/fPsi + sin(3*beta));
    STATE costheta = (I1-k)/(fR*F(k,fPhi));
    STATE sqrtj2 = sigHWCyl[1]/sqrt(2.);
    STATE sintheta = sqrtj2*gamma/(F(k,fPhi)-fN);
    theta = atan2(sintheta, costheta);
    theta = acos(costheta);
    //STATE theta2 = atan((rho*sin(beta))/xi);
#ifdef DEBUG
    STATE err = 1.-sintheta*sintheta-costheta*costheta;
    STATE dist = DistF2(sigproj, theta, beta, k);
    if (fabs(dist) > 1.e-8 || err > 1.e-8) {
        DebugStop();
    }
#endif
    
}

void TPZSandlerExtended::GetRotMatrix(TPZFMatrix<STATE> &Rot)
{
    Rot.Resize(3,3);
    Rot(0,0)=1./sqrt(3.);    Rot(0,1)=1./sqrt(3.);    Rot(0,2)=1./sqrt(3.);
    Rot(1,0)=sqrt(2./3.);    Rot(1,1)=-1./sqrt(6.);   Rot(1,2)=-1./sqrt(6.);
    Rot(2,0)=0;              Rot(2,1)=1./sqrt(2.);    Rot(2,2)=-1./sqrt(2.);
}

STATE TPZSandlerExtended::DistF1(const TPZVec<STATE> &pt,STATE xi,STATE beta) const
{
    TPZFNMatrix<9,STATE> Rot(3,3);
    GetRotMatrix(Rot);
    TPZManVector<STATE,3> cyl(3);
    F1Cyl(xi,beta,cyl);
    TPZManVector<STATE,3> cart(3);
    FromHWCylToHWCart(cyl,cart);
    TPZManVector<STATE,3> carttrial(3);
    FromPrincipalToHWCart(pt,carttrial);
    return ((1./(3.*fK))*(carttrial[0]-cart[0])*(carttrial[0]-cart[0]))
    +(1./(2.*fG))*((carttrial[1]-cart[1])*(carttrial[1]-cart[1])+(carttrial[2]-cart[2])*(carttrial[2]-cart[2]));
    
}

STATE TPZSandlerExtended::DistF2(const TPZVec<STATE> &pt,STATE theta,STATE beta,STATE k) const
{
    TPZFNMatrix<9,STATE> Rot(3,3);
    GetRotMatrix(Rot);
    TPZManVector<STATE,3> cyl(3);
    F2Cyl(theta,beta,k,cyl);
    TPZManVector<STATE,3> cart(3);
    FromHWCylToHWCart(cyl,cart);
    TPZManVector<STATE,3> carttrial(3);
    FromPrincipalToHWCart(pt,carttrial);
    return ((1./(3.*fK))*(carttrial[0]-cart[0])*(carttrial[0]-cart[0]))
    +(1./(2.*fG))*((carttrial[1]-cart[1])*(carttrial[1]-cart[1])+(carttrial[2]-cart[2])*(carttrial[2]-cart[2]));
}


void TPZSandlerExtended::DDistFunc1(const TPZVec<STATE> &pt,STATE xi,STATE beta, TPZFMatrix<STATE> &ddistf1) const
{
    STATE sig1,sig2,sig3,sinbeta,cosbeta,sin3beta,cos3beta;
    sinbeta=sin(beta);
    cosbeta=cos(beta);
    cos3beta=cos(3*beta);
    sin3beta=sin(3*beta);
    sig1 = pt[0];
    sig2 = pt[1];
    sig3 = pt[2];
    STATE ddistf1dxi = (-2*(sig1/sqrt (3) + sig2/sqrt (3) + sig3/sqrt (3) - xi))/(3.*fK) +((-4*sqrt (2)*(-(sqrt (3)*fB*fC*exp (sqrt (3)*fB*xi)) -
                                                                                                        sqrt (3)*fPhi)*cosbeta*(sqrt (0.6666666666666666)*sig1 - sig2/sqrt (6) -sig3/sqrt (6) -(2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -
                                                                                                                                                                                                            sqrt (3)*xi*fPhi)*cosbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta) -
                                                                                          (4*sqrt (2)*(-(sqrt (3)*fB*fC*exp (sqrt (3)*fB*xi)) -sqrt (3)*fPhi)*sinbeta*(sig2/sqrt (2) - sig3/sqrt (2) -(2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -sqrt (3)*xi*fPhi)*sinbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta))/(2.*fG);
    
    STATE ddistf1dbeta=(2*((2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -sqrt (3)*xi*fPhi)*(3*cos3beta - (3*cos3beta)/fPsi)*
                           sinbeta)/pow (1 + (1 - sin3beta)/fPsi + sin3beta,2) - (2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -
                                                                                              sqrt (3)*xi*fPhi)*cosbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta))*(sig2/sqrt (2) -sig3/sqrt (2) - (2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -
                                                                                                                                                                                                           sqrt (3)*xi*fPhi)*sinbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta)) +
                       2*(sqrt (0.6666666666666666)*sig1 - sig2/sqrt (6) -sig3/sqrt (6) - (2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -
                                                                                                       sqrt (3)*xi*fPhi)*cosbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta))*((2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -
                                                                                                                                                                                     sqrt (3)*xi*fPhi)*cosbeta*(3*cos3beta - (3*cos3beta)/fPsi))/pow (1 + (1 - sin3beta)/fPsi + sin3beta,2) +
                                                                                                                                                                        (2*sqrt (2)*(fA - fC*exp (sqrt (3)*fB*xi) - fN -sqrt (3)*xi*fPhi)*sinbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta)))/(2.*fG);
    ddistf1(0,0)=ddistf1dxi;
    ddistf1(1,0)=ddistf1dbeta;
    
}

void TPZSandlerExtended::DDistFunc2new(const TPZVec<STATE> &pt,STATE theta,STATE beta,STATE k,STATE kprev, TPZManVector<STATE> &ddistf2) const
{
    STATE d2distf2dthetatheta,d2distf2dthetabeta,d2distf2dthetak;
    STATE d2distf2dbetatheta,d2distf2dbetabeta,d2distf2dbetak;
    STATE dresktheta,dreskk;
    STATE fR2=fR*fR,Gamma=1 + (1 - sin(3.*beta))/fPsi + sin(3.*beta),stresbeta=sin(3*beta),LPrev=kprev,FFLPrev=F(LPrev,fPhi),
    Gamma2=Gamma*Gamma,FFL=F(k,fPhi),FFL2=FFL*FFL,fN2=fN*fN,
    sb=sin(beta),cb=cos(beta),sb2=sb*sb,cb2=cb*cb,ct=cos(theta),st=sin(theta),
    ct2=ct*ct,st2=st*st,sig1=pt[0],sig2=pt[1],sig3=pt[2],sq2=sqrt(2.),sq3=sqrt(3.),
    L=k,FN=fN,xi=1./sq3*(sig1+sig2+sig3),sig2tiu=(sqrt(2./3.)*sig1)-(sig2/sqrt(6.))-(sig3/sqrt(6.)),
    DGamma=3*cos(3*beta) - (3*cos(3*beta))/fPsi,DGamma2=DGamma*DGamma,Gamma3=Gamma*Gamma*Gamma,Gamma4=Gamma2*Gamma2,expBC=exp(fB*L)*fB*fC + fPhi;

    STATE ddistf2dtheta=(-2*ct*(FFL - fN)*(Gamma*(sq2*cb*sig2tiu + sb*(sig2 - sig3)) + 4*(cb2*(-FFL + FN) + (-FFL + fN)*sb2)*st))/(fG*Gamma2) - (2*FFL*fR*st*(ct*FFL*fR + L - sq3*xi))/(9.*fK);
    STATE ddistf2dbeta =  (2*(FFL - fN)*st*(8*cb2*(-FFL + FN)*st*(3*cos(3*beta) - (3*cos(3*beta))/fPsi) + cb*Gamma*(Gamma*(-sig2 + sig3) + 4*(-fN + FN)*sb*st + 2*sq2*sig2tiu*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)) +sb*(sq2*Gamma2*sig2tiu + 2*(Gamma*(sig2 - sig3) + 4*(-FFL + fN)*sb*st)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi))))/(fG*Gamma3);
    STATE resL =ct*FFL*fR + L - sig1 - sig2 - sig3 + 3*(exp(fD*(-(FFL*fR) + L)) - exp(fD*(-(FFLPrev*fR) + LPrev)))*fK*fW;
    
    ddistf2[0]=ddistf2dtheta;
    ddistf2[1]=ddistf2dbeta;
    ddistf2[2]=resL;
    
}

template<class T>
void TPZSandlerExtended::DDistFunc2(const TPZVec<T> &pt,T theta,T beta,T k,STATE kprev, TPZVec<T> &ddistf2) const
{
    STATE sqrt2=sqrt(2);
    T expfBk=exp(fB*k);
    T costheta, sintheta,c1,c2,c3,c4,c5,c6,c7,sig1,sig2,sig3,sinbeta,cosbeta,sin3beta,cos3beta;
    sinbeta=sin(beta);
    cosbeta=cos(beta);
    cos3beta=cos(3.*beta);
    sin3beta=sin(3.*beta);
    sig1 = pt[0];
    sig2 = pt[1];
    sig3 = pt[2];
    c7=1. + (1. - sin3beta)/fPsi +sin3beta;
    costheta=cos(theta);
    sintheta =sin(theta);
    c1=(sig1+sig2+sig3)/sqrt(3.);
    c2=fR*F(k,fPhi)/sqrt(3.);
    c3= sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c4=2*sqrt2*(F(k,fPhi)-fN)*cos(beta)/c7;
    c5=(1./sqrt(2.))*(sig2-sig3);
    c6=2.*sqrt(2.)*(F(k,fPhi)-fN)*sin(beta)/c7;
    T ddistf2dtheta=(2.*c2*(c1 - k/sqrt(3.) - c2*costheta)*sintheta)/(3.*fK) +
    (-2.*c4*costheta*(c3 - c4*sintheta) - 2.*c6*costheta*(c5 - c6*sintheta))/(2.*fG);
    
    
    c1=pow(sig1/sqrt(3.) + sig2/sqrt(3.) + sig3/sqrt(3.) - (k + fR*(fA - fC*expfBk - k*fPhi)*costheta)/sqrt(3.),2.)/(3.*fK);
    c2=c3;
    c3=2.*sqrt(2.)*(F(k,fPhi)- fN)*sintheta;
    c4=c5;
    c5=c3;
    c6=1./fPsi;
    T ddistf2dbeta =(2.*(c2 - (c3*cosbeta)/(1. + c6*(1. - sin3beta) +sin3beta))*
                        ((c3*cosbeta*(3.*cos3beta - 3.*c6*cos3beta))/pow (1. + c6*(1. - sin3beta) + sin3beta,2.) + (c3*sinbeta)/(1. + c6*(1. - sin3beta) + sin3beta)) +2.*((c5*(3.*cos3beta - 3.*c6*cos3beta)*sinbeta)/pow (1. + c6*(1. - sin3beta) + sin3beta,2.) - (c5*cosbeta)/(1. + c6*(1. - sin3beta) +sin3beta))*(c4 - (c5*sinbeta)/(1 + c6*(1. - sin3beta) +sin3beta)))/(2.*fG);
//    T s1,s2,s3;
//    s1=sig1;
//    s2=sig2;
//    s3=sig3;
//    ddistf2dbeta = (2*((2*sqrt(2)*(fA - fC*exp(fB*k) - fN - k*fPhi)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(beta)*sin(theta))/pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) -
//        (2*sqrt(2)*(fA - fC*exp(fB*k) - fN - k*fPhi)*cos(beta)*sin(theta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)))*
//     (s2/sqrt(2) - s3/sqrt(2) - (2*sqrt(2)*(fA - fC*exp(fB*k) - fN - k*fPhi)*sin(beta)*sin(theta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta))) +
//     2*(sqrt(0.6666666666666666)*s1 - s2/sqrt(6) - s3/sqrt(6) - (2*sqrt(2)*(fA - fC*exp(fB*k) - fN - k*fPhi)*cos(beta)*sin(theta))/(1 + (1 - sin(3*beta))/fPsi + sin(3*beta)))*
//     ((2*sqrt(2)*(fA - fC*exp(fB*k) - fN - k*fPhi)*cos(beta)*(3*cos(3*beta) - (3*cos(3*beta))/fPsi)*sin(theta))/pow(1 + (1 - sin(3*beta))/fPsi + sin(3*beta),2) +
//      (2*sqrt(2)*(fA - fC*exp(fB*k) - fN - k*fPhi)*sin(beta)*sin(theta))/(1 + (1 - sin(3*beta))/fPsi +sin(3*beta))))/(2.*fG);
    
    T resL = ResLF2(pt, theta, beta,k,kprev);
    
    ddistf2[0]=ddistf2dtheta;
    ddistf2[1]=ddistf2dbeta;
    ddistf2[2]=resL;
    
}

void TPZSandlerExtended::D2DistFunc1(const TPZVec<STATE> &pt,STATE xi,STATE beta, TPZFMatrix<STATE> &tangentf1) const
{
    
    STATE sqrt2=sqrt(2);
    STATE sqrt3=sqrt(3);
    STATE expsqrt3fBxi=exp(sqrt3*fB*xi);
    STATE d2distf1dxixi,d2distf1dxibeta;
    STATE d2distf1dbetaxi,d2distf1dbetabeta;
    STATE sinbeta,cos3beta,cosbeta,cos2beta,cos7beta,cos6beta,sin5beta,sin7beta,cos4beta,sin3beta,sig1,sig2,sig3;
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
    sig1 = pt[0];
    sig2 = pt[1];
    sig3 = pt[2];
    d2distf1dxixi=2./(3.*fK) + ((16*pow(cosbeta,2)*pow(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi,2))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) +
                                (16*pow(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi,2)*pow(sinbeta,2))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) +
                                (12*sqrt2*cosbeta*expsqrt3fBxi*pow(fB,2)*fC*(sqrt(0.6666666666666666)*sig1 - sig2/sqrt(6) - sig3/sqrt(6) -
                                                                             (2*sqrt2*cosbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                                (12*sqrt2*expsqrt3fBxi*pow(fB,2)*fC*sinbeta*(sig2/sqrt2 - sig3/sqrt2 - (2*sqrt2*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta))/(2.*fG);
    
    d2distf1dxibeta=((4*sqrt2*cosbeta*(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi)*(3*cos3beta - (3*cos3beta)/fPsi)*
                      (sqrt(0.6666666666666666)*sig1 - sig2/sqrt(6) - sig3/sqrt(6) - (2*sqrt2*cosbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/
                     pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) + (4*sqrt2*(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi)*sinbeta*
                                                                  (sqrt(0.6666666666666666)*sig1 - sig2/sqrt(6) - sig3/sqrt(6) - (2*sqrt2*cosbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/
                     (1 + (1 - sin3beta)/fPsi + sin3beta) - (4*sqrt2*(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi)*sinbeta*
                                                             ((-2*sqrt2*cosbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                                                              (2*sqrt2*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2)))/(1 + (1 - sin3beta)/fPsi + sin3beta) -
                     (4*sqrt2*cosbeta*(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi)*(sig2/sqrt2 - sig3/sqrt2 -
                                                                                  (2*sqrt2*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                     (4*sqrt2*(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi)*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta*
                      (sig2/sqrt2 - sig3/sqrt2 - (2*sqrt2*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) -
                     (4*sqrt2*cosbeta*(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi)*((2*sqrt2*cosbeta*(3*cos3beta - (3*cos3beta)/fPsi)*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/
                                                                                  pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) + (2*sqrt2*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta))/
    (2.*fG);
    
    d2distf1dbetaxi=(2*((2*sqrt2*cosbeta*(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi)*(3*cos3beta - (3*cos3beta)/fPsi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) +
                        (2*sqrt2*(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi)*sinbeta)/(1 + (1 - sin3beta)/fPsi + sin3beta))*
                     (sqrt(0.6666666666666666)*sig1 - sig2/sqrt(6) - sig3/sqrt(6) - (2*sqrt2*cosbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)) -
                     (4*sqrt2*(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi)*sinbeta*((-2*sqrt2*cosbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                                                                                  (2*sqrt2*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2)))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                     2*((-2*sqrt2*cosbeta*(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                        (2*sqrt2*(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi)*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta)/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2))*
                     (sig2/sqrt2 - sig3/sqrt2 - (2*sqrt2*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)) -
                     (4*sqrt2*cosbeta*(-(sqrt3*expsqrt3fBxi*fB*fC) - sqrt3*fPhi)*((2*sqrt2*cosbeta*(3*cos3beta - (3*cos3beta)/fPsi)*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/
                                                                                  pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) + (2*sqrt2*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta)))/(1 + (1 - sin3beta)/fPsi + sin3beta))/
    (2.*fG);
    
    
    d2distf1dbetabeta=(2*(sqrt(0.6666666666666666)*sig1 - sig2/sqrt(6) - sig3/sqrt(6) - (2*sqrt2*cosbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta))*
                       ((-4*sqrt2*cosbeta*pow(3*cos3beta - (3*cos3beta)/fPsi,2)*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,3) +
                        (2*sqrt2*cosbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                        (2*sqrt2*cosbeta*(-9*sin3beta + (9*sin3beta)/fPsi)*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) -
                        (4*sqrt2*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2)) +
                       2*pow((-2*sqrt2*cosbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                             (2*sqrt2*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2),2) +
                       2*pow((2*sqrt2*cosbeta*(3*cos3beta - (3*cos3beta)/fPsi)*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) +
                             (2*sqrt2*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta),2) +
                       2*(sig2/sqrt2 - sig3/sqrt2 - (2*sqrt2*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta))*
                       ((4*sqrt2*cosbeta*(3*cos3beta - (3*cos3beta)/fPsi)*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2) -
                        (4*sqrt2*pow(3*cos3beta - (3*cos3beta)/fPsi,2)*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,3) +
                        (2*sqrt2*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/(1 + (1 - sin3beta)/fPsi + sin3beta) +
                        (2*sqrt2*(-9*sin3beta + (9*sin3beta)/fPsi)*sinbeta*(fA - expsqrt3fBxi*fC - fN - sqrt3*fPhi*xi))/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2)))/(2.*fG);
    
    TPZFMatrix<STATE> JAC(2,2,0);
    tangentf1(0,0)=d2distf1dxixi;
    tangentf1(0,1)=d2distf1dxibeta;
    tangentf1(1,0)=d2distf1dbetaxi;
    tangentf1(1,1)=d2distf1dbetabeta;
    
}

void TPZSandlerExtended::D2DistFunc2(const TPZVec<STATE> &pt,STATE theta,STATE beta,STATE k, TPZFMatrix<STATE> &tangentf2) const
{
  STATE d2distf2dthetatheta,d2distf2dthetabeta,d2distf2dthetak;
  STATE d2distf2dbetatheta,d2distf2dbetabeta,d2distf2dbetak;
  STATE dresktheta,dreskk;
  STATE fR2=fR*fR;
  STATE Gamma=1 + (1 - sin(3.*beta))/fPsi + sin(3.*beta);
  STATE stresbeta=sin(3*beta);
  STATE Gamma2=Gamma*Gamma;
  STATE DGamma=3*cos(3*beta) - (3*cos(3*beta))/fPsi;
  STATE DGamma2=DGamma*DGamma,Gamma3=Gamma*Gamma*Gamma,Gamma4=Gamma2*Gamma2;
  STATE FFL=F(k,fPhi);
  STATE FFL2=FFL*FFL;
  STATE fN2=fN*fN;
  STATE sb=sin(beta),cb=cos(beta),sb2=sb*sb,cb2=cb*cb,ct=cos(theta),st=sin(theta),
  ct2=ct*ct,st2=st*st;
  STATE sig1=pt[0],sig2=pt[1],sig3=pt[2];
  STATE sq2=sqrt(2.),sq3=sqrt(3.);
  
  STATE L=k;
  STATE FN=fN;
  STATE xi=1./sq3*(sig1+sig2+sig3);
  STATE sig2tiu=(sqrt(2./3.)*sig1)-(sig2/sqrt(6.))-(sig3/sqrt(6.));
  STATE expBC=exp(fB*L)*fB*fC + fPhi;
  
  
  d2distf2dthetatheta=(2.*(-(ct2*(FFL2*fG*fR2*Gamma2 - 36.*FFL2*fK*sb2 + 72.*FFL*fK*fN*sb2 - 36.*fK*fN2*sb2)) + 9.*cb*fK*(FFL - fN)*Gamma*sig2tiu*sq2*st +
  st*(FFL2*(fG*fR2*Gamma2 - 36.*fK*sb2)*st - 9.*fK*fN*sb*(Gamma*(sig2 - sig3) + 4.*fN*sb*st) + 9.*FFL*fK*sb*(Gamma*(sig2 - sig3) +
  8.*fN*sb*st)) + 36.*cb2*fK*(FFL - fN)*(ct2*(FFL - fN) + (-FFL + FN)*st2) +
  ct*FFL*fG*fR*Gamma2*(-L + sq3*xi)))/(9.*fG*fK*Gamma2);
  
/*  
  d2distf2dthetabeta=(ct*(FFL - fN)*sq2*(cb*Gamma*(4*DGamma*sig2tiu + sq2*(Gamma*(-sig2 + sig3) + 4*(FFL - fN)*sb*st)) +
  2*(Gamma2*sb*sig2tiu + DGamma*sq2*(Gamma*sb*(sig2 - sig3) - 4*(FFL - fN)*(cb2 + 2*sb2)*st))
  -4*(FFL - FN)*(2*cb*DGamma + Gamma*sb)*sq2*cb*ct))/(fG*Gamma3);
*/

  d2distf2dthetabeta = (2*ct*(FFL - fN)*(8*DGamma*(-FFL + fN)*st + (cb*Gamma - 2*DGamma*sb)*(Gamma*(-sig2 + sig3) + 4*(FFL - fN)*sb*st) - 
       (sq2*Gamma*(2*cb*DGamma + Gamma*sb)*(-((1 + fPsi)*sig2tiu) + 2*sq2*cb*(FFL - FN)*fPsi*st - (-1 + fPsi)*sig2tiu*stresbeta))/
        (1 + fPsi + (-1 + fPsi)*stresbeta)))/(fG*Gamma3);
  
  d2distf2dthetak=(2*(ct*expBC*(2*FFL*fG*fR2*Gamma2*st + 9*fK*(Gamma*(sb*(sig2 - sig3) + cb*sig2tiu*sq2) + 4*(cb2*(-2*FFL + fN + FN) + 2*(-FFL + fN)*sb2)*st)) - fG*fR*Gamma2*st*(FFL - expBC*L + expBC*sq3*xi)))/(9.*fG*fK*Gamma2);
  
  
  
  d2distf2dbetatheta=  -((ct*(FFL - fN)*sq2*(8*cb2*DGamma*(2*FFL - fN - FN)*sq2*st - 2*sb*(Gamma2*sig2tiu + DGamma*sq2*(Gamma*(sig2 - sig3) + 8*(-FFL + fN)*sb*st)) +
  cb*Gamma*(-4*DGamma*sig2tiu + sq2*(Gamma*(sig2 - sig3) + 4*(fN - FN)*sb*st))))/(fG*Gamma3));
  
  STATE one = (4*cb*DGamma*(FFL - fN)*sq2*st)/Gamma2 + (2*(FFL - fN)*sb*sq2*st)/Gamma;
  STATE two = (2*cb*(FFL - fN)*sq2*st)/Gamma - (4*DGamma*(FFL - fN)*sb*sq2*st)/Gamma2;
  d2distf2dbetabeta = (one*one + 
     two*two - 
     (2*(FFL - fN)*sq2*st*(2*cb*(-FFL + FN)*fPsi*sq2*st + sig2tiu*(1 + fPsi - stresbeta + fPsi*stresbeta))*
        (4*DGamma*fPsi*Gamma*Gamma3*sb + cb*(8*DGamma2*fPsi*Gamma*Gamma2 - Gamma3*(fPsi*Gamma2 + 9*Gamma*stresbeta - 9*fPsi*Gamma*stresbeta))))/
      (fPsi*Gamma*Gamma2*Gamma3*(1 + fPsi - stresbeta + fPsi*stresbeta)) + 
     ((FFL - fN)*sq2*st*(Gamma*(sig2 - sig3)*sq2 + 4*(-FFL + fN)*sb*sq2*st)*
        (4*cb*DGamma*fPsi*Gamma*Gamma3 + sb*(-8*DGamma2*fPsi*Gamma*Gamma2 + Gamma3*(fPsi*Gamma2 + 9*Gamma*stresbeta - 9*fPsi*Gamma*stresbeta))))/
      (fPsi*Gamma2*Gamma2*Gamma3))/fG;
/*      
  d2distf2dbetabeta= ((FFL - fN)*sq2*st*(-8*DGamma*Gamma2*sb*sig2tiu*(-1 + stresbeta) - 2*cb*Gamma*(Gamma2*sig2tiu + 8*DGamma2*sig2tiu*(-1 + stresbeta) +
  2*(-(DGamma*sq2*(Gamma*(sig2 - sig3) + 4*(fN - FN)*sb*st)*(-1 + stresbeta)) + Gamma*(-9 + 4*Gamma)*sig2tiu*stresbeta)) +
  sq2*(Gamma3*sb*(-sig2 + sig3) - 8*DGamma2*(Gamma*sb*(sig2 - sig3) + 2*(cb2*(-3*FFL + fN + 2*FN) + 3*(-FFL + fN)*sb2)*st)*(-1 + stresbeta) +
  2*(Gamma*sb*(9*Gamma*(sig2 - sig3 + 2*(FFL - fN)*sb*st) + 4*(Gamma2*(-sig2 + sig3) + 9*(-FFL + fN)*sb*st))*stresbeta -
  2*cb2*st*(fN*Gamma2*(-1 + stresbeta) - 9*FFL*(-2 + Gamma)*Gamma*stresbeta + FN*(Gamma2 + 2*Gamma*(-9 + 4*Gamma)*stresbeta))))))/(fG*Gamma4*(-1 + stresbeta));
*/  
  d2distf2dbetak=(-2*expBC*st*(Gamma2*sb*sig2tiu*sq2 + cb*Gamma*(Gamma*(-sig2 + sig3) + 2*DGamma*sig2tiu*sq2 + 4*(-fN + FN)*sb*st) + 2*DGamma*(Gamma*sb*(sig2 - sig3) + 4*(cb2*(-2*FFL + fN + FN) + 2*(-FFL + fN)*sb2)*st)))/
  (fG*Gamma3);
  
  
  dresktheta= -(FFL*fR*st);
  
  
  dreskk=1 - ct*expBC*fR + 3*fD*exp(fD*(-(FFL*fR) + L))*fK*(1 + expBC*fR)*fW;
  
  tangentf2(0,0)= d2distf2dthetatheta;
  tangentf2(0,1)=d2distf2dthetabeta;
  tangentf2(0,2)=d2distf2dthetak;
  tangentf2(1,0)=d2distf2dbetatheta;
  tangentf2(1,1)=d2distf2dbetabeta;
  tangentf2(1,2)=d2distf2dbetak;
  tangentf2(2,0)=dresktheta;
  tangentf2(2,1)=0;
  tangentf2(2,2)=dreskk;
}


void TPZSandlerExtended::D2DistFunc2new(const TPZVec<STATE> &pt,STATE theta,STATE beta,STATE k, TPZFMatrix<STATE> &tangentf2)const
{
    STATE sqrt3=sqrt(3);
    STATE sqrt2=sqrt(2);
    STATE d2distf2dthetatheta,d2distf2dthetabeta,d2distf2dthetak;
    STATE d2distf2dbetatheta,d2distf2dbetabeta,d2distf2dbetak;
    STATE dresktheta,dreskk;
    STATE sinbeta,cos3beta,sintheta,costheta,cosbeta,cos2beta,cos7beta,cos6beta,sin5beta,sin7beta,cos4beta,sin3beta,cos2theta,c1,c2,c3,c4,c5,c6,c7,c8,c9,gamma,sig1,sig2,sig3;
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
    sig1 = pt[0];
    sig2 = pt[1];
    sig3 = pt[2];
    STATE expfBk=exp(fB*k);
    
    
    c1=(sig1+sig2+sig3)/sqrt(3.);
    c2=(fR*(fA - fC*expfBk - k*fPhi))/sqrt3;
    c3= sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c4=(2*sqrt2*(fA - fC*expfBk - fN - k*fPhi)*cosbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta);
    c5=(1./sqrt(2.))*(sig2-sig3);
    c6=(2*sqrt2*(fA - fC*expfBk - fN - k*fPhi)*sinbeta)/(1 + (1 - sin3beta)/fPsi +sin3beta);
    d2distf2dthetatheta=(2*c2*fG*(3*c1 - sqrt3*k)*costheta + (-6*pow(c2,2)*fG + 9*(pow(c4,2) + pow(c6,2))*fK)*cos2theta + 9*(c3*c4 + c5*c6)*fK*sintheta)/(9.*fG*fK);
    
    
    
    c1=(2*fR*(fA - expfBk*fC - fPhi*k)*(-((k + costheta*fR*(fA - expfBk*fC - fPhi*k))/sqrt3) + sig1/sqrt3 + sig2/sqrt3 + sig3/sqrt3)*sintheta)/(3.*sqrt3*fK);
    c2=4*sqrt2*costheta*(fA - expfBk*fC - fN - fPhi*k);
    c3= sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c4=2*sqrt2*(fA - expfBk*fC - fN - fPhi*k)*sintheta;
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
    c4=2*sqrt2*cosbeta*costheta;
    c5=sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c6=sqrt2*cosbeta*sintheta;
    c7=2*sqrt2*costheta*sinbeta;
    c8=(1./sqrt(2.))*(sig2-sig3);
    c9=sqrt2*sinbeta*costheta;
    gamma=0.5*(1 + (1 - sin3beta)/fPsi +  sin3beta);
    d2distf2dthetak=-(c3*(1 + c2*(-(expfBk*fB*fC) - fPhi))*(fA - expfBk*fC - fPhi*k))/(9.*fK) + (c3*(-(expfBk*fB*fC) - fPhi)*(c1 - (k + c2*(-(fC*expfBk) + fA - fPhi*k))/sqrt3))/(3.*sqrt3*fK) +
    ((c4*c6*(-(expfBk*fB*fC) - fPhi)*(fA - expfBk*fC - fN - fPhi*k))/(gamma*gamma) + (c7*c9*(-(expfBk*fB*fC) - fPhi)*(fA - expfBk*fC - fN - fPhi*k))/(gamma*gamma) -
     (c4*(-(expfBk*fB*fC) - fPhi)*(c5 - (c6*(fA - expfBk*fC - fN - fPhi*k))/gamma))/gamma - (c7*(-(expfBk*fB*fC) - fPhi)*(c8 - (c9*(fA - expfBk*fC - fN - fPhi*k))/gamma))/gamma)/(2.*fG);

    
    
    //////////////
    
    c1=(2*sqrt2*(3*cos3beta - (3*cos3beta)/fPsi)*(fA - expfBk*fC - fN - fPhi*k)*sinbeta)/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2);
    c2=(2*sqrt2*cosbeta*(fA - expfBk*fC - fN - fPhi*k))/(1 + (1 - sin3beta)/fPsi + sin3beta);
    c3=(1./sqrt(2.))*(sig2-sig3);
    c4=(2*sqrt2*(fA - expfBk*fC - fN - fPhi*k)*sinbeta)/(1 + (1 - sin3beta)/fPsi + sin3beta);
    c5=sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c6=(2*sqrt2*(3*cos3beta - (3*cos3beta)/fPsi)*(fA - expfBk*fC - fN - fPhi*k)*cosbeta)/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2);
    c7=(2*sqrt2*sinbeta*(fA - expfBk*fC - fN - fPhi*k))/(1 + (1 - sin3beta)/fPsi + sin3beta);
    d2distf2dbetatheta=(2*(c6*costheta + c7*costheta)*(c5 - c2*sintheta) - 2*c4*costheta*(c1*sintheta - c2*sintheta) + 2*(c1*costheta - c2*costheta)*(c3 - c4*sintheta) - 2*c2*costheta*(c6*sintheta + c7*sintheta))/(2.*fG);

    
    c1=pow(-((k + costheta*fR*(fA - expfBk*fC - fPhi*k))/sqrt3) + sig1/sqrt3 + sig2/sqrt3 + sig3/sqrt3,2)/(3.*fK);
    c2=sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c3=2*sqrt2*sintheta*(fA - expfBk*fC - fN - k*fPhi);
    c4=(1./sqrt(2.))*(sig2-sig3);
    c5=2*sqrt2*(fA -expfBk*fC - fN - fPhi*k)*sintheta;
    c6=1/fPsi;
    d2distf2dbetabeta=(2*pow(c5,2)*pow(cosbeta,2)*pow(1 + c6 + 2*(-1 + c6)*sin3beta - 6*(-1 + c6)*sinbeta,2) + 2*pow(c3*(-1 + c6)*(2*cos2beta + cos4beta) - c3*(1 + c6)*sinbeta,2) +2*c3*cosbeta*(-(c2*(1 + c6)) + c3*cosbeta + c2*(-1 + c6)*sin3beta)*(15 - 34*c6 + 15*pow(c6,2) - 6*pow(-1 + c6,2)*
                                                                                                                                                                                                                                                        cos2beta + 6*pow(-1 + c6,2)*cos4beta + 2*pow(-1 + c6,2)*cos6beta -13*(-1 + pow(c6,2))*sin3beta + 12*(-1 + pow(c6,2))*sinbeta) + c5*(-(c4*(1 + c6)) + c4*(-1 + c6)*sin3beta + c5*sinbeta)*((1 - pow(c6,2))*cos2beta + 13*(-1 + pow(c6,2))*cos4beta + 2*(pow(-1 + c6,2)*(-4*sin5beta + sin7beta) + 4*(3 + c6*(-7 + 3*c6))*sinbeta)))/(2.*fG*pow(1 + c6 - (-1 + c6)*sin3beta,4));

    
    
    c1=(2.*sqrt(2.)*(3*cos3beta - (3*cos3beta)/fPsi)*sinbeta*sintheta)/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2);
    c2=(2*sqrt2*cosbeta*sintheta)/(1 + (1 - sin3beta)/fPsi + sin3beta);
    c3=(2*sqrt2*sinbeta*sintheta)/(1 + (1 - sin3beta)/fPsi + sin3beta);
    c4=(1./sqrt(2.))*(sig2-sig3);
    c5=c2;
    c6=sqrt(2./3.)*sig1 - sig2/sqrt(6.) -sig3/sqrt(6.);
    c7=(2*sqrt2*(3*cos3beta - (3*cos3beta)/fPsi)*cosbeta*sintheta)/pow(1 + (1 - sin3beta)/fPsi + sin3beta,2);
    c8=c3;
    d2distf2dbetak=(-2*c3*(-(expfBk*fB*fC) - fPhi)*(c1*(fA - expfBk*fC - fN - fPhi*k) - c2*(fA - expfBk*fC - fN - fPhi*k)) +
                    2*(c1*(-(expfBk*fB*fC) - fPhi) - c2*(-(expfBk*fB*fC) - fPhi))*(c4 - c3*(fA - expfBk*fC - fN - fPhi*k)) +
                    2*(c7*(-(expfBk*fB*fC) - fPhi) + c8*(-(expfBk*fB*fC) - fPhi))*(c6 - c5*(fA - expfBk*fC - fN - fPhi*k)) -
                    2*c5*(-(expfBk*fB*fC) - fPhi)*(c7*(fA - expfBk*fC - fN - fPhi*k) + c8*(fA - expfBk*fC - fN - fPhi*k)))/(2.*fG);
    
    
    dreskk=1. + costheta*(-(expfBk*fB*fC) - fPhi)*fR +
    3.*exp(fD*(-((fA - expfBk*fC)*fR) + k))*fD*fK*
    (1. + expfBk*fB*fC*fR)*fW;

    
    dresktheta=-fR*(fA - fC * expfBk - k*fPhi)* sintheta;

    tangentf2(0,0)= d2distf2dthetatheta;
    tangentf2(0,1)=d2distf2dthetabeta;
    tangentf2(0,2)=d2distf2dthetak;
    tangentf2(1,0)=d2distf2dbetatheta;
    tangentf2(1,1)=d2distf2dbetabeta;
    tangentf2(1,2)=d2distf2dbetak;
    tangentf2(2,0)=dresktheta;
    tangentf2(2,1)=0;
    tangentf2(2,2)=dreskk;
}



void TPZSandlerExtended::YieldFunction(const TPZVec<STATE> &sigma, STATE kprev, TPZVec<STATE> &yield) const
{
    
    STATE II1,JJ2,JJ3,ggamma,temp1,temp3,f2,sqrtj2,f1,beta;
    TPZManVector<STATE,3> cylstress(3);
    FromPrincipalToHWCyl(sigma,cylstress);
    beta=cylstress[2];
    TPZTensor<STATE> sigten;
    sigten.XX() = sigma[0];
    sigten.YY() = sigma[1];
    sigten.ZZ() = sigma[2];
    II1 = sigten.I1();
    JJ2 = sigten.J2();
    JJ3 = sigten.J3();
    if (JJ2<1.e-6) {
        JJ2=1.e-6;
    }
    sqrtj2=sqrt(JJ2);
    ggamma = 0.5*(1. + (1. - sin(3.*beta))/fPsi + sin(3.*beta));
    
    temp1=(-II1+kprev)/(-fR*F(kprev,fPhi));
    temp3=(ggamma*sqrtj2)/(F(kprev,fPhi));
    
    f1=sqrtj2-F(II1,fPhi);
    f2=temp1*temp1+temp3*temp3-1;

    yield[0]=f1;
    yield[1]=f2;
    
}

void TPZSandlerExtended::ProjectF1(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigproj, STATE &kproj) const
{
//#ifdef LOG4CXX
//    {
//        std::stringstream outfile;
//        outfile << "\n projection over F1 " <<endl;
//        LOGPZ_DEBUG(logger,outfile.str());
//        
//    }
//#endif
    
    //ofstream convergenceF1("convergenceF1.txt");
    STATE xi,resnorm,beta,distxi,distnew;
    distxi=1.e8;
    for (STATE xiguess=-M_PI; xiguess <= M_PI; xiguess += M_PI/20.)
    {
        for(STATE betaguess=0;betaguess<=2*M_PI;betaguess+=M_PI/20.)
        {
            distnew=DistF1(sigmatrial, xiguess,betaguess);
            if (fabs(distnew) < fabs(distxi))
            {
                xi = xiguess;
                beta=betaguess;
                distxi = distnew;
            }
        }
    }
    
    resnorm=1.;
    long counter=1;
    TPZFMatrix<STATE> xn1(2,1,0.),xn(2,1,0.),jac,invjac,sol(2,1,0.),fxn(2,1,0.),diff(2,1,0.);
    xn(0,0)=xi;
    xn(1,0)=beta;
    while (resnorm > 10e-12 && counter < 30)
    {
        
        TPZFNMatrix<4,STATE> jac(2,2);
        D2DistFunc1(sigmatrial, xn[0],xn[1],jac);
        DDistFunc1(sigmatrial, xn[0],xn[1],fxn);
        sol = fxn;
        jac.Solve_LU(&sol);
        xn1=xn-sol;
        diff=xn1-xn;
        resnorm=Norm(diff);
        //convergenceF1 << counter << " "<<resnorm <<endl;
        
//#ifdef LOG4CXX
//        {
//            std::stringstream outfile;//("convergencF1.txt");
//            outfile<< "\n" <<counter << " "<<resnorm <<endl;
//            LOGPZ_DEBUG(logger,outfile.str());
//            
//        }
//#endif
        xn=xn1;
        counter++;

        
    }
    
    TPZManVector<STATE,3> sigprojcyl(3);
    F1Cyl(xn[0], xn[1], sigprojcyl);
    FromHWCylToPrincipal(sigprojcyl, sigproj);
    
  
    STATE kguess = sigproj[0]+sigproj[1]+sigproj[2];
    //Firstk(kguess, kproj);
    STATE resl = ResLF1(sigmatrial, sigproj, kguess, kprev);
    
    int count =0;
    while (fabs(resl) > 1.e-14 && count < 30) {
        STATE dresl = DResLF1(sigmatrial, sigproj, kguess, kprev);
        kguess -= resl/dresl;
        resl = ResLF1(sigmatrial, sigproj, kguess, kprev);
        count++;
    }
    
    kproj = kguess;
   
}

void TPZSandlerExtended::ProjectF2(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigproj, STATE &kproj) const
{
//#ifdef LOG4CXX
//    {
//        std::stringstream outfile;
//        outfile << "\n projection over F2 " <<endl;
//        LOGPZ_DEBUG(logger,outfile.str());
//        
//    }
//#endif

    STATE theta,beta=0.,distnew;
    STATE resnorm,disttheta;
    disttheta=1.e8;
    TPZManVector<STATE,3> vectempcyl(3);
    FromPrincipalToHWCyl(sigmatrial, vectempcyl);
    STATE betaguess=vectempcyl[2];
    for (STATE thetaguess=M_PI/2.; thetaguess <= M_PI; thetaguess += M_PI/20.) {
            distnew=DistF2(sigmatrial,thetaguess,betaguess,kprev);
            if (fabs(distnew) < fabs(disttheta)) {
                theta = thetaguess;
                beta=betaguess;
                disttheta = distnew;
            }
    }
    
    resnorm=1;
    int counter=1;
    TPZFNMatrix<3,STATE> xn1(3,1,0.),xn(3,1,0.),sol(3,1,0.),fxn(3,1,0.),diff(3,1,0.);
    xn(0,0)=theta;
    xn(1,0)=beta;
    xn(2,0)=kprev;
    while (resnorm > 10.e-15 && counter < 30)
    {

        TPZFNMatrix<9,STATE> jac(3,3);
        D2DistFunc2(sigmatrial, xn(0),xn(1),xn(2),jac);
        TPZManVector<STATE> fxnvec(3);
        DDistFunc2(sigmatrial, xn(0),xn(1),xn(2),kprev,fxnvec);
//#ifdef LOG4CXX
//        {
//            std::stringstream outfile;//("convergencF1.txt");
//            outfile<< "\n" <<counter << " "<<resnorm <<endl;
//            jac.Print(outfile);
//            outfile<< "\n xn " << " "<<fxnvec <<endl;
//            //outfile<< "\n res " << " "<<fxnvec <<endl;
//            
//            
//            LOGPZ_DEBUG(logger,outfile.str());
//            
//        }
//#endif
        for(int k=0; k<3; k++) sol(k,0) = fxnvec[k];
        jac.Solve_LU(&sol);
        xn1=xn-sol;
        diff=xn1-xn;
        resnorm=Norm(diff);

        xn=xn1;
        counter++;
        
    }
    if(counter == 30) cout << "resnorm = " << resnorm << std::endl;
//    cout<< "\n resnorm = "<<resnorm <<endl;
//    cout<< "\n counter = "<<counter <<endl;
//    cout<< "\n k = "<<(xn1(2,0)-kprev) <<endl;
    STATE thetasol,betasol,ksol;
    
    thetasol=xn1(0);
    betasol=xn1(1);
    ksol=xn1(2);
    kproj=ksol;

    TPZManVector<STATE,3> f2cyl(3);
    F2Cyl(thetasol, betasol, ksol, f2cyl);
    FromHWCylToPrincipal(f2cyl,sigproj);
    
}

void TPZSandlerExtended::ProjectRing(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigproj,STATE &kproj) const
{
    STATE theta,beta=0.,distnew;
    STATE resnorm,disttheta;
    disttheta=1.e8;
    
    for (STATE betaguess=0; betaguess <= 2*M_PI; betaguess += M_PI/20.) {
        distnew=DistF2(sigmatrial,M_PI/2,betaguess,kprev);
        if (fabs(distnew) < fabs(disttheta)) {
            theta = M_PI/2;
            beta=betaguess;
            disttheta = distnew;
        }
    }
    resnorm=1;
    long counter=1;
    TPZFMatrix<STATE> xn1(3,1,0.),xn(3,1,0.),sol(3,1,0.),fxn(3,1,0.),diff(3,1,0.);
    xn(0,0)=M_PI/2;
    xn(1,0)=beta;
    xn(2,0)=kprev;
    while (resnorm > 10e-15 && counter < 30)
    {
        TPZFNMatrix<9,STATE> jac(3,3);
        D2DistFunc2(sigmatrial,xn[0],xn[1],xn[2],jac);
        TPZManVector<STATE> fxnvec(3);
        DDistFunc2new(sigmatrial, xn(0),xn(1),xn(2),kprev,fxnvec);
        for(int k=0; k<3; k++) fxn(k,0) = fxnvec[k];
        for (int i=0; i<3; i++) {
            jac(i,0) = 0.;
            jac(0,i) = 0.;
        }
        jac(0,0) = 1.;
        fxn(0,0) = 0.;
        sol = fxn;
        jac.Solve_LU(&sol);
        
        xn1(0,0)=xn(0,0);
        xn1(1,0)=xn(1,0)-sol(1,0);
        xn1(2,0)=xn(2,0)-sol(2,0);
        
        diff=xn1-xn;
        resnorm=Norm(diff);
        xn=xn1;
        counter++;

    }
//    cout<< "\n resnorm = "<<resnorm <<endl;
//    cout<< "\n counter = "<<xn1 <<endl;
//    cout<< "\n k = "<<xn1[2] <<endl;
    STATE thetasol,betasol,ksol;
    
    thetasol=xn1[0];
    betasol=xn1[1];
    ksol=xn1[2];

    TPZManVector<STATE,3> f2cyl(3);
    F2Cyl(thetasol, betasol, ksol, f2cyl);
    FromHWCylToPrincipal(f2cyl,sigproj);

    kproj = ksol;
    
}

void TPZSandlerExtended::ProjectBetaConstF2(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigproj,STATE &kproj) const
{
    //#ifdef LOG4CXX
    //    {
    //        std::stringstream outfile;
    //        outfile << "\n projection over F2 " <<endl;
    //        LOGPZ_DEBUG(logger,outfile.str());
    //
    //    }
    //#endif
    
    STATE theta,beta=0.,distnew;
    STATE resnorm,disttheta;
    disttheta=1.e8;
    STATE betaconst =0;
    for (STATE thetaguess=0; thetaguess <= M_PI; thetaguess += M_PI/20.) {
            distnew=DistF2(sigmatrial,thetaguess,betaconst,kprev);
            if (fabs(distnew) < fabs(disttheta)) {
                theta = thetaguess;
                beta=betaconst;
                disttheta = distnew;
            }
    }
    
    resnorm=1;
    int counter=1;
    TPZFNMatrix<3,STATE> xn1(3,1,0.),xn(3,1,0.),sol(3,1,0.),fxn(3,1,0.),diff(3,1,0.);
    xn(0,0)=theta;
    xn(1,0)=beta;
    xn(2,0)=kprev;
    while (resnorm > 10.e-15 && counter < 30)
    {
        TPZFNMatrix<9,STATE> jac(3,3);
        D2DistFunc2(sigmatrial, xn(0),xn(1),xn(2),jac);
        TPZManVector<STATE> fxnvec(3);
        DDistFunc2(sigmatrial, xn(0),xn(1),xn(2),kprev,fxnvec);
        for(int k=0; k<3; k++) fxn(k,0) = fxnvec[k];
        for (int i=0; i<3; i++) {
            jac(i,1) = 0.;
            jac(1,i) = 0.;
        }
        jac(1,1) = 1.;
        fxn(1,0) = 0.;
        sol = fxn;
        jac.Solve_LU(&sol);
        
        xn1(0,0)=xn(0,0)-sol(0,0);
        xn1(1,0)=xn(1,0);
        xn1(2,0)=xn(2,0)-sol(2,0);
        
        diff=xn1-xn;
        resnorm=Norm(diff);
        xn=xn1;
        counter++;
        
    }
    if(counter == 30) cout << "resnorm = " << resnorm << std::endl;
    STATE thetasol,betasol,ksol;
    
    thetasol=xn1(0);
    betasol=xn1(1);
    ksol=xn1(2);
    kproj=ksol;
    
    TPZManVector<STATE,3> f2cyl(3);
    F2Cyl(thetasol, betasol, ksol, f2cyl);
    FromHWCylToPrincipal(f2cyl,sigproj);
}

void TPZSandlerExtended::ComputeI1(TPZVec<STATE> stress, STATE &I1)const
{
    STATE sig1,sig2,sig3;
    sig1 = stress[0];
    sig2 = stress[1];
    sig3 = stress[2];
    I1=sig1+sig2+sig3;
    
}

void TPZSandlerExtended::ComputeJ2(TPZVec<STATE> stress,STATE &J2)const
{
    STATE sig1,sig2,sig3;
    sig1 = stress[0];
    sig2 = stress[1];
    sig3 = stress[2];
    J2=(2.*sig1*sig2 +pow(sig1 + (-sig1 - sig2 - sig3)/3.,2.) +
     pow(sig2 + (-sig1 - sig2 - sig3)/3.,2.) + 2*sig1*sig3 + 2.*sig2*sig3 +
        pow((-sig1 - sig2 - sig3)/3. + sig3,2.))/2.;
}


void TPZSandlerExtended::ApplyStrainComputeElasticStress(TPZVec<STATE> &strain,TPZVec<STATE> &stress)const
{
    STATE sig1,sig2,sig3,s1,s2,s3;
    sig1 = strain[0];
    sig2 = strain[1];
    sig3 = strain[2];
    
    s1=sig1-(1./3.)*(sig1+sig2+sig3);
    s2=sig2-(1./3.)*(sig1+sig2+sig3);
    s3=sig3-(1./3.)*(sig1+sig2+sig3);
    
    stress[0]=s1*(2*fG)+fK*(sig1+sig2+sig3);
    stress[1]=s2*(2*fG)+fK*(sig1+sig2+sig3);
    stress[2]=s3*(2*fG)+fK*(sig1+sig2+sig3);
}

void TPZSandlerExtended::ApplyStressComputeElasticStrain(TPZVec<STATE> &stress,TPZVec<STATE> &strain)const
{
    STATE sig1,sig2,sig3,s1,s2,s3;
    sig1 = stress[0];
    sig2 = stress[1];
    sig3 = stress[2];
    
    s1=sig1-(1./3.)*(sig1+sig2+sig3);
    s2=sig2-(1./3.)*(sig1+sig2+sig3);
    s3=sig3-(1./3.)*(sig1+sig2+sig3);
    
    strain[0]=s1/(2.*fG)+(sig1+sig2+sig3)/(9.*fK);
    strain[1]=s2/(2.*fG)+(sig1+sig2+sig3)/(9.*fK);
    strain[2]=s3/(2.*fG)+(sig1+sig2+sig3)/(9.*fK);
    
}

/**
 * Imposes the specified strain tensor and returns the correspondent stress state.
 *
 * @param[in] epsTotal Imposed total strain tensor
 * @param[out] sigma Resultant stress
 */
void TPZSandlerExtended::ApplyStrainComputeSigma(TPZVec<STATE> &epst,TPZVec<STATE> &epsp,STATE & kprev,TPZVec<STATE> &epspnext,TPZVec<STATE> &stressnext,STATE & knext) const
{

    STATE I1tr,I1proj;
    
    TPZManVector<STATE,3> stresstrial(3),yield(2),deltastress(3),delepsp(3),epsT(epst);
    epsT[0]-=epsp[0];
    epsT[1]-=epsp[1];
    epsT[2]-=epsp[2];
    ApplyStrainComputeElasticStress(stresstrial,epsT);
    YieldFunction(stresstrial, kprev, yield);
    ComputeI1(stresstrial,I1tr);
    
    
    if ((yield[1]<=0 && I1tr<=kprev)||(yield[0]<=0 && I1tr>kprev))
    {

        epspnext=epsp;
        knext=kprev;
        stressnext=stresstrial;
        cout<<"\n elastic "<<endl;
        
    }
    else
    {
        cout<<"\n plastic "<<endl;
        if (yield[1]>0 && I1tr<kprev)
        {
            cout<<"\n F2 "<<endl;
            ProjectF2(stresstrial,kprev,stressnext,knext);
        }
        else
        {
            cout<<"\n F1 "<<endl;
            ProjectF1(stresstrial,kprev,stressnext,knext);
            ComputeI1(stressnext,I1proj);
            if (I1proj<knext)
            {
                cout<<"\n Ring "<<endl;
                ProjectRing(stresstrial,kprev,stressnext,knext);
            }
        }
        
        for (int i=0;i<3;i++)
        {
            deltastress[i]=stresstrial[i]-stressnext[i];
        }
        
        ApplyStressComputeElasticStrain(deltastress,delepsp);
        
        for (int i=0;i<3;i++)
        {
            epspnext[i]=epsp[i]+delepsp[i];
        }
    

    }
    
    
    

}

void TPZSandlerExtended::ProjectSigma(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &sigproj,STATE &kproj) const
{
    STATE I1;
    //Firstk(epspv,k0);
    TPZManVector<STATE,2> yield(2);
    I1 = sigtrial[0]+sigtrial[1]+sigtrial[2];
    
    YieldFunction(sigtrial,kprev,yield);
    
    if (I1<kprev)
    {
        if (yield[1]>0.)
        {
            ProjectF2(sigtrial,kprev,sigproj,kproj);
#ifdef DEBUG
            {
                TPZManVector<STATE> cyltr(3), cylproj(3);
                FromPrincipalToHWCyl(sigtrial, cyltr);
                FromPrincipalToHWCyl(sigproj, cylproj);
                std::cout << "cyltr " << cyltr << " cylpr " << cylproj << std::endl;
            }
#endif
        }
        else
        {
            sigproj = sigtrial;
            kproj = kprev;
        }
    }
    else
    {
        if (yield[0]>0.)
        {
            ProjectF1(sigtrial,kprev,sigproj,kproj);
            // this is a wrong condition!!
            I1 = 0.;
            for (int i=0; i<3; i++) {
                I1 += sigproj[i];
            }
            if (I1<kproj)
            {
                ProjectRing(sigtrial,kprev,sigproj,kproj);
            }
            
        }
        else
        {
            // elastic behaviour
            sigproj = sigtrial;
            kproj=kprev;
        }
    }

}

void TPZSandlerExtended::ProjectSigmaDep(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &sigproj,STATE &kproj, TPZFMatrix<STATE> &GradSigma) const
{
    STATE I1;
    //Firstk(epspv,k0);
    TPZManVector<STATE,2> yield(2);
    I1 = sigtrial[0]+sigtrial[1]+sigtrial[2];
    
    YieldFunction(sigtrial,kprev,yield);
    bool treeEigEqual = false;
    STATE tol=1.e-8;
    if (fabs(sigtrial[0]-sigtrial[1])<tol && fabs(sigtrial[1]-sigtrial[2])<tol) {
        treeEigEqual=true;
    }
    
    
    if (I1<kprev)
    {
        if (yield[1]>0. && treeEigEqual==false)
        {
            ProjectF2(sigtrial,kprev,sigproj,kproj);
            // we can compute the tangent matrix
            TPZFNMatrix<9,STATE> dbetadsigtrial(3,3), jacF2(3,3), DF2cart(3,3);
            STATE theta,beta;
            SurfaceParamF2(sigproj, kproj, theta, beta);
            GradF2SigmaTrial(sigtrial, theta, beta, kproj, kprev, dbetadsigtrial);
            D2DistFunc2(sigtrial, theta, beta, kproj, jacF2);
            jacF2.Solve_LU(&dbetadsigtrial);
            DF2Cart(theta, beta, kproj, DF2cart);
            DF2cart.Multiply(dbetadsigtrial, GradSigma);
            GradSigma *= -1.;
        }
        else if (yield[1]>0. && treeEigEqual==true)
        {
            ProjectBetaConstF2(sigtrial,kprev,sigproj,kproj);
            // we can compute the tangent matrix
            TPZFNMatrix<9,STATE> dbetadsigtrial(3,3), jacF2(3,3), DF2cart(3,3);
            STATE theta,beta;
            SurfaceParamF2(sigproj, kproj, theta, beta);
            beta=0;
//#ifdef DEBUG
//            if(fabs(sigproj[1]) > tol)
//            {
//                DebugStop();
//            }
//#endif
            GradF2SigmaTrial(sigtrial, theta, beta, kproj, kprev, dbetadsigtrial);
            for(int i=0; i<3; i++) dbetadsigtrial(1,i) = 0.;
            D2DistFunc2(sigtrial, theta, beta, kproj, jacF2);
            for (int i=0; i<3; i++) {
                jacF2(i,1) = 0.;
                jacF2(1,i) = 0.;
            }
            jacF2(1,1) = 1.;
            jacF2.Solve_LU(&dbetadsigtrial);
            DF2Cart(theta, beta, kproj, DF2cart);
            for(int i=0; i<3; i++) DF2cart(i,1) = 0.;
            DF2cart.Multiply(dbetadsigtrial, GradSigma);
            GradSigma *= -1.;
            
        }
        else
        {
            sigproj = sigtrial;
            kproj = kprev;
            GradSigma.Identity();
        }
    }
    else
    {
        if (yield[0]>0.)
        {
            ProjectF1(sigtrial,kprev,sigproj,kproj);
            
            I1 = 0.;
            for (int i=0; i<3; i++) {
                I1 += sigproj[i];
            }
            if (I1<kproj)
            {
                ProjectRing(sigtrial,kprev,sigproj,kproj);
                
                // we can compute the tangent matrix
                TPZFNMatrix<9,STATE> dbetadsigtrial(3,3), jacF2(3,3), DF2cart(3,3);
                STATE theta,beta;
                SurfaceParamF2(sigproj, kproj, theta, beta);
#ifdef DEBUG
                if(fabs(theta - M_PI_2) > 1.e-8)
                {
                    DebugStop();
                }
#endif
                GradF2SigmaTrial(sigtrial, theta, beta, kproj, kprev, dbetadsigtrial);
                for(int i=0; i<3; i++) dbetadsigtrial(0,i) = 0.;
                D2DistFunc2(sigtrial, theta, beta, kproj, jacF2);
                for (int i=0; i<3; i++) {
                    jacF2(i,0) = 0.;
                    jacF2(0,1) = 0.;
                }
                jacF2(0,0) = 1.;
                jacF2.Solve_LU(&dbetadsigtrial);
                DF2Cart(theta, beta, kproj, DF2cart);
                for(int i=0; i<3; i++) DF2cart(i,0) = 0.;
                DF2cart.Multiply(dbetadsigtrial, GradSigma);
                GradSigma *= -1.;
            }
            else
            {
                // we can compute the tangent matrix
                TPZFNMatrix<9,STATE> dbetadsigtrial(2,3), jacF1(2,2), DF1cart(3,2);
                STATE xi,beta;
                SurfaceParamF1(sigproj, xi, beta);
                GradF1SigmaTrial(sigtrial, xi, beta, dbetadsigtrial);
                D2DistFunc1(sigtrial, xi, beta, jacF1);
                jacF1.Solve_LU(&dbetadsigtrial);
                DF1Cart(xi, beta, DF1cart);
                DF1cart.Multiply(dbetadsigtrial, GradSigma);
                GradSigma *= -1.;
            }
            
        }
        else
        {
            // elastic behaviour
            sigproj = sigtrial;
            kproj=kprev;
            GradSigma.Identity();
        }
    }

}


#define Cos cos
#define Sin sin
#define Sqrt sqrt


/// Compute the derivative of the residual with respect to sigtrial
void TPZSandlerExtended::GradF1SigmaTrial(const TPZVec<STATE> &sigtrial, STATE xi, STATE beta, TPZFMatrix<STATE> &deriv) const
{
    STATE s3beta = sin(3.*beta);
    STATE c3beta = cos(3.*beta);
    STATE Gamma = (1+s3beta)+(1.-s3beta)/fPsi;
    STATE DGamma = 3.*c3beta*(1.-1./fPsi);
    STATE SQR3 = sqrt(3.);
    STATE FFI= fA - fPhi*SQR3*xi - fC*exp(fB*SQR3*xi);
    STATE NN=fN;
    deriv.Redim(2,3);
    deriv(0,0) = (-2*(fG*Gamma - 6*fK*(SQR3*fA*fB - SQR3*fB*FFI + SQR3*fPhi - 3*fB*fPhi*xi)*Cos(beta)))/(3.*SQR3*fG*fK*Gamma);
    deriv(1,0) = (4*(FFI - NN)*(DGamma*Cos(beta) + Gamma*Sin(beta)))/(SQR3*fG*Gamma*Gamma);
    
    deriv(0,1) = (-2*(SQR3*fG*Gamma + 9*fK*(fA*fB + fPhi - fB*(FFI + SQR3*fPhi*xi))*Cos(beta) - 
                      9*fK*(SQR3*fA*fB + SQR3*fPhi - fB*(SQR3*FFI + 3*fPhi*xi))*Sin(beta)))/(9.*fG*fK*Gamma);
    deriv(1,1) = (-2*(FFI - NN)*((SQR3*DGamma + 3*Gamma)*Cos(beta) + (-3*DGamma + SQR3*Gamma)*Sin(beta)))/(3.*fG*Gamma*Gamma);
    
    deriv(0,2) = (-2*(SQR3*fG*Gamma + 9*fK*(fA*fB + fPhi - fB*(FFI + SQR3*fPhi*xi))*Cos(beta) + 
                      9*fK*(SQR3*fA*fB + SQR3*fPhi - fB*(SQR3*FFI + 3*fPhi*xi))*Sin(beta)))/(9.*fG*fK*Gamma);
    deriv(1,2) = (-2*(FFI - NN)*((SQR3*DGamma - 3*Gamma)*Cos(beta) + (3*DGamma + SQR3*Gamma)*Sin(beta)))/(3.*fG*Gamma*Gamma);
}

/// Compute the derivative of the F2 residual with respecto do sigtrial
void TPZSandlerExtended::GradF2SigmaTrial(const TPZVec<STATE> &sigtrial, STATE theta, STATE beta, STATE k, STATE kprev, TPZFMatrix<STATE> &deriv) const
{
    STATE s3beta = sin(3.*beta);
    STATE c3beta = cos(3.*beta);
    STATE Gamma = (1+s3beta)+(1.-s3beta)/fPsi;
    STATE DGamma = 3.*c3beta*(1.-1./fPsi);
    STATE SQR3 = sqrt(3.);
    STATE FFK= fA - fPhi*k - fC*exp(fB*k);
    deriv.Redim(3,3);
    deriv(0,0) = (-4*(FFK - fN)*Cos(beta)*Cos(theta))/(SQR3*fG*Gamma) + (2*FFK*fR*Sin(theta))/(9.*fK);
    deriv(1,0) = (4*(FFK - fN)*(DGamma*Cos(beta) + Gamma*Sin(beta))*Sin(theta))/(SQR3*fG*Gamma*Gamma);
    deriv(2,0) = -1.;
    
    deriv(0,1) = (2*(3*SQR3*fK*(FFK - fN)*Cos(beta)*Cos(theta) - 9*fK*(FFK - fN)*Cos(theta)*Sin(beta) + FFK*fG*fR*Gamma*Sin(theta)))/(9.*fG*fK*Gamma);
    deriv(1,1) = (-2*(FFK - fN)*((SQR3*DGamma + 3*Gamma)*Cos(beta) + (-3*DGamma + SQR3*Gamma)*Sin(beta))*Sin(theta))/(3.*fG*Gamma*Gamma);
    deriv(2,1) = -1.;
    
    deriv(0,2) = (2*(3*SQR3*fK*(FFK - fN)*Cos(beta)*Cos(theta) + 9*fK*(FFK - fN)*Cos(theta)*Sin(beta) + FFK*fG*fR*Gamma*Sin(theta)))/(9.*fG*fK*Gamma);
    deriv(1,2) = (-2*(FFK - fN)*((SQR3*DGamma - 3*Gamma)*Cos(beta) + (3*DGamma + SQR3*Gamma)*Sin(beta))*Sin(theta))/(3.*fG*Gamma*Gamma);
    deriv(2,2) = -1.;
    
}

void TPZSandlerExtended::TaylorCheckDistF1(const TPZVec<STATE> &sigmatrial, STATE xi, STATE beta, TPZVec<STATE> &xnorm,
                       TPZVec<STATE> &errnorm) const
{
    STATE deltaxi = 0.4;
    STATE deltabeta = 0.05;
    STATE dist0 = DistF1(sigmatrial, xi, beta);
    TPZFNMatrix<4,STATE> jac(2,1);
    DDistFunc1(sigmatrial, xi, beta, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i=1; i<=10; i++) {
        STATE diffxi = deltaxi*i/10.;
        STATE xinext = xi+diffxi;
        STATE diffbeta = deltabeta*i/10.;
        STATE betanext = beta + diffbeta;
        STATE distnext = DistF1(sigmatrial, xinext, betanext);
        STATE distguess = dist0+jac(0,0)*diffxi+jac(1,0)*diffbeta;
        xnorm[i-1] = sqrt(diffxi*diffxi+diffbeta*diffbeta);
        errnorm[i-1] = fabs(distnext-distguess);
    }
    
}
void TPZSandlerExtended::TaylorCheckDDistF1(const TPZVec<STATE> &sigmatrial, STATE xi, STATE beta, TPZVec<STATE> &xnorm,
                        TPZVec<STATE> &errnorm) const
{
    STATE deltaxi = 0.4;
    STATE deltabeta = 0.05;
    TPZFNMatrix<2,STATE> res0(2,1),resid(2,1),residguess(2,1),diff(2,1);
    TPZFNMatrix<4,STATE> jac(2,2);
    DDistFunc1(sigmatrial, xi, beta, res0);
    D2DistFunc1(sigmatrial, xi, beta, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i=1; i<=10; i++) {
        STATE diffxi = deltaxi*i/10.;
        STATE xinext = xi+diffxi;
        STATE diffbeta = deltabeta*i/10.;
        diff(0) = diffxi;
        diff(1) = diffbeta;
        STATE betanext = beta + diffbeta;
        DDistFunc1(sigmatrial, xinext, betanext,resid);
        jac.Multiply(diff, residguess);
        residguess += res0;
        xnorm[i-1] = Norm(diff);
        errnorm[i-1] = Norm(resid-residguess);
    }
    
}
void TPZSandlerExtended::TaylorCheckDDistF1DSigtrial(const TPZVec<STATE> &sigmatrial, STATE xi, STATE beta, TPZVec<STATE> &xnorm,
                                 TPZVec<STATE> &errnorm) const
{
    TPZManVector<STATE,3> deltasigma(3,0.3);
    deltasigma[1] *= -1.;
    deltasigma[2] *= 0.7;
    
    TPZFNMatrix<2,STATE> res0(2,1),resid(2,1),residguess(2,1),diff(3,1);
    TPZFNMatrix<6,STATE> jac(2,3);
    DDistFunc1(sigmatrial, xi, beta, res0);
    GradF1SigmaTrial(sigmatrial, xi, beta, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i=1; i<=10; i++) {
        for(int j=0; j<3; j++) diff(j) = deltasigma[j]*i/10.;
        TPZManVector<STATE,3> sigmanext(3);
        for(int j=0; j<3; j++) sigmanext[j] = sigmatrial[j]+diff(j);
        DDistFunc1(sigmanext, xi, beta,resid);
        jac.Multiply(diff, residguess);
        residguess += res0;
        xnorm[i-1] = Norm(diff);
        errnorm[i-1] = Norm(resid-residguess);
    }
    
}

void TPZSandlerExtended::ConvergenceRate(TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm, TPZVec<STATE> &convergence)
{
    convergence.resize(xnorm.size()-1);
    for (int i=1; i<xnorm.size(); i++) {
        convergence[i-1] = log(errnorm[i]/errnorm[i-1])/log(xnorm[i]/xnorm[i-1]);
    }
}

void TPZSandlerExtended::TaylorCheckDistF2(const TPZVec<STATE> &sigmatrial, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &xnorm,
                                           TPZVec<STATE> &errnorm) const
{
    STATE deltatheta = 0.4;
    STATE deltabeta = 0.05;
    STATE deltak = 0.;
    STATE dist0 = DistF2(sigmatrial, theta, beta, k);
    TPZFNMatrix<4,STATE> jac(3,1);
    TPZManVector<STATE> fxnvec(3);
    DDistFunc2new(sigmatrial, theta,beta,k,kprev,fxnvec);
    for(int k=0; k<3; k++) jac(k,0) = fxnvec[k];
//    DDistFunc2(sigmatrial, theta, beta, k, kprev, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i=1; i<=10; i++) {
        STATE difftheta = deltatheta*i/10.;
        STATE thetanext = theta+difftheta;
        STATE diffbeta = deltabeta*i/10.;
        STATE betanext = beta + diffbeta;
        STATE diffk = deltak*i/10.;
        STATE knext = k+deltak;
        STATE distnext = DistF2(sigmatrial, thetanext, betanext,knext);
        STATE distguess = dist0+jac(0,0)*difftheta+jac(1,0)*diffbeta+jac(2,0)*diffk;
        xnorm[i-1] = sqrt(difftheta*difftheta+diffbeta*diffbeta+diffk*diffk);
        errnorm[i-1] = fabs(distnext-distguess);
    }
    
}

void TPZSandlerExtended::TaylorCheckDDistF2(const TPZVec<STATE> &sigmatrial, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &xnorm,
                                            TPZVec<STATE> &errnorm) const
{
    STATE deltatheta = 0.4;
    STATE deltabeta = 0.05;
    STATE deltak = 0.5;
    TPZFNMatrix<3,STATE> res0(3,1),resid(3,1),residguess(3,1),diff(3,1);
    TPZFNMatrix<9,STATE> jac(3,3);
    TPZManVector<STATE> fxnvec(3);
    DDistFunc2new(sigmatrial, theta,beta,k,kprev,fxnvec);
    for(int k=0; k<3; k++) res0(k,0) = fxnvec[k];
//    DDistFunc2(sigmatrial, theta, beta, k, kprev, res0);
    D2DistFunc2(sigmatrial, theta, beta, k, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i=1; i<=10; i++) {
        STATE difftheta = deltatheta*i/10.;
        STATE thetanext = theta+difftheta;
        STATE diffbeta = deltabeta*i/10.;
        STATE betanext = beta+diffbeta;
        STATE diffk = deltak*i/10.;
        STATE knext = k+diffk;
        diff(0) = difftheta;
        diff(1) = diffbeta;
        diff(2) = diffk;
        TPZManVector<STATE> fxnvec(3);
        DDistFunc2new(sigmatrial, thetanext,betanext,knext,kprev,fxnvec);
        for(int k=0; k<3; k++) resid(k,0) = fxnvec[k];
//        DDistFunc2(sigmatrial, thetanext, betanext, knext, kprev, resid);
        jac.Multiply(diff, residguess);
        residguess += res0;
        xnorm[i-1] = Norm(diff);
        errnorm[i-1] = Norm(resid-residguess);
    }
    
}

/// teste da derivada D(ResF2)/D(sigtrial)

void TPZSandlerExtended::TaylorCheckDDistF2DSigtrial(const TPZVec<STATE> &sigmatrial, STATE theta, STATE beta, STATE k, STATE kprev, TPZVec<STATE> &xnorm,
                                                     TPZVec<STATE> &errnorm) const
{
    TPZManVector<STATE,3> deltasigma(3,0.3);
    deltasigma[1] *= -1.;
    deltasigma[2] *= 0.7;
    
    TPZFNMatrix<3,STATE> res0(3,1),resid(3,1),residguess(3,1),diff(3,1);
    TPZFNMatrix<9,STATE> jac(3,3);
    TPZManVector<STATE> fxnvec(3);
    DDistFunc2new(sigmatrial, theta,beta,k,kprev,fxnvec);
    for(int k=0; k<3; k++) res0(k,0) = fxnvec[k];
//    DDistFunc2(sigmatrial, theta, beta, k, kprev, res0);
    GradF2SigmaTrial(sigmatrial, theta, beta, k, kprev, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i=1; i<=10; i++) {
        for(int j=0; j<3; j++) diff(j) = deltasigma[j]*i/10.;
        TPZManVector<STATE,3> sigmanext(3);
        for(int j=0; j<3; j++) sigmanext[j] = sigmatrial[j]+diff(j);
        TPZManVector<STATE> fxnvec(3);
        DDistFunc2new(sigmanext, theta,beta,k,kprev,fxnvec);
        for(int k=0; k<3; k++) resid(k,0) = fxnvec[k];
//        DDistFunc2(sigmanext, theta, beta,k,kprev,resid);
        jac.Multiply(diff, residguess);
        residguess += res0;
        xnorm[i-1] = Norm(diff);
        errnorm[i-1] = Norm(resid-residguess);
    }
    
}

#include "pzvec_extras.h"

void TPZSandlerExtended::CheckCoordinateTransformation(TPZVec<STATE> &cart)
{
    TPZManVector<STATE,3> HWCart(3),HWCyl(3),HWCart2(3),Cart2(3);
    FromPrincipalToHWCart(cart, HWCart);
    FromPrincipalToHWCyl(cart, HWCyl);
    FromHWCylToHWCart(HWCyl, HWCart2);
    FromHWCylToPrincipal(HWCyl, Cart2);
    REAL dist1 = dist(cart,Cart2);
    REAL dist2 = dist(HWCart,HWCart2);
    cout << __FUNCTION__ <<  " dist1 = " << dist1 << " dist2 = " << dist2 << endl;
    
}

/// Compute the derivative of the stress (principal s;tresses) as a function of xi and beta
void TPZSandlerExtended::DF1Cart(STATE xi, STATE beta, TPZFMatrix<STATE> &DF1) const
{
    STATE s3beta = sin(3.*beta);
    STATE c3beta = cos(3.*beta);
    STATE Gamma = (1+s3beta)+(1.-s3beta)/fPsi;
    STATE DGamma = 3.*c3beta*(1.-1./fPsi);
    STATE SQR3 = sqrt(3.);
    STATE FFI= fA - fPhi*SQR3*xi - fC*exp(fB*SQR3*xi);
    DF1.Resize(3,2);
    DF1(0,0) = 1/SQR3 + (4*(-(SQR3*fPhi) - SQR3*fB*(fA - FFI - SQR3*fPhi*xi))*Cos(beta))/
    (SQR3*Gamma);
    DF1(1,0) = 1/SQR3 + (4*(-(SQR3*fPhi) -
                                         SQR3*fB*(fA - FFI - SQR3*fPhi*xi))*Sin(beta - M_PI/6.))/(SQR3*Gamma);
    DF1(2,0)=
         1/SQR3 - (4*(-(SQR3*fPhi) - SQR3*fB*(fA - FFI - SQR3*fPhi*xi))*
                      Sin(beta + M_PI/6.))/(SQR3*Gamma);
    
    DF1(0,1) = (-4*DGamma*(FFI - fN)*Cos(beta))/(SQR3*Gamma*Gamma) -
         (4*(FFI - fN)*Sin(beta))/(SQR3*Gamma);
    DF1(1,1) = (4*(FFI - fN)*Cos(beta - M_PI/6.))/(SQR3*Gamma) -
    (4*DGamma*(FFI - fN)*Sin(beta - M_PI/6.))/(SQR3*Gamma*Gamma);
    
    DF1(2,1) = (-4*(FFI - fN)*Cos(beta + M_PI/6.))/(SQR3*Gamma) +
    (4*DGamma*(FFI - fN)*Sin(beta + M_PI/6.))/(SQR3*Gamma*Gamma);
}

/// Compute the derivative of the stress (principal s;tresses) as a function of xi and beta
void TPZSandlerExtended::DF2Cart(STATE theta, STATE beta, STATE k, TPZFMatrix<STATE> &DF2) const
{
    STATE s3beta = sin(3.*beta);
    STATE c3beta = cos(3.*beta);
    STATE Gamma = (1+s3beta)+(1.-s3beta)/fPsi;
    STATE DGamma = 3.*c3beta*(1.-1./fPsi);
    STATE SQR3 = sqrt(3.);
    STATE FFK= fA - fPhi*k - fC*exp(fB*k);
    DF2.Resize(3, 3);
    DF2(0,0) = (4*(FFK - fN)*Cos(beta)*Cos(theta))/(SQR3*Gamma) - (FFK*fR*Sin(theta))/3.;
    DF2(1,0) = (4*(FFK - fN)*Cos(theta)*Sin(beta - M_PI/6.))/(SQR3*Gamma) - (FFK*fR*Sin(theta))/3.;
    DF2(2,0) = (4*(-FFK + fN)*Cos(theta)*Sin(beta + M_PI/6.))/(SQR3*Gamma) - (FFK*fR*Sin(theta))/3.;
    
    DF2(0,1) = (-4*(FFK - fN)*(DGamma*Cos(beta) + Gamma*Sin(beta))*Sin(theta))/(SQR3*Gamma*Gamma);
    DF2(1,1) = (4*(FFK - fN)*(Gamma*Cos(beta - M_PI/6.) - DGamma*Sin(beta - M_PI/6.))*Sin(theta))/(SQR3*Gamma*Gamma);
    DF2(2,1) = (-4*(FFK - fN)*(Gamma*Cos(beta + M_PI/6.) - DGamma*Sin(beta + M_PI/6.))*Sin(theta))/(SQR3*Gamma*Gamma);
    
    DF2(0,2) = (Gamma - (fA*fB + fPhi - fB*(FFK + fPhi*k))*
                (fR*Gamma*Cos(theta) + 4*SQR3*Cos(beta)*Sin(theta)))/(3.*Gamma);
    DF2(1,2) = (Gamma - (fA*fB + fPhi - fB*(FFK + fPhi*k))*
                (fR*Gamma*Cos(theta) + 4*SQR3*Sin(beta - M_PI/6.)*Sin(theta)))/(3.*Gamma);
    DF2(2,2) = (Gamma - (fA*fB + fPhi - fB*(FFK + fPhi*k))*
                (fR*Gamma*Cos(theta) - 4*SQR3*Sin(beta + M_PI/6.)*Sin(theta)))/(3.*Gamma);
}


void TPZSandlerExtended::TaylorCheckDF1Cart(STATE xi, STATE beta,TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const
{
    STATE deltaxi = 0.4;
    STATE deltabeta = 0.05;
    TPZFNMatrix<3,STATE> res0(3,1),resid(3,1),residguess(3,1),diff(2,1);
    TPZFNMatrix<9,STATE> jac(3,2);
    TPZManVector<STATE> sigHWCyl(3),sigCart(3);
    F1Cyl(xi, beta, sigHWCyl);
    FromHWCylToPrincipal(sigHWCyl, sigCart);
    for (int i=0; i<3; i++) {
        res0(i) = sigCart[i];
    }
    DF1Cart(xi, beta, jac);
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i=1; i<=10; i++) {
        STATE diffxi = deltaxi*i/10.;
        STATE xinext = xi+diffxi;
        STATE diffbeta = deltabeta*i/10.;
        diff(0) = diffxi;
        diff(1) = diffbeta;
        STATE betanext = beta + diffbeta;
        F1Cyl(xinext, betanext, sigHWCyl);
        FromHWCylToPrincipal(sigHWCyl, sigCart);
        for (int ii=0; ii<3; ii++) {
            resid(ii) = sigCart[ii];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        xnorm[i-1] = Norm(diff);
        errnorm[i-1] = Norm(resid-residguess);
    }
}

void TPZSandlerExtended::TaylorCheckDF2Cart(STATE theta, STATE beta, STATE k, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const
{
    STATE deltatheta = 0.4;
    STATE deltabeta = 0.05;
    STATE deltak = 0.5;
    TPZFNMatrix<3,STATE> res0(3,1),resid(3,1),residguess(3,1),diff(3,1);
    TPZFNMatrix<9,STATE> jac(3,3);
    DF2Cart(theta, beta, k, jac);
    TPZManVector<STATE> sigHWCyl(3),sigCart(3);
    F2Cyl(theta, beta, k, sigHWCyl);
    FromHWCylToPrincipal(sigHWCyl, sigCart);
    for (int i=0; i<3; i++) {
        res0(i) = sigCart[i];
    }
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i=1; i<=10; i++) {
        STATE difftheta = deltatheta*i/10.;
        STATE thetanext = theta+difftheta;
        STATE diffbeta = deltabeta*i/10.;
        STATE betanext = beta+diffbeta;
        STATE diffk = deltak*i/10.;
        STATE knext = k+diffk;
        diff(0) = difftheta;
        diff(1) = diffbeta;
        diff(2) = diffk;
        F2Cyl(thetanext, betanext, knext, sigHWCyl);
        FromHWCylToPrincipal(sigHWCyl, sigCart);
        for (int i=0; i<3; i++) {
            resid(i) = sigCart[i];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        xnorm[i-1] = Norm(diff);
        errnorm[i-1] = Norm(resid-residguess);
    }
    
}

void TPZSandlerExtended::TaylorCheckProjectSigma(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const
{
    TPZManVector<STATE,3> deltasigma(3,-0.000012), sigproj(3);
    deltasigma[1] = -4.e-6;
    deltasigma[2] = -4.e-6;
//    TPZManVector<STATE,3> deltasigma(3,-0.01), sigproj(3);
//    deltasigma[1] *= 3;
//    deltasigma[2] *= -2;

    TPZFNMatrix<3,STATE> res0(3,1),diff(3,1),resid(3,1),residguess(3,1);
    TPZFNMatrix<9,STATE> jac(3,3);
    STATE kproj;
    ProjectSigmaDep(sigtrial, kprev, sigproj, kproj, jac);
    for(int j=0; j<3; j++) res0(j) = sigproj[j];
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i=1; i<=10; i++) {
        TPZManVector<STATE,3> diffsigma(3),nextsigma(3);
        for (int j=0; j<3; j++) {
            diffsigma[j] = deltasigma[j]*i/10.;
            nextsigma[j] = sigtrial[j]+diffsigma[j];
            diff(j) = diffsigma[j];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        ProjectSigma(nextsigma, kprev, sigproj, kproj);
        for(int j=0; j<3; j++) resid(j) = sigproj[j];
        xnorm[i-1] = Norm(diff);
        errnorm[i-1] = Norm(resid-residguess);
    }
}

/// verifies the validity of dxi/dsigtrial and dbeta/dsigtrial
void TPZSandlerExtended::TaylorCheckParamF1Sigtrial(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const
{
    TPZManVector<STATE,3> deltasigma(3,-0.01), sigproj(3);
    deltasigma[1] *= 3.;
    deltasigma[2] *= -2.;
    TPZFNMatrix<3,STATE> res0(2,1),diff(3,1),resid(2,1),residguess(2,1);
    TPZFNMatrix<9,STATE> jac(2,3), jacF1(2,2), gradF1(2,3);
    STATE kproj;
    ProjectF1(sigtrial, kprev, sigproj, kproj);
    STATE xi,beta;
    SurfaceParamF1(sigproj, xi, beta);
    res0(0) = xi;
    res0(1) = beta;
    D2DistFunc1(sigtrial, xi, beta, jacF1);
    GradF1SigmaTrial(sigtrial, xi, beta, gradF1);
    //gradF1.Transpose();
    jacF1.Solve_LU(&gradF1);
    jac = -1.*gradF1;
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i=1; i<=10; i++) {
        TPZManVector<STATE,3> diffsigma(3),nextsigma(3);
        for (int j=0; j<3; j++) {
            diffsigma[j] = deltasigma[j]*i/10.;
            nextsigma[j] = sigtrial[j]+diffsigma[j];
            diff(j) = diffsigma[j];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        ProjectF1(nextsigma, kprev, sigproj, kproj);
        SurfaceParamF1(sigproj, xi, beta);
        resid(0) = xi;
        resid(1) = beta;
        xnorm[i-1] = Norm(diff);
        errnorm[i-1] = Norm(resid-residguess);
    }
    
}


void TPZSandlerExtended::TaylorCheckProjectF1(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const
{
    TPZManVector<STATE,3> deltasigma(3,-0.01), sigproj(3);
    deltasigma[1] *= 3.;
    deltasigma[2] *= -2.;
    TPZFNMatrix<3,STATE> res0(3,1),diff(3,1),resid(3,1),residguess(3,1);
    TPZFNMatrix<9,STATE> jac(3,3), jacF1(2,2), gradF1(2,3), DF1cart(2,3), GradSigma(3,3);
    STATE kproj;
    ProjectF1(sigtrial, kprev, sigproj, kproj);
    STATE xi,beta;
    SurfaceParamF1(sigproj, xi, beta);
    res0(0) = sigproj[0];;
    res0(1) = sigproj[1];
    res0(2) = sigproj[2];
    D2DistFunc1(sigtrial, xi, beta, jacF1);
    GradF1SigmaTrial(sigtrial, xi, beta, gradF1);
    //gradF1.Transpose();
    jacF1.Solve_LU(&gradF1);
    DF1Cart(xi, beta, DF1cart);
    DF1cart.Multiply(gradF1, GradSigma);
    GradSigma *= -1.;

    jac = GradSigma;
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i=1; i<=10; i++) {
        TPZManVector<STATE,3> diffsigma(3),nextsigma(3);
        for (int j=0; j<3; j++) {
            diffsigma[j] = deltasigma[j]*i/10.;
            nextsigma[j] = sigtrial[j]+diffsigma[j];
            diff(j) = diffsigma[j];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        ProjectF1(nextsigma, kprev, sigproj, kproj);
        resid(0) = sigproj[0];
        resid(1) = sigproj[1];
        resid(2) = sigproj[2];
        xnorm[i-1] = Norm(diff);
        errnorm[i-1] = Norm(resid-residguess);
    }
    
}

void TPZSandlerExtended::TaylorCheckProjectF2(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const
{
//    TPZManVector<STATE,3> deltasigma(3,-0.01), sigproj(3);
//    deltasigma[1] *= 3.;
//    deltasigma[2] *= -2.;
    TPZManVector<STATE,3> deltasigma(3,-0.000012), sigproj(3);
    deltasigma[1] = -4.e-6;
    deltasigma[2] = -4.e-6;
    TPZFNMatrix<3,STATE> res0(3,1),diff(3,1),resid(3,1),residguess(3,1);
    TPZFNMatrix<9,STATE> jac(3,3), jacF2(3,3), gradF2(3,3), DF2cart(3,3), GradSigma(3,3);
    STATE kproj;
    ProjectF2(sigtrial, kprev, sigproj, kproj);
    STATE theta,beta;
    SurfaceParamF2(sigproj, kproj, theta, beta);
    res0(0) = sigproj[0];;
    res0(1) = sigproj[1];
    res0(2) = sigproj[2];
    D2DistFunc2(sigtrial, theta, beta, kproj, jacF2);
    
    TFad<3,STATE> thetafad(theta,0),betafad(beta,1),kprojfad(kproj,2);
    TPZManVector<TFad<3,STATE>, 3> sigtrialfad(3), ddistf2(3);
    for (int m=0; m<3; m++) {
        sigtrialfad[m] = sigtrial[m];
    }
    DDistFunc2(sigtrialfad, thetafad, betafad, kprojfad, kprev, ddistf2);
    
    TPZFNMatrix<9,STATE> diffjac(3,3);
    for (int m=0; m<3; m++) {
        for (int n=0; n<3; n++) {
            diffjac(m,n) = jacF2(m,n)-ddistf2[m].fastAccessDx(n);
            jacF2(m,n) = ddistf2[m].fastAccessDx(n);
        }
    }

    
    GradF2SigmaTrial(sigtrial, theta, beta, kproj, kprev, gradF2);
    //gradF1.Transpose();
    jacF2.Solve_LU(&gradF2);
    DF2Cart(theta, beta, kproj, DF2cart);
    DF2cart.Multiply(gradF2, GradSigma);
    GradSigma *= -1.;
    
    jac = GradSigma;
    xnorm.resize(10);
    errnorm.resize(10);
    for (int i=1; i<=10; i++) {
        TPZManVector<STATE,3> diffsigma(3),nextsigma(3);
        for (int j=0; j<3; j++) {
            diffsigma[j] = deltasigma[j]*i/10.;
            nextsigma[j] = sigtrial[j]+diffsigma[j];
            diff(j) = diffsigma[j];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
        ProjectF2(nextsigma, kprev, sigproj, kproj);
        resid(0) = sigproj[0];
        resid(1) = sigproj[1];
        resid(2) = sigproj[2];
        xnorm[i-1] = Norm(diff);
        errnorm[i-1] = Norm(resid-residguess);
    }
    
}

void TPZSandlerExtended::TaylorCheckDtbkDsigtrial(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const
{
    //    TPZManVector<STATE,3> deltasigma(3,-0.01), sigproj(3);
    //    deltasigma[1] *= 3.;
    //    deltasigma[2] *= -2.;
    TPZManVector<STATE,3> deltasigma(3,-0.000012), sigproj(3);
    deltasigma[1] = -4.e-6;
    deltasigma[2] = -4.e-6;
    TPZFNMatrix<3,STATE> res0(3,1),diff(3,1),resid(3,1),residguess(3,1);
    TPZFNMatrix<9,STATE> jac(3,3), jacF2(3,3), gradF2(3,3), GradTBK(3,3);
    STATE kproj;
    ProjectF2(sigtrial, kprev, sigproj, kproj);
    STATE theta,beta;
    SurfaceParamF2(sigproj, kproj, theta, beta);
    res0(0) = theta;
    res0(1) = beta;
    res0(2) = kproj;
    D2DistFunc2(sigtrial, theta, beta, kproj, jacF2);
    TFad<3,STATE> thetafad(theta,0),betafad(beta,1),kprojfad(kproj,2);
    TPZManVector<TFad<3,STATE>, 3> sigtrialfad(3), ddistf2(3);
    for (int m=0; m<3; m++) {
        sigtrialfad[m] = sigtrial[m];
    }
    DDistFunc2(sigtrialfad, thetafad, betafad, kprojfad, kprev, ddistf2);
    
    TPZFNMatrix<9,STATE> diffjac(3,3);
    for (int m=0; m<3; m++) {
        for (int n=0; n<3; n++) {
            diffjac(m,n) = jacF2(m,n)-ddistf2[m].fastAccessDx(n);
            jacF2(m,n) = ddistf2[m].fastAccessDx(n);
        }
    }
    diffjac.Print("DiffMatrix");
    
    GradF2SigmaTrial(sigtrial, theta, beta, kproj, kprev, gradF2);
    //gradF1.Transpose();
    TPZFMatrix<STATE> jacF2temp(jacF2),gradF2temp(gradF2);
    jacF2temp.Solve_LU(&gradF2temp);
    GradTBK = gradF2temp;
    GradTBK *= -1.;
    
    jac = GradTBK;
    
        TPZFNMatrix<3,STATE> tbk(3,1,0.), rhs;
        tbk(1,0) = 0.1;
        rhs = tbk;
        GradTBK.Solve_LU(&rhs);
        for (int i=0; i<3; i++) {
            deltasigma[i] = rhs(i,0);
        }
    
    
    xnorm.resize(10);
    errnorm.resize(10);
    TPZFNMatrix<30,STATE> erros(3,10);
    for (int i=1; i<=10; i++) {
        TPZManVector<STATE,3> diffsigma(3),nextsigma(3);
        TPZFNMatrix<3,STATE> difftbk(3,1);
        for (int j=0; j<3; j++) {
            diffsigma[j] = deltasigma[j]*i/10.;
            difftbk(j) = tbk[j]*i/10.;
            nextsigma[j] = sigtrial[j]+diffsigma[j];
            diff(j) = diffsigma[j];
        }
        jac.Multiply(diff, residguess);
        residguess += res0;
//        TPZFNMatrix<3,STATE> multsig(3,0),multbk(3,0);
//        jacF2.Multiply(difftbk, multbk);
//        gradF2.Multiply(diff, multsig);
//        residguess = multbk;
//        residguess = multsig;
//        TPZFNMatrix<3,STATE> fxn(3,1);
//        DDistFunc2(nextsigma, theta+difftbk[0],beta+difftbk[1],kproj+difftbk[2],kprev,fxn);
//        TPZManVector<STATE> fxnvec(3);
//        DDistFunc2<STATE>(sigtrial,theta+difftbk[0],beta+difftbk[1],kproj+difftbk[2],kprev,fxnvec);
//        for(int k=0; k<3; k++) fxn(k,0) = fxnvec[k];
//        DDistFunc2(sigtrial, theta+difftbk[0],beta+difftbk[1],kproj+difftbk[2],kprev,fxn);
//        DDistFunc2(nextsigma, theta,beta,kproj,kprev,fxn);
            ProjectF2(nextsigma, kprev, sigproj, kproj);
            STATE theta,beta;
            SurfaceParamF2(sigproj, kproj, theta, beta);
        resid(0) = theta;
        resid(1) = beta;
        resid(2) = kproj;
        xnorm[i-1] = Norm(diff);
        errnorm[i-1] = Norm(resid-residguess);
        for(int a=0; a<3; a++) erros(a,i-1) = fabs(resid[a]-residguess[a]);
    }
    erros.Print(cout);
    
}

void TPZSandlerExtended::MCormicRanchSand(TPZSandlerExtended &mat)//em ksi
{
    STATE E=100,nu=0.25,A=0.25,B=0.67,C=0.18,D=0.67,R=2.5,W=0.066,N=0.,phi=0,psi=1.0;
    STATE G=E/(2.*(1.+nu));
    STATE K=E/(3.*(1.-2*nu));
    mat.fA=A;
    mat.fB=B;
    mat.fC=C;
    mat.fD=D;
    mat.fK=K;
    mat.fG=G;
    mat.fW=W;
    mat.fR=R;
    mat.fPhi=phi;
    mat.fN=N;
    mat.fPsi=psi;
    mat.fE=E;
    mat.fnu=nu;

}

void TPZSandlerExtended::ReservoirSandstone(TPZSandlerExtended &mat)//em ksi
{
    STATE E=1305,nu=0.25,A=2.61,B=0.169,C=2.57,D=0.05069,R=1.5,W=0.0908,N=0.,phi=0,psi=1.0;
    STATE G=E/(2.*(1.+nu));
    STATE K=E/(3.*(1.-2*nu));
    mat.fA=A;
    mat.fB=B;
    mat.fC=C;
    mat.fD=D;
    mat.fK=K;
    mat.fG=G;
    mat.fW=W;
    mat.fR=R;
    mat.fPhi=phi;
    mat.fN=N;
    mat.fPsi=psi;
    mat.fE=E;
    mat.fnu=nu;
    
    
}

void TPZSandlerExtended::SalemLimestone(TPZSandlerExtended &mat)// em MPa
{
    STATE E=22547.,nu=0.2524,A=689.2,
    B=3.94e-4,C=675.2,D=1.47e-3,R=28,W=0.08,N=6.,phi=0,psi=1.0;
    STATE G=E/(2.*(1.+nu));
    STATE K=E/(3.*(1.-2*nu));
    mat.fA=A;
    mat.fB=B;
    mat.fC=C;
    mat.fD=D;
    mat.fK=K;
    mat.fG=G;
    mat.fW=W;
    mat.fR=R;
    mat.fPhi=phi;
    mat.fN=N;
    mat.fPsi=psi;
    mat.fE=E;
    mat.fnu=nu;
}

void TPZSandlerExtended::PreSMat(TPZSandlerExtended &mat)// em MPa
{
    
    STATE E=29269,nu=0.203,A=116.67,
    B=0.0036895,C=111.48,D=0.018768,R=0.91969,W=0.006605,N=0.,phi=0,psi=1.0;
    STATE G=E/(2.*(1.+nu));
    STATE K=E/(3.*(1.-2*nu));
    mat.fA=A;
    mat.fB=B;
    mat.fC=C;
    mat.fD=D;
    mat.fK=K;
    mat.fG=G;
    mat.fW=W;
    mat.fR=R;
    mat.fPhi=phi;
    mat.fN=N;
    mat.fPsi=psi;
    mat.fE=E;
    mat.fnu=nu;
}
//REAL E = 29269,
//poisson = 0.203;
//
//material.fER.SetUp(E, poisson);
//
//REAL A = 116.67,
//B = 0.0036895,
//C = 111.48,
//D = 0.018768,
//R = 0.91969,
//W = 0.006605;
