/**
 * @file
 * @brief DEPRECATED FILE. Contains the implementations of the TPZPlaca2 methods.
 */
#include "pzmaterial.h"
//#include "pztempmat.h"

#include "pzbndcond.h"
#include <math.h>
#include <fstream>
using namespace std;
#include "pzvec.h"
#include "pzerror.h"


TPZPlaca2::TPZPlaca2(int num, REAL h, REAL f, REAL E1 , REAL E2 ,
		     REAL ni1 , REAL ni2 , REAL G12 , REAL G13 ,
		     REAL G23 , TPZFMatrix &naxes, TPZVec<REAL> &xf) : 
  TPZMaterial(num), 
  fnaxes(naxes),
  fE1(E1), fE2(E2), fG12(G12), fG13(G13), fG23(G23),
  fh(h),ff(f),fmi(1./(-1.+ni1*ni2)),fni1(ni1),fni2(ni2),
  fRmat(6,6,0.),fRmatT(6,6,0.),
  fKxxR(6,6,0.),fKyyR(6,6,0.),fKxyR(6,6,0.),fKyxR(6,6,0.),fBx0R(6,6,0.),fXF(xf)
{
  
  
   TPZFMatrix Kxx(6,6,0.),Kxy(6,6,0.),Kyx(6,6,0.),Kyy(6,6,0.),
              Bx0(6,6,0.),B0x(6,6,0.),By0(6,6,0.),B0y(6,6,0.),B00(6,6,0.),
              B0xR(6,6,0.),By0R(6,6,0.),B0yR(6,6,0.),B00R(6,6,0.);
   TPZFMatrix Kn1n1(6,6,0.),Kn1n2(6,6,0.),Kn2n1(6,6,0.),Kn2n2(6,6,0.),
              Bn10(6,6,0.),B0n1(6,6,0.),Bn20(6,6,0.),B0n2(6,6,0.),B000(6,6,0.);
//   TPZFMatrix fRmat(6,6),fRmatT(6,6);
   REAL Small , k, mi;
   Small = E1*REAL(1.E-5);
   k = 5./6.; // coeficiente de cisalhamento
   mi = 1.0/(-1.0 + ni1 * ni2);

   fRmat(0,0) = fnaxes(0,0); fRmat(0,1) = fnaxes(0,1); fRmat(0,2) = fnaxes(0,2);
   fRmat(1,0) = fnaxes(1,0); fRmat(1,1) = fnaxes(1,1); fRmat(1,2) = fnaxes(1,2);
   fRmat(2,0) = fnaxes(2,0); fRmat(2,1) = fnaxes(2,1); fRmat(2,2) = fnaxes(2,2);

   fRmat(3,3) = fnaxes(0,0); fRmat(3,4) = fnaxes(0,1); fRmat(3,5) = fnaxes(0,2);
   fRmat(4,3) = fnaxes(1,0); fRmat(4,4) = fnaxes(1,1); fRmat(4,5) = fnaxes(1,2);
   fRmat(5,3) = fnaxes(2,0); fRmat(5,4) = fnaxes(2,1); fRmat(5,5) = fnaxes(2,2);

   fRmat.Transpose(&fRmatT);



   Kxx(0,0) = -E1*h*mi;
   Kxx(0,4) = -E1*f*h*mi;
   Kxx(1,1) =  G12*h/2. + h*Small/4.;
   Kxx(1,3) = -f*G12*h/2.;
   Kxx(2,2) =  G13*h*k/2.;
   Kxx(3,1) = -f*G12*h/2.;
   Kxx(3,3) =  f*f*G12*h/2. + G12*h*h*h/24.;
   Kxx(4,0) = -E1*f*h*mi;
   Kxx(4,4) = -E1*f*f*h*mi - E1*h*h*h*mi/12.;

   Kxy(0,1) = -E1*h*mi*ni2;
   Kxy(0,3) =  E1*f*h*mi*ni2;
   Kxy(1,0) =  G12*h/2. - h*Small/4.;
   Kxy(1,4) =  f*G12*h/2.;
   Kxy(3,0) = -f*G12*h/2.;
   Kxy(3,4) = -f*f*G12*h/2. - G12*h*h*h/24.;
   Kxy(4,1) = -E1*f*h*mi*ni2;
   Kxy(4,3) =  E1*f*f*h*mi*ni2 + E1*h*h*h*mi*ni2/12.;

   Kxy.Transpose(&Kyx);

   Kyy(0,0) =  G12*h/2. + h*Small/4.;
   Kyy(0,4) =  f*G12*h/2.;
   Kyy(1,1) = -E2*h*mi;
   Kyy(1,3) =  E2*f*h*mi;
   Kyy(2,2) =  G23*h*k/2.;
   Kyy(3,1) =  E2*f*h*mi;
   Kyy(3,3) = -E2*f*f*h*mi - E2*h*h*h*mi/12.;
   Kyy(4,0) =  f*G12*h/2.;
   Kyy(4,4) =  f*f*G12*h/2. + G12*h*h*h/24.;

   B0x(4,2) =  G13*h*k/2.;
   B0x(5,1) = -h*Small/2.;

   B0x.Transpose(&Bx0);

   B0y(3,2) = -G23*h*k/2.;
   B0y(5,0) =  h*Small/2.;

   B0y.Transpose(&By0);

   B00(3,3) =  G23*h*k/2.;
   B00(4,4) =  G13*h*k/2.;
   B00(5,5) =  h*Small;

   fKxxR = fRmatT * (Kxx * fRmat);
   fKyxR = fRmatT * (Kyx * fRmat);
   fKxyR = fRmatT * (Kxy * fRmat);
   fKyyR = fRmatT * (Kyy * fRmat);
   fB0xR = fRmatT * (B0x * fRmat);
   fB0yR = fRmatT * (B0y * fRmat);
   fBx0R = fRmatT * (Bx0 * fRmat);
   fBy0R = fRmatT * (By0 * fRmat);
   fB00R = fRmatT * (B00 * fRmat);

}

static ofstream placatest("placatest.dat");
void TPZPlaca2::Contribute(TPZMaterialData &data,
                           REAL weight,
                           TPZFMatrix &ek,
                           TPZFMatrix &ef) {
  // this method adds the contribution of the material to the stiffness
  // matrix and right hand side

  // check on the validity of the arguments
  //rows x cols

TPZFMatrix &dphi = data.dphix;
// TPZFMatrix &dphiL = data.dphixl;
// TPZFMatrix &dphiR = data.dphixr;
TPZFMatrix &phi = data.phi;
// TPZFMatrix &phiL = data.phil;
// TPZFMatrix &phiR = data.phir;
// TPZManVector<REAL,3> &normal = data.normal;
TPZManVector<REAL,3> &x = data.x;
// int &POrder=data.p;
// int &LeftPOrder=data.leftp;
// int &RightPOrder=data.rightp;
// TPZVec<REAL> &sol=data.sol;
// TPZVec<REAL> &solL=data.soll;
// TPZVec<REAL> &solR=data.solr;
// TPZFMatrix &dsol=data.dsol;
// TPZFMatrix &dsolL=data.dsoll;
// TPZFMatrix &dsolR=data.dsolr;
// REAL &faceSize=data.HSize;
// TPZFMatrix &daxesdksi=data.daxesdksi;
TPZFMatrix &axes=data.axes;

  if(phi.Cols() != 1 || dphi.Rows() != 2 || phi.Rows() != dphi.Cols()){
    PZError << "TPZPlaca2.contr, inconsistent input data : phi.Cols() = "
	    << phi.Cols() << " dphi.Cols + " << dphi.Cols() <<
      " phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
      dphi.Rows() << "\n";
  }
  if(fForcingFunction) {
    fForcingFunction(x,fXF);//fXf = xfloat
  }


   REAL Dax1n1, Dax1n2, Dax2n1, Dax2n2,intern33;

   Dax1n1 = axes(0,0)* fnaxes(0,0) + axes(0,1)* fnaxes(0,1) + axes(0,2)* fnaxes(0,2);
   Dax1n2 = axes(0,0)* fnaxes(1,0) + axes(0,1)* fnaxes(1,1) + axes(0,2)* fnaxes(1,2);
   Dax2n1 = axes(1,0)* fnaxes(0,0) + axes(1,1)* fnaxes(0,1) + axes(1,2)* fnaxes(0,2);
   Dax2n2 = axes(1,0)* fnaxes(1,0) + axes(1,1)* fnaxes(1,1) + axes(1,2)* fnaxes(1,2);
   intern33 = axes(2,0)* fnaxes(2,0) + axes(2,1)* fnaxes(2,1) + axes(2,2)* fnaxes(2,2);
   if(fabs(fabs(intern33)-1.) > 1.e-6) {
	   PZError << "third axes of the material not parallel to element axis " << intern33 << endl;
	   PZError << axes(2,0) << ' ' << axes(2,1) << ' ' << axes(2,2) << endl;
	   PZError << fnaxes(2,0) << ' ' << fnaxes(2,1) << ' ' << fnaxes(2,2) << endl;
   }

   TPZFMatrix Kn1n1(6,6,0.),Kn1n2(6,6,0.),Kn2n1(6,6,0.),Kn2n2(6,6,0.),
              Bn10(6,6,0.),B0n1(6,6,0.),Bn20(6,6,0.),B0n2(6,6,0.),B000(6,6,0.);


   Kn1n1 = fKxxR * Dax1n1 * Dax1n1 + fKyyR * Dax1n2 * Dax1n2 +
           fKxyR * Dax1n1 * Dax1n2 + fKyxR * Dax1n1 * Dax1n2 ;

   Kn2n2 = fKxxR * Dax2n1 * Dax2n1 + fKyyR * Dax2n2 * Dax2n2 +
           fKxyR * Dax2n2 * Dax2n1 + fKyxR * Dax2n2 * Dax2n1;

   Kn1n2 = fKxxR * Dax1n1 * Dax2n1 + fKyyR * Dax1n2 * Dax2n2 +
           fKxyR * Dax1n1 * Dax2n2 + fKyxR * Dax1n2 * Dax2n1;

   Kn2n1 = fKxxR * Dax2n1 * Dax1n1 + fKyyR * Dax2n2 * Dax1n2 +
           fKxyR * Dax2n1 * Dax1n2 + fKyxR * Dax2n2 * Dax1n1;


   B0n1  = fB0xR * Dax1n1 + fB0yR * Dax1n2 ;

   Bn10  = fBx0R * Dax1n1 + fBy0R * Dax1n2 ;

   B0n2  = fB0xR * Dax2n1 + fB0yR * Dax2n2;

   Bn20  = fBx0R * Dax2n1 + fBy0R * Dax2n2;

   B000  = fB00R ;

   /*
Kn1n1.Print("Kn1n1 = ",placatest);
Kn1n2.Print("Kn1n2 = ",placatest);
Kn2n2.Print("Kn2n2 = ",placatest);
Kn2n1.Print("Kn2n1 = ",placatest);
B0n1.Print("B0n1 = ",placatest);
Bn10.Print("Bn10 = ",placatest);
B0n2.Print("B0n2 = ",placatest);
Bn20.Print("B0n1 = ",placatest);
B000.Print("B00 = ",placatest);
placatest.flush();
*/
  int nshape = phi.Rows();
  TPZFMatrix KIJ(6,6);



  int idf,jdf,i,j;
  for(i=0; i<nshape; i++) {
  	for(idf=0; idf<6; idf++) ef(6*i+idf) += weight*fXF[idf]*phi(i,0);
  	for(j=0; j<nshape; j++) {
     KIJ = weight*(dphi(0,i)*dphi(0,j)*Kn1n1+
                          dphi(0,i)*dphi(1,j)*Kn1n2+
                          dphi(1,i)*dphi(0,j)*Kn2n1+
                          dphi(1,i)*dphi(1,j)*Kn2n2+
                          dphi(0,i)*phi(j)   *Bn10 +
                          dphi(1,i)*phi(j)   *Bn20 +
                          phi(i)   *dphi(0,j)*B0n1 +
                          phi(i)   *dphi(1,j)*B0n2 +
                          phi(i)   *phi(j)   *B000 );
      for(idf=0; idf<6; idf++) for(jdf=0; jdf<6; jdf++)
      	ek(i*6+idf,j*6+jdf) += KIJ(idf,jdf);
   }
  }
}

void TPZPlaca2::ContributeBC(TPZMaterialData &data,
                             REAL weight,
                             TPZFMatrix &ek,
                             TPZFMatrix &ef,
                             TPZBndCond &bc) {

// TPZFMatrix &dphi = data.dphix;
// TPZFMatrix &dphiL = data.dphixl;
// TPZFMatrix &dphiR = data.dphixr;
TPZFMatrix &phi = data.phi;
// TPZFMatrix &phiL = data.phil;
// TPZFMatrix &phiR = data.phir;
// TPZManVector<REAL,3> &normal = data.normal;
// TPZManVector<REAL,3> &x = data.x;
// int &POrder=data.p;
// int &LeftPOrder=data.leftp;
// int &RightPOrder=data.rightp;
// TPZVec<REAL> &sol=data.sol;
// TPZVec<REAL> &solL=data.soll;
// TPZVec<REAL> &solR=data.solr;
// TPZFMatrix &dsol=data.dsol;
// TPZFMatrix &dsolL=data.dsoll;
// TPZFMatrix &dsolR=data.dsolr;
// REAL &faceSize=data.HSize;
// TPZFMatrix &daxesdksi=data.daxesdksi;
// TPZFMatrix &axes=data.axes;


  if(bc.Material().operator ->() != this){
    PZError << "TPZMat1dLin.apply_bc warning : this material didn't create the boundary condition!\n";
  }

  if(bc.Type() < 0 && bc.Type() > 2){
    PZError << "TPZMat1dLin.aplybc, unknown boundary condition type :"  <<
      bc.Type() << " boundary condition ignored\n";
  }

  int numdof = NStateVariables();
  int numnod = ek.Rows()/numdof;
  int r = numdof;

  int idf,jdf,in,jn;
  switch(bc.Type()){

     case 0:
         for(in=0 ; in<numnod ; ++in){
            for(idf = 0;idf<r;idf++) {
               (ef)(in*r+idf,0) += gBigNumber*phi(in,0)*bc.Val2()(idf,0)*weight;
            }
            for(jn=0 ; jn<numnod ; ++jn) {
               for(idf = 0;idf<r;idf++) {
                 ek(in*r+idf,jn*r+idf) += gBigNumber*phi(in,0)*phi(jn,0)*weight;
               }
            }
         }
         break;

     case 1:
         for(in=0 ; in<numnod ; ++in){
            for(idf = 0;idf<r;idf++) {
               (ef)(in*r+idf,0) += phi(in,0)*bc.Val2()(idf,0)*weight;
            }
         }
         break;

     case 2:
         for(in=0 ; in<numnod ; ++in){
            for(idf = 0;idf<r;idf++) {
               (ef)(in*r+idf,0) += phi(in,0)*bc.Val2()(idf,0)*weight;
            }
            for(jn=0 ; jn<numnod ; ++jn) {
               for(idf = 0;idf<r;idf++) {
                 for(jdf = 0;jdf<r;jdf++) {
                   ek(in*r+idf,jn*r+jdf) += bc.Val1()(idf,jdf)*phi(in,0)*phi(jn,0)*weight;
                 }
               }
            }
         }//fim switch
  }
}

int TPZPlaca2::NFluxes() {return 1;}

void TPZPlaca2::Flux(TPZVec<REAL> &/*x*/,TPZVec<REAL> &/*u*/,TPZFMatrix &/*dudx*/,TPZFMatrix &/*axes*/,TPZVec<REAL> &/*fl*/) {
  PZError << "TPZPlaca2::Flux is called\n";
}

void TPZPlaca2::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &/*u*/,TPZFMatrix &/*dudx*/,TPZFMatrix &/*axes*/,TPZVec<REAL> &/*flux*/,
			 TPZVec<REAL> &/*u_exact*/,TPZFMatrix &/*du_exact*/,TPZVec<REAL> &/*values*/) {
  PZError << "TPZPlaca2::Errors is called\n";
}

void TPZPlaca2::Print(ostream & out) {
  //out << "Material type TPZPlaca2 -- number = " << Id() << "\n";
  //out << "Matrix xk ->  "; fXk.Print("fXk",out);
  //out << "Matrix xC ->  "; fCk.Print("fCf",out);
  //out << "Matrix xf ->  "; fXf.Print("fXf",out);
}

/*TPZBndCond *TPZPlaca2::CreateBC(int num,int typ,TPZFMatrix &val1,TPZFMatrix &val2) {
  PZError << "TPZPlaca2::CreateBC is called\n";
  return 0;
}        */

  /**returns the variable index associated with the name*/
int TPZPlaca2::VariableIndex(const std::string &name){
	if(!strcmp(name,"Displacement6")) return 0;
	if(!strcmp(name,"Deslocx")) return 2;// Desloc. eixo x global
	if(!strcmp(name,"Deslocy")) return 3;// Desloc. eixo y global
	if(!strcmp(name,"Deslocz")) return 4;// Desloc. eixo z global
   if(!strcmp(name,"Mn1"))     return 5;// Mom. fletor eixo n1 da fibra
   if(!strcmp(name,"Mn2"))     return 6;// Mom. fletor eixo n2 da fibra
   if(!strcmp(name,"Mn1n2"))   return 7;// Mom. volvente eixos n1 e n2 da fibra
   if(!strcmp(name,"Vn1"))     return 8;// forca cortante Vn1 (positiva se antihorario)
   if(!strcmp(name,"Vn2"))     return 9;// forca cortante Vn2 (positiva se antihorario)
   if(!strcmp(name,"Sign1"))   return 10;// tens� normal na dire�o n1
   if(!strcmp(name,"Sign2"))   return 11;// tens� normal na dire�o n2
   if(!strcmp(name,"Taun1n2")) return 12;// tens� cisalhamento eixos n1 e n2
   if(!strcmp(name,"Taun1n3")) return 13;// tens� cisalhamento eixos n1 e n3
   if(!strcmp(name,"Taun2n3")) return 14;// tens� cisalhamento eixos n2 e n2
   if(!strcmp(name,"Displacement")) return 15;// tens� cisalhamento eixos n2 e n2



   return TPZMaterial::VariableIndex(name);
}

  /** returns the number of variables associated with the variable indexed by var.
      var is obtained by calling VariableIndex*/
int TPZPlaca2::NSolutionVariables(int var){
  	if(var == 2) return 1;
  	if(var == 3) return 1;
  	if(var == 4) return 1;
  	if(var == 5) return 1;
  	if(var == 6) return 1;
  	if(var == 7) return 1;
  	if(var == 8) return 1;
  	if(var == 9) return 1;
  	if(var == 10) return 1;
  	if(var == 11) return 1;
  	if(var == 12) return 1;
	if(var == 13) return 1;
	if(var == 14) return 1;
	if(var == 15) return 3;

   return TPZMaterial::NSolutionVariables(var);
}

  /**returns the solution associated with the var index based on the finite element approximation*/
void TPZPlaca2::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,
               TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){

	REAL k = 5./6.;

	if(var == 2) {
   	Solout.Resize(1);
   	Solout[0] = Sol[0];
      return;
   }
	if(var == 3) {
   	Solout.Resize(1);
   	Solout[0] = Sol[1];
      return;
   }
	if(var == 4) {
   	Solout.Resize(1);
   	Solout[0] = Sol[2];
      return;
   }
	if(var == 15) {
	  Solout.Resize(3);
	  Solout[0] = Sol[0];
	  Solout[1] = Sol[1];
	  Solout[2] = Sol[2];
	}


   if(var > 4) {
   	  Solout.Resize(1);
      TPZVec<REAL> Soln(6);
      TPZFMatrix DSolnax(2,6),DSolnn(2,6);
      int idf,jdf;
      for(idf=0; idf<6; idf++) {
         Soln[idf] = 0;
         DSolnax(0,idf) = 0;
         DSolnax(1,idf) = 0.;
         for(jdf=0; jdf<6; jdf++) {
            Soln[idf] += fRmat(idf,jdf)*Sol[jdf];
            DSolnax(0,idf) += fRmat(idf,jdf)*DSol(0,jdf);
            DSolnax(1,idf) += fRmat(idf,jdf)*DSol(1,jdf);
         }
      }
      REAL Dax1n1, Dax1n2, Dax2n1, Dax2n2;

      Dax1n1 = axes(0,0)* fnaxes(0,0) + axes(0,1)* fnaxes(0,1) + axes(0,2)* fnaxes(0,2);
      Dax1n2 = axes(0,0)* fnaxes(1,0) + axes(0,1)* fnaxes(1,1) + axes(0,2)* fnaxes(1,2);
      Dax2n1 = axes(1,0)* fnaxes(0,0) + axes(1,1)* fnaxes(0,1) + axes(1,2)* fnaxes(0,2);
      Dax2n2 = axes(1,0)* fnaxes(1,0) + axes(1,1)* fnaxes(1,1) + axes(1,2)* fnaxes(1,2);

      for(idf=0;idf<6;idf++) {
         DSolnn(0,idf) = Dax1n1*DSolnax(0,idf)+Dax2n1*DSolnax(1,idf);
         DSolnn(1,idf) = Dax1n2*DSolnax(0,idf)+Dax2n2*DSolnax(1,idf);
      }

      if (var == 5) {
          REAL Mn1;
          Mn1 = -fE1*fh*fh*fh*fmi*(-fni2*DSolnn(1,3)+DSolnn(0,4))/12.;
          Mn1 += fE1*ff*fh*fmi*(ff*fni2*DSolnn(1,3)-fni2*DSolnn(1,1)
                -ff*DSolnn(0,4)-DSolnn(0,0));
          Solout[0] = Mn1;
          return;
      }

      if (var == 6) {
          REAL Mn2;
          Mn2  =-(fE2*fh*fh*fh*fmi*(-DSolnn(1,3) + fni1*DSolnn(0,4)))/12.;
          Mn2 +=  fE2*ff*fh*fmi*(ff*DSolnn(1,3)  - DSolnn(1,1) -
                  ff*fni1*DSolnn(0,4) - fni1*DSolnn(0,0));
          Solout[0] = Mn2;
          return;
      }

      if (var == 7 ) {
          REAL Mn1n2;
          Mn1n2  = fG12*fh*fh*fh*(DSolnn(1,4) - DSolnn(0,3))/24.;
          Mn1n2 += (ff*fG12*fh*(ff*DSolnn(1,4) + DSolnn(1,0) -
                    ff*DSolnn(0,3) + DSolnn(0,1)))/2.;
          Solout[0] = Mn1n2;
          return;
      }
      if (var == 8) {
          REAL Vn1;
          Vn1  =  (fG12*fh*k*(Soln[4]+DSolnn(0,2)));
          Solout[0] = Vn1;
          return;
      }
      if (var == 9) {
         REAL Vn2;
         Vn2  = (fG12*fh*k*(-Soln[3]+ DSolnn(1,2)))/2.;
         Solout[0] = Vn2;
         return;
	  }

      if (var == 10 ) {
        REAL Sign1;
        REAL z;
        z=0.;
        Sign1 = -fE1*fmi*(fni2*(-(ff + z)*DSolnn(1,3) + DSolnn(1,1))+
                (ff + z)*DSolnn(0,4) + DSolnn(0,0));
        Solout[0] = Sign1;
        return;
	  }
      if (var == 11 ) {
        REAL Sign2;
        REAL z;
        z=0.;
        Sign2=-fE2*fmi*(-(ff + z)*DSolnn(1,3) + DSolnn(1,1) +
               fni1*((ff + z)*DSolnn(0,4) + DSolnn(0,0)));
        Solout[0] = Sign2;
        return;
	  }
      if (var == 12 ) {
       REAL Taun1n2;
       REAL z;
       z=0.;
       Taun1n2=(fG12*((ff + z)*DSolnn(1,4) + DSolnn(1,0) -
               (ff + z)*DSolnn(0,3) + DSolnn(0,1)))/2.;
       Solout[0] = Taun1n2;
       return;
	  }
      if (var == 13 ) {
       REAL Taun1n3;
       Taun1n3=fG13*(Soln[4] + DSolnn(0,2))/2.;
       Solout[0] = Taun1n3;
       return;
	  }
      if (var == 14 ) {
       REAL Taun2n3;
       Taun2n3=fG23*(-Soln[3] + DSolnn(1,2))/2.;
       Solout[0] = Taun2n3;
       return;
	  }
   }

   TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}


