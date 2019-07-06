// $Id: pzonedref.cpp,v 1.9.2.2 2007-09-25 13:06:00 longhin Exp $

#include "pzonedref.h"
#include "pzquad.h"
#include "pztrnsform.h"
#include "pzvec.h"
#include "pzshapelinear.h"

#include "pzlog.h"

using namespace std;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.adapt.pzonedref"));
#endif

//ofstream TPZOneDRef::fLogFile("onedref.txt");

int TPZOneDRef::gMaxP = 10;

TPZOneDRef::TPZOneDRef(int nstate) : fMS1S1(gMaxP+1,gMaxP+1,0.),
				     fMS1B(gMaxP+1,gMaxP+1,0.),
				     fMS2B(gMaxP+1,gMaxP+1,0.),
				     fMBB(gMaxP+1,gMaxP+1,0.),
				     fKS1S1(gMaxP+1,gMaxP+1,0.),
				     fKS1B(gMaxP+1,gMaxP+1,0.),
				     fKS2B(gMaxP+1,gMaxP+1,0.),
				     fKBB(gMaxP+1,gMaxP+1,0.),
				     fRhs1(gMaxP+1,nstate),
				     fRhs2(gMaxP+1,nstate),
				     fRhsb(gMaxP+1,nstate),
				     fM(gMaxP+1,gMaxP+1,0.),
				     fNState(nstate)  /* , fDest1(11), fDest2(11), fNodes1(11), fNodes2(11)*/ {
  fTryP = 1;
  IntegrateMatrices();
}

void TPZOneDRef::IntegrateMatrices(){
  TPZInt1d integr(20);
  TPZVec<int64_t> ids(2);
  ids[0] = 0;
  ids[1] = 1;
  TPZTransform t1(1),t2(1);
  t1.Mult()(0,0) = 0.5;
  t2.Mult()(0,0) = 0.5;
  t1.Sum()(0,0) = -0.5;
  t2.Sum()(0,0) = 0.5;
  TPZFMatrix<REAL>phis(gMaxP+1,1),dphis(1,gMaxP+1),phib1(gMaxP+1,1),dphib1(1,gMaxP+1),phib2(gMaxP+1,1),dphib2(1,gMaxP+1);
  int ip;
  REAL weight;
  TPZVec<REAL> pos(1),posbig1(1),posbig2(1);
  TPZManVector<int,1> maxvec(1,gMaxP);
  for(ip=0; ip<integr.NPoints(); ip++) {
    integr.Point(ip,pos,weight);
    pzshape::TPZShapeLinear::Shape(pos,ids,maxvec,phis,dphis);
    t1.Apply(pos,posbig1);
    t2.Apply(pos,posbig2);
    pzshape::TPZShapeLinear::Shape(posbig1,ids,maxvec,phib1,dphib1);
    pzshape::TPZShapeLinear::Shape(posbig2,ids,maxvec,phib2,dphib2);
    int i,j;
    for(i=0; i<gMaxP+1; i++) {
      for(j=0; j<gMaxP+1; j++) {
	fMS1S1(i,j) += weight*phis(i,0)*phis(j,0);
	fMS1B(i,j) += weight*phis(i,0)*phib1(j,0);
	fMS2B(i,j) += weight*phis(i,0)*phib2(j,0);
	fMBB(i,j) += 2.*weight*phis(i,0)*phis(j,0);
	fKS1S1(i,j) += weight*dphis(0,i)*dphis(0,j);
	fKS1B(i,j) += 0.5*weight*dphis(0,i)*dphib1(0,j);
	fKS2B(i,j) += 0.5*weight*dphis(0,i)*dphib2(0,j);
	fKBB(i,j) += 0.5*weight*dphis(0,i)*dphis(0,j);
      }
    }
  }
}

void TPZOneDRef::TransformU(TPZFMatrix<REAL>&U, TPZVec<int> &id, int p1, int p2) {
  int i,st;
  if(id[0] > id[1]) {
    for(st=0; st<fNState; st++) {
      for(i=1; i<p1-1;i+=2) U((i+1)*fNState+st,0) = -U((i+1)*fNState+st,0);
    }
  }
  if(id[1] > id[2]) {
    for(st=0; st<fNState; st++) {
      for(i=1; i<p2-1;i+=2) U((p1+1+i)*fNState+st,0) = - U((p1+1+i)*fNState+st,0);
    }
  }
}

void TPZOneDRef::LoadU(TPZFMatrix<REAL>&U, int p1, int p2, REAL delx) {

  //Cesar 2003-02-20 
  // Moved to BestPattern (public interface used)====
//   if(U.Rows() != (p1+p2+1)*fNState) {
//     cout << "TPZOneDRef inconsistent input data\n";
//     return;
//   }
//   fDelx = delx;
  //=================================================
  //   fp1 = p1;
  //   fp2 = p2;
  //   fpb = fp2 + fp1 - 2;
  //   if(fpb > gMaxP){
  //     // Cesar 2003-01-09++++++
  //     fp1 = (int)(((float)gMaxP / (float)fpb) * p1);
  //     fp2 = (int)(((float)gMaxP / (float)fpb) * p2);
  //     // -------------
  //     fpb = gMaxP;
  //   }
  // fU.Redim(2*fpb-1,fNState);
  //==================================================
  int i,j;
  for(i=0; i<fpb+1; i++) {
    for(j=0; j<fpb+1; j++) {
      //       fM(i,j) = 0.5*fDelx*fMS1S1(i,j)+2./fDelx*fKS1S1(i,j);
      fM(i,j) = 2./fDelx*fKS1S1(i,j);
    }
  }
  //  fM.Print("M Matrix");
  int st;
  //Cesar 2003-01-09 s�alterei o loop de p1 -> fp1 e p2 ->fp2
  for(st=0; st<fNState; st++) {
    for(i=0; i< fp1-1; i++) fU(i,st) = U((i+1)*fNState+st,0);
    fU(fpb-1,st) = U(p1*fNState+st,0);
    for(i=0; i< fp2-1; i++) fU(fpb+i,st) = U((p1+i+1)*fNState+st,0);
    fU(fpb-1,st) -= (U(st,0)+U((p1+p2)*fNState+st,0))*0.5;
  }
  fRhs1.Zero();
  fRhs2.Zero();
  fRhsb.Zero();
  TPZFMatrix<REAL>uloc(gMaxP+1,1,0.);
  for(st=0; st<fNState; st++) {
    uloc.Zero();
    uloc(1,0) = fU(fpb-1,st);
    for(i=0; i<p1-1; i++) uloc(i+2,0) = fU(i,st);
    
    //  uloc.Print("uloc for first element");
    for(j=0; j<fpb+1; j++) {
      for(i=0; i<fpb-1; i++) {
	fRhs1(i,st) += fM(i+2,j)*uloc(j,0);
      }
      fRhs1(fpb-1,st) += fM(1,j)*uloc(j,0);
    }
    for(j=0; j<p1+1; j++) {
      for(i=0; i<fpb-1; i++) {
	//       fRhsb(i,0) += (0.5*fDelx*fMS1B(j,i+2)+2./fDelx*fKS1B(j,i+2))*uloc(j,0);
	fRhsb(i,st) += 2./fDelx*fKS1B(j,i+2)*uloc(j,0);
      }
    }
    uloc.Zero();
    uloc(0,0) = fU(fpb-1,st);
    for(i=0; i<p2-1; i++) uloc(i+2,0) = fU(i+fpb,st);
    
    //  uloc.Print("uloc for second element");
    for(j=0; j<p2+1; j++) {
      for(i=0; i<fpb-1; i++) {
	fRhs2(i+1,st) += fM(i+2,j)*uloc(j,0);
      }
      fRhs2(0,st) += fM(0,j)*uloc(j,0);
    }
    for(j=0; j<p2+1; j++) {
      for(i=0; i<fpb-1; i++) {
	//       fRhsb(i,0) += (0.5*fDelx*fMS2B(j,i+2)+2./fDelx*fKS2B(j,i+2))*uloc(j,0);
	fRhsb(i,st) += 2./fDelx*fKS2B(j,i+2)*uloc(j,0);
      }
    }
  }
  
  // Compute the stiffness for computing the energy error
  BuildStiffness(fpb,fpb,fStiffU);
  
  //  Print("After loading the solution");
}

REAL TPZOneDRef::Error(int p1, int p2) {
  if (p1 == 0 || p2== 0){
    PZError << "TPZOneDRef::Error ERROR trying to set p1 = " << p1 << " p2 = " << p2 << endl;
    PZError << "The null value will be set equal 1" << endl;
    p1 = p1==0 ? 1 : p1;
    p2 = p2==0 ? 1 : p2;  
  }
  TPZFMatrix<REAL>stiff;
  BuildStiffness(p1,p2,stiff);
  TPZFMatrix<REAL>rhs(p1+p2-1,fNState,0.);
  int i,j,st;
  for(st=0; st<fNState; st++) {
    for(i=0; i<p1-1; i++) {
      rhs(i,st) = fRhs1(i,st);
    }
    for(i=0; i<p2; i++) {
      rhs(p1-1+i,st) = fRhs2(i,st);
    }
    rhs(p1-1,st) += fRhs1(fpb-1,st);
  }
  //  rhs.Print("right hand side for error computation");
  
  TPZFMatrix<REAL>u(fU);
  stiff.SolveDirect(rhs,ELDLt);
  
  //  rhs.Print("solution from error computation");
  
  for(st=0; st<fNState; st++) {
    for(i=0; i<p1-1; i++) u(i,st) -= rhs(i,st);
    for(i=0; i<p2; i++) u(i+fpb-1,st) -= rhs(i+p1-1,st);
  }
  REAL error = 0.;
  for(st=0; st<fNState; st++) {
    for(i=0; i<fpb+fpb-1; i++) {
      for(j=0; j<fpb+fpb-1; j++) {
        error += u(i,st)*fStiffU(i,j)*u(j,st);
      }
    }
  }
  return error;
}

void TPZOneDRef::BuildStiffness(int p1, int p2, TPZFMatrix<REAL>&stiff){

  stiff.Redim(p1+p2-1,p1+p2-1);
  
  int i,j;
  for(i=0; i< p1-1; i++) {
    stiff(p1-1,i) = fM(1,i+2);
    stiff(i,p1-1) = fM(i+2,1);
    for(j=0; j< p1-1; j++) {
      stiff(i,j) = fM(i+2,j+2);
    }
  }
  stiff(p1-1,p1-1) = fM(1,1);
  for(i=0; i<p2-1; i++) {
    stiff(p1-1,i+p1) = fM(0,i+2);
    stiff(i+p1,p1-1) = fM(i+2,0);
    for(j=0; j<p2-1; j++) {
      stiff(p1+i,p1+j) = fM(i+2,j+2);
    }
  }
  stiff(p1-1,p1-1) += fM(0,0);

  //  stiff.Print("After building the stiffness matrix");
}

REAL TPZOneDRef::Error(int pb) {

  TPZFMatrix<REAL>stiff(pb-1,pb-1);
  TPZFMatrix<REAL>rhsb(pb-1,fNState);
  int i,j,st;
  for(i=0; i<pb-1; i++) {
    for(st=0; st<fNState; st++) {
      rhsb(i,st) = fRhsb(i,st);
    }
    for(j=0; j<pb-1; j++) {
//       stiff(i,j) = fDelx*fMS1S1(i+2,j+2)+1./fDelx*fKS1S1(i+2,j+2);
      stiff(i,j) = 1./fDelx*fKS1S1(i+2,j+2);
    }
  }

  //  stiff.Print("stiffness matrix for large element");
  stiff.SolveDirect(rhsb,ELDLt);

  //  rhsb.Print("Solution for large element");

  TPZFMatrix<REAL>u1(pb,fNState,0.),MS1S1(pb,pb);
  for(i=0; i<pb; i++) {
    for(j=0; j<pb-1; j++) {
      for(st=0; st<fNState; st++) {
        u1(i,st) += fMS1B(i+1,j+2)*rhsb(j,st);
      }
    }
    for(j=0; j<pb; j++) {
      MS1S1(i,j) = fMS1S1(i+1,j+1);
    }
  }
  //****************** CHECK FROM HERE
  MS1S1.SolveDirect(u1,ELDLt);

  //  u1.Print("Solution for the left element");
  TPZFMatrix<REAL>u2(pb,fNState,0.),MS2S2(pb,pb);
  for(j=0; j<pb-1; j++) {
    for(st=0; st<fNState; st++) {
      u2(0,st) += fMS2B(0,j+2)*rhsb(j,st);
    }
  }
  for(i=1; i<pb; i++) {
    for(j=0; j<pb-1; j++) {
      for(st=0; st<fNState; st++) {
        u2(i,st) += fMS2B(i+1,j+2)*rhsb(j,st);
      }
    }
    for(j=1; j<pb; j++) {
      MS2S2(i,j) = fMS1S1(i+1,j+1);
    }
    MS2S2(i,0) = fMS1S1(i+1,0);
    MS2S2(0,i) = fMS1S1(0,i+1);
  }
  MS2S2(0,0) = fMS1S1(0,0);
  MS2S2.SolveDirect(u2,ELDLt);

  //  u2.Print("Solution for the right element");
  TPZFMatrix<REAL>u(fU);
  for(st=0; st<fNState; st++) {
    for(i=0; i<pb-1; i++) {
      u(i,st) -= u1(i+1,st);
      u(i+fpb,st) -= u2(i+1,st);
    }
    u(fpb-1,st) -= u1(0,st);
  }
  REAL error = 0.;
  for(st=0; st<fNState; st++) {
    for(i=0; i<fpb+fpb-1; i++) {
      for(j=0; j<fpb+fpb-1; j++) {
        error += u(i,st)*fStiffU(i,j)*u(j,st);
      }
    }
  }
  return error;
}

REAL TPZOneDRef::BestPattern(TPZFMatrix<REAL>&U, TPZVec<int> &id, int &p1, int &p2, int &hp1, int &hp2, REAL & hperror, REAL delx) {
  if(p1<2 && p2<2) {
    cout << "TPZOneDRef called with wrong parameters\n";
    return -1.;
  }
  //Cesar 2003-02-20
  //====================================================
  if(U.Rows() != (p1+p2+1)*fNState) {
    cout << "TPZOneDRef inconsistent input data\n";
    return -1.;
  }
  fDelx = delx;
  fp1 = p1;
  fp2 = p2;
  fpb = fp2 + fp1 - 2;
  if(fpb > gMaxP){
    // Cesar 2003-01-09++++++
    fp1 = (int)(((float)gMaxP / (float)fpb) * p1);
    fp2 = (int)(((float)gMaxP / (float)fpb) * p2);
    // -------------
    fpb = gMaxP;
    fTryP = 0;
  }
  fU.Redim(2*fpb-1,fNState);
  TransformU(U,id,fp1,fp2);
  LoadU(U,fp1,fp2,delx);
  TransformU(U,id,fp1,fp2);
  //=======================================================

  // numdof indicates the number of degrees of freedom permitted by the refinement pattern
  //int numdof = p1 < p2 ? p2 : p1;
  int numdof = fp1 < fp2 ? fp2 : fp1;
  numdof--;
  int pl1 = 1, pl2 = numdof+1-pl1;
  if(pl2 > gMaxP) pl2 = gMaxP;
  int bestp1,bestp2;
  REAL besterror;

  //Cesar 2004-13-03-->>>
  pl1 = pl1==0 ? 1 : pl1;
  pl2 = pl2==0 ? 1 : pl2;  
  //<<<---
  
  bestp1 = pl1;
  bestp2 = pl2;
  besterror = Error(pl1,pl2);
  int maxpl1 = numdof+1;
  //  if(maxpl1 > maxp+1) maxpl1 = gMaxP+1;
  
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "p1 = " << pl1 << " p2 = " << pl2 << " error = " << besterror << endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    //  fLogFile  << "p1 = " << pl1 << " p2 = " << pl2 << " error = " << besterror << endl;
//#endif

  int pl1try;
//   cout << "maxpl1 = " << maxpl1 << endl;
  for(pl1try=2; pl1try<maxpl1; pl1try++) {
    pl1 = pl1try;
    pl2 = numdof-pl1try+1;
    if(pl2 > gMaxP) pl2 = gMaxP;
    if(pl1 > gMaxP) pl1 = gMaxP;
    
    REAL error = Error(pl1,pl2);
#ifdef LOG4CXX
      {
          std::stringstream sout;
          sout << "p1 = " << pl1 << " p2 = " << pl2 << " error = " << error << endl;
          LOGPZ_DEBUG(logger,sout.str())
      }
#endif
    if(error < besterror) {
      bestp1 = pl1;
      bestp2 = pl2;
      besterror = error;
    }
//    fLogFile << "p1 = " << pl1 << " p2 = " << pl2 << " error = " << error << endl;
  }
  hp1 = bestp1;
  hp2 = bestp2;
  hperror = besterror;

  // Verifies the only p refinement. 
  if (fTryP){
    REAL error = Error(fp1);
#ifdef LOG4CXX
      {
          std::stringstream sout;
          sout << "pref = " << fp1 << " error = " << error << endl;
          LOGPZ_DEBUG(logger,sout.str())
      }
#endif
    if(error < besterror) {
      //bestp1 = p1;
      bestp1 = fp1;
      bestp2 = -1;
      besterror = error;
    }
    //fLogFile  << "pb = " << fp1 << " error = " << error << endl;
  }
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << " input p1 = " << fp1 << " input p2 = " << fp2 << " bestp1 = " << bestp1 << " bestp2 = " << bestp2 << " error = " << besterror << endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
  //fLogFile << " p1 = " << fp1 << " p2 = " << fp2 << " bestp1 = " << bestp1 << " bestp2 = " << bestp2 << " error = " << besterror << endl;
  p1 = bestp1;
  p2 = bestp2;
  return besterror;
}

void TPZOneDRef::BuildIndexVectors(int p1, int p2, int pb) {
  //???
}

void TPZOneDRef::Print(char *msg, ostream &out) {
  out << "TPZOneDRef output for object \"";
  if(msg) out << msg;
  out << "\"" << endl;
  fU.Print("Solution U",out);
  fRhs1.Print("Rhs1",out);
  fRhs2.Print("Rhs2",out);
  fRhsb.Print("Rhsb",out);
  fStiffU.Print("fStiffU",out);
}

int TPZOneDRef::main() {
  int nstate = 2;
  cout << "input number of state variables ";
  cin >> nstate;
  TPZOneDRef f(nstate);

  int poly = 10;
  int st;
  TPZFMatrix<REAL>solcubic(poly+1,nstate,1.);
  for(st=0; st<nstate; st++) {
    solcubic(0,st) = 1.;
    solcubic(1,st) = 0.;
  }
  TPZFMatrix<REAL>sol1cubic(poly+1,nstate,0.),sol2cubic(poly+1,nstate,0.);
  TPZFMatrix<REAL>M1(poly+1,poly+1);
//   f.fKS1B.Print("fKS1B");
//   f.fKS2B.Print("fKS2B");
//   f.fMS1B.Print("fMS1B");
//   f.fMS2B.Print("fMS2B");

  int i,j;
  for(i=0; i<poly+1; i++) {
    for(j=0; j<poly+1; j++) {
      M1(i,j) = f.fMS1S1(i,j)+f.fKS1S1(i,j);
      for(st=0; st<nstate; st++) {
	sol1cubic(i,st) += (f.fMS1B(i,j)+f.fKS1B(i,j))*solcubic(j,st);
	sol2cubic(i,st) += (f.fMS2B(i,j)+f.fKS2B(i,j))*solcubic(j,st);
      }
    }
  }

//   sol1cubic.Print("right hand side left element");
//   sol2cubic.Print("right hand side right element");

//   M1.Print("stiffness matrix of one element projection");
  M1.SolveDirect(sol1cubic,ELDLt);
  M1.SolveDirect(sol2cubic,ELDLt);

  sol1cubic.Print("left solution");
  sol2cubic.Print("right solution");
  
  TPZFMatrix<REAL>sol((2*poly+1)*nstate,1,0.);

  for(st=0; st<nstate; st++) {
    sol(st,0) = sol1cubic(0,st);
    for(i=0; i< poly-1; i++) sol((i+1)*nstate+st,0) = sol1cubic(i+2,st);
    sol(poly*nstate+st,0) = sol2cubic(0,st);
    for(i=0; i< poly-1; i++) sol((i+poly+1)*nstate+st,0) = sol2cubic(i+2,st);
    sol((2*poly)*nstate+st,0) = sol2cubic(1,st);
  }

  TPZVec<int> ids(3);
  for(i=0; i<3; i++) ids[i] = i;
  int p1 = poly, p2=poly;

  for(i=0; i<4; i++) {
    p1 = poly;
    p2=poly;
    int h1,h2;
    REAL herror;
    REAL error = f.BestPattern(sol,ids,p1,p2,h1,h2,herror,2.0);
    cout << "The best pattern found is p1 = " << p1 << " p2 = " << p2 << " error " << error << " h1 = " << h1 << " h2 = " << h2 << " herror " << herror << endl;
    for(st=0; st<nstate; st++) {
      for(j=0; j<poly; j++) {
	sol((poly+j)*nstate+st,0) += 0.1;
      }
    }
  }
  return 0;
}
