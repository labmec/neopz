//$Id: pzelasAXImat.cpp,v 1.5 2009-01-19 23:43:39 erick Exp $
// -*- c++ -*-
#include "pzelasAXImat.h" 
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>
#include "pzaxestools.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.axisymetric"));
#endif

#include <fstream>
using namespace std;

TPZElasticityAxiMaterial::TPZElasticityAxiMaterial() : TPZMaterial(0), f_AxisR(3,0.), f_AxisZ(3,0.),f_Origin(3,0.) {
  f_AxisZ[1] = 1.;
  f_AxisR[0] = 1.;
  fE	= -1.;  // Young modulus
  fnu	= -1.;   // poisson coefficient
  ff[0]	= 0.; // X component of the body force
  ff[1]	= 0.; // Y component of the body force
  ff[2] = 0.; // Z component of the body force - not used for this class
  fEover1MinNu2 = -1.;  //G = E/2(1-nu);
  fEover21PlusNu = -1.;//E/(1-nu)

  f_c = 0.;
  f_phi = 0.;
}

TPZElasticityAxiMaterial::TPZElasticityAxiMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy) : TPZMaterial(num) {

  fE	= E;  // Young modulus
  fnu	= nu;   // poisson coefficient
  ff[0]	= fx; // X component of the body force
  ff[1]	= fy; // Y component of the body force
  ff[2] = 0.; // Z component of the body force - not used for this class
  fEover1MinNu2 = E/(1-fnu*fnu);  //G = E/2(1-nu);
  fEover21PlusNu = E/(2.*(1+fnu));//E/(1-nu)
  f_c = 0.;
  f_phi = 0.;
}

void TPZElasticityAxiMaterial::SetOrigin(vector<REAL> &Orig, vector<REAL> &AxisZ, vector<REAL> &AxisR)
{
  if(Orig.size() == 3 && AxisZ.size() == 3 && AxisR.size() == 3)
  {
    TPZFNMatrix<6> Vecs(3,2,0.), VecsOrt(3,2,0.), JacVecsOrth(2,2,0.);
    for(int i = 0; i < 3; i++)
    {
      Vecs.PutVal(i,0,AxisR[i]);
      Vecs.PutVal(i,1,AxisZ[i]);
    }
    Vecs.GramSchmidt(VecsOrt,JacVecsOrth);

    for(int j = 0; j < 3; j++)
    {
      AxisR[j] = VecsOrt.GetVal(j,0);
      AxisZ[j] = VecsOrt.GetVal(j,1);
    }
    f_Origin = Orig;
    f_AxisZ = AxisZ;
    f_AxisR = AxisR;
  }
  else
  {
    cout << "Invalid Origin and/or Axis vector on TPZElasticityAxiMaterial()!\n";
    DebugStop();
  }
}

REAL TPZElasticityAxiMaterial::ComputeR(TPZVec<REAL> &x)
{
  return (x[0] - f_Origin[0])*f_AxisR[0] + (x[1] - f_Origin[1])*f_AxisR[1] + (x[2] - f_Origin[2])*f_AxisR[2];
}

vector<REAL> TPZElasticityAxiMaterial::GetAxisR()
{
  return f_AxisR;
}

vector<REAL> TPZElasticityAxiMaterial::GetAxisZ()
{
  return f_AxisZ;
}

vector<REAL> TPZElasticityAxiMaterial::GetOrigin()
{
  return f_Origin;
}

TPZElasticityAxiMaterial::~TPZElasticityAxiMaterial() {
}

int TPZElasticityAxiMaterial::NStateVariables() {
  return 2;
}

void TPZElasticityAxiMaterial::Print(std::ostream &out) {
  out << "name of material : " << Name() << "\n";
  out << "properties : \n";
  out << "\tE   = " << fE   << endl;
  out << "\tnu   = " << fnu   << endl;
  out << "\tF   = " << ff[0] << ' ' << ff[1]   << endl;
}

void TPZElasticityAxiMaterial::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix &ek,TPZFMatrix &ef)
{
  TPZFMatrix &dphi = data.dphix;
  TPZFMatrix &phi  = data.phi;
  TPZFMatrix &axes = data.axes;

  int phc,phr,dphc,dphr,efr,efc,ekr,ekc;
  phc = phi.Cols();
  phr = phi.Rows();
  dphc = dphi.Cols();
  dphr = dphi.Rows();
  efr = ef.Rows();
  efc = ef.Cols();
  ekr = ek.Rows();
  ekc = ek.Cols();
  if(phc != 1 || dphr != 2 || phr != dphc || ekr != phr*2 || ekc != phr*2 || efr != phr*2 || efc != 1)
  {
    PZError << "\nTPZElasticityMaterial.contr, inconsistent input data : \n" <<
      "phi.Cols() = " << phi.Cols() << " dphi.Cols() = " << dphi.Cols() <<
      " phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
      dphi.Rows() << "\nek.Rows() = " << ek.Rows() << " ek.Cols() = "
	    << ek.Cols() <<
      "\nef.Rows() = " << ef.Rows() << " ef.Cols() = "
	    << ef.Cols() << "\n";
    return;
  }
  if(fForcingFunction)
  {
    TPZManVector<REAL> res(3);
    fForcingFunction(data.x,res);
    ff[0] = res[0];
    ff[1] = res[1];
    ff[2] = res[2];
  }

  /// R = Dot[{data.x - origin},{AxisR}]   ***because AxisR is already normalized!
  REAL R = (data.x[0] - f_Origin[0])*f_AxisR[0] + (data.x[1] - f_Origin[1])*f_AxisR[1] + (data.x[2] - f_Origin[2])*f_AxisR[2];

  int s = (R > 0)? 1:-1;
  R = fabs(R);

  /**
   * Plain strain materials values
   */
  REAL nu1 = 1 - fnu;//(1-nu)
  REAL nu2 = (1-2*fnu)/2;
  REAL F = fE/((1+fnu)*(1-2*fnu));

  TPZFNMatrix<4> dphiRZi(2,1), dphiRZj(2,1);

  double axis0DOTr = 0., axis1DOTr = 0., axis0DOTz = 0., axis1DOTz = 0.;
  for(int pos = 0; pos < 3; pos++)
  {
    axis0DOTr += axes.GetVal(0,pos) * f_AxisR[pos] * s;
    axis1DOTr += axes.GetVal(1,pos) * f_AxisR[pos] * s;
    axis0DOTz += axes.GetVal(0,pos) * f_AxisZ[pos];
    axis1DOTz += axes.GetVal(1,pos) * f_AxisZ[pos];
  }

  double R2PI = 2. * M_PI * R;
  for( int in = 0; in < phr; in++ )
  {
    ef(2*in, 0)   += weight * R2PI * (ff[0] * phi(in,0)); // direcao x
    ef(2*in+1, 0) += weight * R2PI * (ff[1] * phi(in,0)); // direcao y

    //dphi_i/dr = dphi_i/axis0 <axes0,f_AxisR> + dphi_i/axis1 <axes1,f_AxisR>
    dphiRZi.PutVal(0,0, dphi.GetVal(0,in)*axis0DOTr + dphi.GetVal(1,in)*axis1DOTr );

    //dphi_i/dz = dphi_i/axis0 <axes0,f_AxisZ> + dphi_i/axis1 <axes1,f_AxisZ>
    dphiRZi.PutVal(1,0, dphi.GetVal(0,in)*axis0DOTz + dphi.GetVal(1,in)*axis1DOTz );

    for( int jn = 0; jn < phr; jn++ )
    {

      //dphi_j/dr = dphi_j/axis0 <axes0,f_AxisR> + dphi_j/axis1 <axes1,f_AxisR>
      dphiRZj.PutVal(0,0, dphi.GetVal(0,jn)*axis0DOTr + dphi.GetVal(1,jn)*axis1DOTr );

      //dphi_j/dz = dphi_j/axis0 <axes0,f_AxisZ> + dphi_j/axis1 <axes1,f_AxisZ>
      dphiRZj.PutVal(1,0, dphi.GetVal(0,jn)*axis0DOTz + dphi.GetVal(1,jn)*axis1DOTz );

      double lambda = -((fE*fnu)/((1. + fnu)*(2.*fnu-1.)));
      double mi =  fE/(2.*(1. + fnu));

      double term00 = dphiRZi(0,0) * (lambda + 2.*mi) * dphiRZj(0,0) +
                      dphiRZi(1,0) * mi * dphiRZj(1,0) +
                      phi(in,0) * lambda/R * dphiRZj(0,0) +
                      dphiRZi(0,0) * lambda/R * phi(jn,0) +
                      phi(in,0) * (lambda+2.*mi)/(R*R) * phi(jn,0);

      double term01 = dphiRZi(1,0) * mi * dphiRZj(0,0) +
                      dphiRZi(0,0) * lambda * dphiRZj(1,0) +
                      phi(in,0) * lambda/R * dphiRZj(1,0);

      double term10 = dphiRZi(1,0) * lambda * dphiRZj(0,0) +
                      dphiRZi(0,0) * mi * dphiRZj(1,0) +
                      dphiRZi(1,0) * lambda/R * phi(jn,0);

      double term11 = dphiRZi(0,0) * mi * dphiRZj(0,0) +
                      dphiRZi(1,0) * (lambda + 2.*mi) * dphiRZj(1,0);

      ek(2*in,2*jn)     += weight * R2PI * term00;
      ek(2*in,2*jn+1)   += weight * R2PI * term01;
      ek(2*in+1,2*jn)   += weight * R2PI * term10;
      ek(2*in+1,2*jn+1) += weight * R2PI * term11;
    }
  }
}

void TPZElasticityAxiMaterial::ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc)
{
  TPZFMatrix &phi = data.phi;

  const REAL BIGNUMBER  = TPZMaterial::gBigNumber;

  int phr = phi.Rows();
  short in,jn;
  REAL v2[2];
  v2[0] = bc.Val2()(0,0);
  v2[1] = bc.Val2()(1,0);

  /// R = Dot[{data.x - origin},{AxisR}]   ***because AxisR is already normalized!
  REAL R = (data.x[0] - f_Origin[0])*f_AxisR[0] + (data.x[1] - f_Origin[1])*f_AxisR[1] + (data.x[2] - f_Origin[2])*f_AxisR[2];

  int s = (R > 0) ? 1:-1;
  R = fabs(R);
  double R2PI = 2. * M_PI * R;
  
  static REAL accum1 = 0., accum2 = 0.;

  switch (bc.Type())
  {
      case 0 :// Dirichlet condition
      {
          for(in = 0 ; in < phr; in++)
          {
            ef(2*in,0)   += BIGNUMBER * v2[0]*  // x displacement
                phi(in,0) * R2PI * weight;        // forced v2 displacement

            ef(2*in+1,0) += BIGNUMBER * v2[1]*// x displacement
                phi(in,0) * R2PI * weight;        // forced v2 displacement

            for (jn = 0 ; jn < phi.Rows(); jn++)
            {
              ek(2*in,2*jn)     += BIGNUMBER * phi(in,0) * phi(jn,0) * R2PI * weight;
              ek(2*in+1,2*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * R2PI * weight;
            }
          }
      }
      break;

      case 1 :// Neumann condition
      {
          for(in = 0 ; in < phi.Rows(); in++)
          {           // componentes da tracao na direcao de v2
            ef(2*in,0)   += v2[0] * phi(in,0) * R2PI * weight;   // tracao em x  (ou pressao)
            ef(2*in+1,0) += v2[1] * phi(in,0) * R2PI * weight; // tracao em y (ou pressao) , nula se n� h
          }      // ou deslocamento nulo  v2 = 0
          accum1 += v2[1] * R2PI *weight;
#ifdef LOG4CXX
          {
            std::stringstream sout;
            sout << "Accumulated force 1 " << accum1;
            LOGPZ_DEBUG(logger,sout.str());
          }
#endif
      }
      break;

      case 2 :// condicao mista
      {
          for(in = 0 ; in < phi.Rows(); in++)
          {
            ef(2*in, 0)   += v2[0] * phi(in, 0) * R2PI * weight;   // Neumann , Sigmaij
            ef(2*in+1, 0) += v2[1] * phi(in, 0) * R2PI * weight; // Neumann

            for (jn = 0 ; jn < phi.Rows(); jn++)
            {
              ek(2*in,2*jn)     += bc.Val1()(0,0) * phi(in,0) * phi(jn,0) * R2PI * weight; // peso de contorno => integral de contorno
              ek(2*in+1,2*jn)   += bc.Val1()(1,0) * phi(in,0) * phi(jn,0) * R2PI * weight;
              ek(2*in+1,2*jn+1) += bc.Val1()(1,1) * phi(in,0) * phi(jn,0) * R2PI * weight;
              ek(2*in,2*jn+1)   += bc.Val1()(0,1) * phi(in,0) * phi(jn,0) * R2PI * weight;
            }
          }   // este caso pode reproduzir o caso 0 quando o deslocamento
      }
      break;

      case 3 :// Neumann condition - Normal rotacionada 90 graus no sentido horário em relacao ao vetor axes (1D)
      {
          REAL Nxy[2],Nrz[2];
          Nxy[0] =  data.axes(0,1);
          Nxy[1] = -data.axes(0,0);
          Nrz[0] = s*(Nxy[0]*f_AxisR[0] + Nxy[1]*f_AxisR[1]);
          Nrz[1] = Nxy[0]*f_AxisZ[0] + Nxy[1]*f_AxisZ[1];

//           #ifdef LOG4CXX
//           {
//               std::stringstream sout;
//               sout << "CoordX: " << data.x << " R= " << R << " Nrz = " << Nrz[0] << "," << Nrz[1] << endl;
//               LOGPZ_DEBUG(logger,sout.str());
//           }
//           #endif

          for(in = 0 ; in < phi.Rows(); in++)
          {           // componentes da tracao normal ao contorno
            ef(2*in,0)   += v2[0] * Nrz[0] * phi(in,0) * R2PI * weight;   // tracao em x  (ou pressao)
            ef(2*in+1,0) += v2[0] * Nrz[1] * phi(in,0) * R2PI * weight; // tracao em y (ou pressao) , nula se n� h
          }      // ou deslocamento nulo  v2 = 0
          accum2 += v2[0] * Nrz[1] * R2PI *weight;
#ifdef LOG4CXX
          {
            std::stringstream sout;
            sout << "Accumulated force 2 " << accum2;
            LOGPZ_DEBUG(logger,sout.str());
          }
#endif
      }
      break;
  }      // �nulo introduzindo o BIGNUMBER pelos valores da condicao
}         // 1 Val1 : a leitura �00 01 10 11

/** returns the variable index associated with the name*/
int TPZElasticityAxiMaterial::VariableIndex(const std::string &name)
{
  if(!strcmp("Eigenvector1",name.c_str()))     return 9;
  if(!strcmp("Eigenvector2",name.c_str()))     return 1;
  if(!strcmp("Eigenvector3",name.c_str()))     return 2;
  if(!strcmp("Sigmarr",name.c_str()))       return 3;
  if(!strcmp("Sigmazz",name.c_str()))       return 4;
  if(!strcmp("Sigmatt",name.c_str()))       return 5;
  if(!strcmp("Taurz",name.c_str()))         return 6;
  if(!strcmp("displacement",name.c_str()))  return 7;
  if(!strcmp("MohrCoulomb",name.c_str()))  return 8;

  return TPZMaterial::VariableIndex(name);
  return -1;
}

/**returns the number of variables associated with the
   variable indexed by var, var is obtained by calling VariableIndex*/
int TPZElasticityAxiMaterial::NSolutionVariables(int var)
{
  switch(var) {
  case 9:
    return 3;
  case 1:
    return 3;
  case 2:
    return 3;
  case 3:
    return 1;
  case 4:
    return 1;
  case 5:
    return 1;
  case 6:
    return 1;
  case 7:
    return 3;
    case 8:
      return 1;
  default:
    return TPZMaterial::NSolutionVariables(var);
    return 0;
  }
}

/** returns the solution associated with the var index based
    on the finite element approximation*/
void TPZElasticityAxiMaterial::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
  if(var == 0) 
  {
    TPZMaterial::Solution(data,var,Solout);
    return;
  }
  TPZFMatrix &axes = data.axes;
  TPZVec<REAL> &SolAxes = data.sol;
  TPZFMatrix &DSolAxes = data.dsol;

  double R = (data.x[0] - f_Origin[0])*f_AxisR[0] + (data.x[1] - f_Origin[1])*f_AxisR[1] + (data.x[2] - f_Origin[2])*f_AxisR[2];
  int s = (R>0) ? 1:-1;
  R = fabs(R);

  double axis0DOTr = 0., axis1DOTr = 0., axis0DOTz = 0., axis1DOTz = 0.;
  for(int pos = 0; pos < 3; pos++)
  {
    axis0DOTr += axes.GetVal(0,pos) * f_AxisR[pos] * s;
    axis1DOTr += axes.GetVal(1,pos) * f_AxisR[pos] * s;
    axis0DOTz += axes.GetVal(0,pos) * f_AxisZ[pos];
    axis1DOTz += axes.GetVal(1,pos) * f_AxisZ[pos];
  }
  TPZFNMatrix<9> DSolrz(2,2,0.);
  DSolrz.PutVal(0,0, DSolAxes(0,0)*axis0DOTr + DSolAxes(1,0)*axis1DOTr );
  DSolrz.PutVal(0,1, DSolAxes(0,0)*axis0DOTz + DSolAxes(1,0)*axis1DOTz );
  DSolrz.PutVal(1,0, DSolAxes(0,1)*axis0DOTr + DSolAxes(1,1)*axis1DOTr );
  DSolrz.PutVal(1,1, DSolAxes(0,1)*axis0DOTz + DSolAxes(1,1)*axis1DOTz );


#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << "Point " << data.x << std::endl;
    sout << "Solution " << data.sol << std::endl;
    DSolrz.Print("Derivatives of the solution\n",sout);
    sout << "Radius " << R << std::endl;
    LOGPZ_DEBUG(logger,sout.str());
  }
#endif

  ///Infinitesimal Tensor
  TPZFNMatrix<9> Einf(3,3,0.);
  Einf.PutVal(0,0,DSolrz(0,0)); 
  Einf.PutVal(0,1,0.5*(DSolrz(0,1) + DSolrz(1,0)));
  Einf.PutVal(1,0,0.5*(DSolrz(0,1) + DSolrz(1,0))); 
  Einf.PutVal(1,1,DSolrz(1,1));
  Einf.PutVal(2,2,data.sol[0]/R);

#ifdef LOG4CXX
  {
    std::stringstream sout;
    Einf.Print("Deformation tensor",sout);
    LOGPZ_DEBUG(logger,sout.str());
  }
#endif

  double lambda = -((fE*fnu)/((1. + fnu)*(2.*fnu-1.)));
  double mi =  fE/(2.*(1. + fnu));
  double trE = Einf(0,0) + Einf(1,1) + Einf(2,2);

#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << "E = " << fE << " fnu = " << fnu << " mi = " << mi << " lambda = " << lambda << " trE = " << trE;
    LOGPZ_DEBUG(logger,sout.str());
  }
#endif

  ///Stress Tensor
  TPZFNMatrix<9> T(3,3,0.);
  double cte = lambda*trE;

  //T = lambda.tr(E).I + 2.mi.E
  T.PutVal(0,0,cte + 2.*mi*Einf(0,0)); 
  T.PutVal(0,1,      2.*mi*Einf(0,1)); 
  T.PutVal(0,2,      2.*mi*Einf(0,2));
  T.PutVal(1,0,      2.*mi*Einf(1,0)); 
  T.PutVal(1,1,cte + 2.*mi*Einf(1,1)); 
  T.PutVal(1,2,      2.*mi*Einf(1,2));
  T.PutVal(2,0,      2.*mi*Einf(2,0)); 
  T.PutVal(2,1,      2.*mi*Einf(2,1)); 
  T.PutVal(2,2,cte + 2.*mi*Einf(2,2));

#ifdef LOG4CXX
  {
    std::stringstream sout;
    T.Print("Stress tensor",sout);
    LOGPZ_DEBUG(logger,sout.str());
  }
#endif

  switch(var)
  {
      case 9: //Solout = 1stEigenvalue * {1stEigenvector}
      {
        int NumIt = 1000;
        REAL tol = 1.E-5;
        TPZVec<REAL> EigValues(3,0.);
        TPZFNMatrix<9> EigVectors(3,3,0.);
        bool EigenWorks;
        EigenWorks = T.SolveEigensystemJacobi(NumIt, tol, EigValues, EigVectors);
        if(EigenWorks)
        {
          Solout.Resize(3);
          for(int i = 0; i < 3; i++) Solout[i] = EigValues[0] * EigVectors(0,i);
        }
        else cout << "TPZElasticityAxiMaterial::Solution Error -> case 0\n";
#ifdef LOG4CXX
        {
          std::stringstream sout;
          sout << "First eigenvector " << std::endl;
          sout << Solout;
          LOGPZ_DEBUG(logger,sout.str());
        }
#endif
      }
      break;

      case 1: //Solout = 2ndEigenvalue * {2ndEigenvector}
      {
        int NumIt = 1000;
        REAL tol = 1.E-5;
        TPZVec<REAL> EigValues(3,0.);
        TPZFNMatrix<9> EigVectors(3,3,0.);
        bool EigenWorks;
        EigenWorks = T.SolveEigensystemJacobi(NumIt, tol, EigValues, EigVectors);
        if(EigenWorks)
        {
          Solout.Resize(3);
          for(int i = 0; i < 3; i++) Solout[i] = EigValues[1] * EigVectors(1,i);
        }
        else cout << "TPZElasticityAxiMaterial::Solution Error -> case 1\n";
#ifdef LOG4CXX
        {
          std::stringstream sout;
          sout << "Second eigenvector " << std::endl;
          sout << Solout;
          LOGPZ_DEBUG(logger,sout.str());
        }
#endif
      }
      break;

      case 2: //Solout = 3rdEigenvalue * {3rdEigenvector}
      {
        int NumIt = 1000;
        REAL tol = 1.E-5;
        TPZVec<REAL> EigValues(3,0.);
        TPZFNMatrix<9> EigVectors(3,3,0.);
        bool EigenWorks;
        EigenWorks = T.SolveEigensystemJacobi(NumIt, tol, EigValues, EigVectors);
        if(EigenWorks)
        {
          Solout.Resize(3);
          for(int i = 0; i < 3; i++) Solout[i] = EigValues[2] * EigVectors(2,i);
        }
        else cout << "TPZElasticityAxiMaterial::Solution Error -> case 2\n";
#ifdef LOG4CXX
        {
          std::stringstream sout;
          sout << "Third eigenvector " << std::endl;
          sout << Solout;
          LOGPZ_DEBUG(logger,sout.str());
        }
#endif
      }
      break;

      case 3: //Solout = Sigma_r,r
      {
          Solout.Resize(1);
          Solout[0] = T(0,0);
#ifdef LOG4CXX
          {
            std::stringstream sout;
            sout << "Sigma r " << T(0,0);
            LOGPZ_DEBUG(logger,sout.str());
          }
#endif
      }
      break;

      case 4: //Solout = Sigma_z,z
      {
        Solout.Resize(1);
        Solout[0] = T(1,1);
#ifdef LOG4CXX
        {
          std::stringstream sout;
          sout << "Sigma z " << T(1,1);
          LOGPZ_DEBUG(logger,sout.str());
        }
#endif
      }
      break;

      case 5: //Solout = Sigma_theta,theta
      {
        Solout.Resize(1);
        Solout[0] = T(2,2);
#ifdef LOG4CXX
        {
          std::stringstream sout;
          sout << "Sigma theta " << T(2,2);
          LOGPZ_DEBUG(logger,sout.str());
        }
#endif
      }
      break;

      case 6: //Solout = Tau_r,z
      {
        Solout.Resize(1);
        Solout[0] = T(0,1);
#ifdef LOG4CXX
        {
          std::stringstream sout;
          sout << "Sigma rz " << T(0,1);
          LOGPZ_DEBUG(logger,sout.str());
        }
#endif
      }
      break;

      case 7: //Solout = Displacement
      {
        Solout.Resize(3);
        Solout[0] = SolAxes[0];
        Solout[1] = SolAxes[1];
        Solout[2] = SolAxes[2];
      }
      break;

      case 8: //MohrCoulomb plasticity criteria
      {
        Solout.Resize(1);
        double i1, i2, i3, j1, j2, j3;

        i1 = T(0,0) + T(1,1) + T(2,2);
        TPZFMatrix T2(3,3,0.);
        T.Multiply(T,T2);

        i2 = 0.5*( i1*i1 - (T2(0,0) + T2(1,1) + T2(2,2)) );

        i3 = T(0,0)*T(1,1)*T(2,2) + T(0,1)*T(1,2)*T(2,0) + T(0,2)*T(1,0)*T(2,1) - T(0,2)*T(1,1)*T(2,0) - T(1,2)*T(2,1)*T(0,0) - T(2,2)*T(0,1)*T(1,0);

        j1 = 0.;
        j2 = 1./3.*(i1*i1 - 3.*i2);
        j3 = 1./27.*(2.*i1*i1*i1 - 9.*i1*i2 + 27.*i3);

        double cos3theta = 3.*sqrt(3.)/2. * j3/(pow(j2,1.5));

        double theta = acos(cos3theta)/3.;

        Solout[0] = 1./3.*i1*sin(f_phi) + sqrt(i2)*sin(theta + M_PI/3.) + sqrt(i2/3.)*cos(theta + M_PI/3.)*sin(f_phi) - f_c*cos(f_phi);
#ifdef LOG4CXX
        {
          std::stringstream sout;
          sout << "Criterio de Mohr Coulomb \ni1 " << i1 << " i2 " << i2 <<
              " i3 " << i3 << " \nj1 " << j1 << " j2 " << j2 << " j3 " << j3
              << " \nf_phi " << f_phi << 
              "  f_c " << f_c << " theta " << theta << " MohrCoul " << Solout[0];
          LOGPZ_DEBUG(logger,sout.str());
        }
#endif
      }
      break;

      default:
      {
        cout << "TPZElasticityAxiMaterial::Solution Error -> default\n";
        TPZMaterial::Solution(data.sol,data.dsol,data.axes,var,Solout);
      }
      break;
  }
}

void TPZElasticityAxiMaterial::Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux)
{
  if(fabs(axes(2,0)) >= 1.e-6 || fabs(axes(2,1)) >= 1.e-6)
  {
    cout << "TPZElasticityAxiMaterial::Flux only serves for xy configuration\n";
    axes.Print("axes");
  }
}

void TPZElasticityAxiMaterial::Errors(TPZVec<REAL> &x,TPZVec<REAL> &u, TPZFMatrix &dudaxes, TPZFMatrix &axes, TPZVec<REAL> &flux,
				      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values)
{
  values[0] = 0.;
  TPZManVector<REAL> sigma(3,0.),sigma_exact(3,0.);
  REAL sigx,sigy,sigxy,gamma;

  TPZFNMatrix<9> du;///du = dudx
  TPZAxesTools::Axes2XYZ(dudaxes, du, axes);

  //tensoes aproximadas : uma forma
  gamma = du(1,0)+du(0,1);
  sigma[0] = fEover1MinNu2*(du(0,0)+fnu*du(1,1));
  sigma[1] = fEover1MinNu2*(fnu*du(0,0)+du(1,1));
  sigma[2] = fE*0.5/(1.+fnu)*gamma;

  TPZMaterialData mydata;
  mydata.sol  = u;
  mydata.dsol = du;
  mydata.axes = axes;
  mydata.x = x;

  //tensoes aproximadas : outra forma
  TPZVec<REAL> sol(1);
  Solution(mydata,5,sol);
  sigma[0] = sol[0];
  Solution(mydata,6,sol);
  sigma[1] = sol[0];
  Solution(mydata,8,sol);
  sigma[2] = sol[0];

  //exata
  gamma = du_exact(1,0)+du_exact(0,1);
  sigma_exact[0] = fEover1MinNu2*(du_exact(0,0)+fnu*du_exact(1,1));
  sigma_exact[1] = fEover1MinNu2*(fnu*du_exact(0,0)+du_exact(1,1));
  sigma_exact[2] = fE*0.5/(1.+fnu)*gamma;
  sigx  = (sigma[0] - sigma_exact[0]);
  sigy  = (sigma[1] - sigma_exact[1]);
  sigxy = (sigma[2] - sigma_exact[2]);
  //values[0] = calculo do erro estimado em norma Energia
  values[0] = fE*(sigx*sigx + sigy*sigy + 2*fnu*sigx*sigy)/(1-fnu*fnu);
  values[0] = (values[0] + .5*fE*sigxy*sigxy/(1+fnu));

  //values[1] : erro em norma L2 em tensoes
  //values[1] = sigx*sigx + sigy*sigy + sigxy*sigxy;

  //values[1] : erro em norma L2 em deslocamentos
  values[1] = pow(fabs(u[0] - u_exact[0]),2.0)+pow(fabs(u[1] - u_exact[1]),2.0);

  //values[2] : erro estimado
  values[2] = 0.;
}


TPZElasticityAxiMaterial::TPZElasticityAxiMaterial(const TPZElasticityAxiMaterial &copy) : 
		TPZMaterial(copy), fE(copy.fE),
        fnu(copy.fnu), fEover21PlusNu(copy.fEover21PlusNu),
        fEover1MinNu2(copy.fEover1MinNu2),f_c(copy.f_c),f_phi(copy.f_phi),
		f_Origin(copy.f_Origin),f_AxisZ(copy.f_AxisZ),f_AxisR(copy.f_AxisR)
{
	ff[0] = copy.ff[0];
	ff[1] = copy.ff[1];
	ff[2] = copy.ff[2];
}

int TPZElasticityAxiMaterial::ClassId() const
{
  return TPZELASTICITYMATERIALID;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZElasticityAxiMaterial,TPZELASTICITYMATERIALID>;
#endif

void TPZElasticityAxiMaterial::Read(TPZStream &buf, void *context)
{
  TPZMaterial::Read(buf,context);
  buf.Read(&fE,1);
  buf.Read(&fnu,1);
  buf.Read(&fEover21PlusNu,1);
  buf.Read(&fEover1MinNu2,1);

  buf.Read(ff,3);
}

void TPZElasticityAxiMaterial::Write(TPZStream &buf, int withclassid)
{
  TPZMaterial::Write(buf,withclassid);
  buf.Write(&fE,1);
  buf.Write(&fnu,1);
  buf.Write(&fEover21PlusNu,1);
  buf.Write(&fEover1MinNu2,1);

  buf.Write(ff,3);
}
