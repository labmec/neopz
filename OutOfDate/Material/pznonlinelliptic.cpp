/**
 * @file
 * @brief DEPRECATED FILE. Contains the implementations of the TPZNonLinElliptic methods.
 */
// 
// C++ Implementation: pznonlinelliptic
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pznonlinelliptic.h"
#include "pzbndcond.h"

TPZNonLinElliptic::TPZNonLinElliptic(int id, int dimension) : TPZMaterial(id)
{
  this->fDim = dimension;
  this->fConvDir.Resize( dimension );
  this->fConvDir.Fill( 0.0 );
}


TPZNonLinElliptic::~TPZNonLinElliptic(){}

void TPZNonLinElliptic::SetParameters(REAL D, TPZVec<REAL> &V, REAL Sigma, REAL LambdaDivK, REAL F)
{
  this->fCoeffD = D;
  this->fConvDir = V;
  this->fSigma = Sigma;
  this->fLambdaDivK = LambdaDivK;
  this->fSource = F;
}
    
void TPZNonLinElliptic::GetParameters(REAL &D, TPZVec<REAL> &V, REAL &Sigma, REAL &LambdaDivK, REAL &F)
{
  D = this->fCoeffD;
  V = this->fConvDir;
  Sigma = this->fSigma;
  LambdaDivK = this->fLambdaDivK;
  F = this->fSource;
}

void TPZNonLinElliptic::Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,
                            TPZVec<REAL> &sol, TPZFMatrix &dsol,
                            REAL weight,TPZFMatrix &axes,
                            TPZFMatrix &phi, TPZFMatrix &dphi,
                            TPZFMatrix &ek, TPZFMatrix &ef)
{
  int phr = phi.Rows();
  if(fForcingFunction)            // phi(in, 0) = phi_in
  {
    TPZManVector<REAL, 10> res(1);
    fForcingFunction(x,res);       // dphi(i,j) = dphi_j/dxi
    fSource = res[0];
  }
  TPZManVector<REAL, 3> ConvDirAx(3, 0.0);
    switch(fDim)
    {
      case 1:
        ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
        break;
      case 2:
        ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
        ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
        break;
      case 3:
        ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
        ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
        ConvDirAx[2] = axes(2,0)*fConvDir[0]+axes(2,1)*fConvDir[1]+axes(2,2)*fConvDir[2];
        break;
      default:
        PZError << "TPZNonLinElliptic::Contribute dimension error " << fDim << endl;
    }
  for( int in = 0; in < phr; in++ )
  {
    int kd;
    ef(in, 0) += weight * fSource*phi(in,0);
    for( int jn = 0; jn < phr; jn++ )
    {
      for(kd=0; kd<fDim; kd++)
      {
        ek(in,jn) += weight * (
            fCoeffD * ( dphi(kd,in) * dphi(kd,jn) )
          + ConvDirAx[kd]* dphi(kd,jn) * phi(in) );
      }
      ek(in,jn) += weight * ( fSigma * phi(jn) * phi(in) );
    }
  }
  
  //Termos n� lineares
  for( int in = 0; in < phr; in++ )
  {
    int kd;
    for(kd = 0; kd < fDim; kd++)
    {
      ef(in, 0) += weight * (fCoeffD * ( dphi(kd,in) * dsol(kd) )
            + ConvDirAx[kd]* dsol(kd) * phi(in) );
    }
    ef(in,0) += weight * ( ( fSigma * sol[0] * phi(in) ) + fLambdaDivK * sol[0] * sol[0] * phi(in) );
          
          
    for( int jn = 0; jn < phr; jn++ )
    {
      ek(in,jn) += weight * fLambdaDivK * 2.0 * sol[0] * phi(in)*  phi(jn);
    }
  }  
  
  
}

void TPZNonLinElliptic::Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,
                            TPZVec<REAL> &sol, TPZFMatrix &dsol,
                            REAL weight,TPZFMatrix &axes,
                            TPZFMatrix &phi, TPZFMatrix &dphi,
                            TPZFMatrix &ef)
{
  TPZFNMatrix<100> pseudo_ek( ef.Rows(), ef.Rows() );
  this->Contribute(x, jacinv, sol, dsol, weight, axes, phi, dphi, pseudo_ek, ef);
}

void TPZNonLinElliptic::ContributeBC(TPZVec<REAL> &x, TPZVec<REAL> &sol,
                            REAL weight, TPZFMatrix &axes,
                            TPZFMatrix &phi, TPZFMatrix &ek,
                            TPZFMatrix &ef, TPZBndCond &bc)
{
  int phr = phi.Rows();
  int in,jn;
  REAL v2[1], v1[1];
  v2[0] = bc.Val2()(0,0);
  v1[0] = bc.Val1()(0,0);

//linear  
  switch (bc.Type())
  {
  
  case 0 :			// Dirichlet condition
  for(in = 0 ; in < phr; in++)
  {
    ef(in,0) += gBigNumber * v2[0] * phi(in,0) * weight;
    for (jn = 0 ; jn < phr; jn++)
    {
      ek(in,jn) += gBigNumber * phi(in,0) * phi(jn,0) * weight;
    }   
  }
  break;
  case 1 :			// Neumann condition
  for(in = 0 ; in < phi.Rows(); in++)
  {
    ef(in,0) += v2[0] * phi(in,0) * weight;
  }
  break;  
  
  case 2 :		// condi�o mista
  //v1[0] = p
  //v2[0] = gamma
  for(in = 0 ; in < phi.Rows(); in++)
  {
    ef(in, 0) += v2[0] * phi(in, 0) * weight;
    for (jn = 0 ; jn < phi.Rows(); jn++)
    {
      ek(in,jn) += v1[0] * phi(in,0) * phi(jn,0) * weight;
    }
  }
  break;
  }

// termos n� lineares
  REAL aux = 0.;
  switch (bc.Type())
  {
  
  case 0 :  // Dirichlet condition
//Nothing to be done!
  break; 
  
  case 1 :			// Neumann condition
//Nothing to be done!  
  break;
  
  
  case 2 :		// condi�o mista
  for(in = 0 ; in < phi.Rows(); in++)
  {
    ef(in, 0) += v1[0] * sol[0] * phi(in, 0) * weight;
  }
  break;
  }
}
  
int TPZNonLinElliptic::VariableIndex(const std::string &name){
  if(!strcmp("Solution",name)) return 1;
  cout << "TPZNonLinElliptic::VariableIndex Error\n";
  return -1;
}

int TPZNonLinElliptic::NSolutionVariables(int var){
  if(var == 1) return 1;
  cout << "TPZNonLinElliptic::NSolutionVariables Error\n";
  return 0;
}

void TPZNonLinElliptic::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &/*axes*/,int var,TPZVec<REAL> &Solout){
  if(var == 1){ 
    Solout[0] = Sol[0];
    return;
  }
  cout << "TPZNonLinElliptic::Solution Error\n";  
}

void TPZNonLinElliptic::Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
                               TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
                               TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values) {

  TPZVec<REAL> sol(1);
  Solution(u,dudx,axes,1,sol);
    
  //values[1] : eror em norma L2
  values[1]  = (sol[0] - u_exact[0])*(sol[0] - u_exact[0]);
  
  //values[2] : erro em semi norma H1
  values[2] = 0.;

  //values[0] : erro em norma H1 <=> norma Energia
  values[0]  = 0.;
}
