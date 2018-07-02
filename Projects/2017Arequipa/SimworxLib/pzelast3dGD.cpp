// $Id: pzelast3d.cpp,v 1.16 2010-09-06 14:50:47 phil Exp $

#include "pzmaterialid.h"
#include "pzelast3dGD.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
// #include "pztempmat.h"
#include "pzmanvector.h"
#include <math.h>
#include <fstream>
//#include "Dialogs.hpp"

REAL TPZElasticity3DGD::gTolerance = 1.e-11;

TPZElasticity3DGD::TPZElasticity3DGD(int nummat, REAL E, REAL poisson,
                                 TPZVec<REAL> &force,
                                 REAL preStressXX, REAL preStressYY, REAL preStressZZ) : TPZDiscontinuousGalerkin(nummat) {
  this->fE = E;
  this->fPoisson = poisson;
#ifdef DEBUG
  if (force.NElements() != 3)
    PZError << __PRETTY_FUNCTION__ << " - error!" << std::endl;
#endif
  int i;
  this->fForce.Resize(3);
  for (i = 0; i < 3; i++) {
    this->fForce[i] = force[i];
  }
  // Default directions is {1,0,0}
  this->fPostProcessDirection.Resize(3);
  this->fPostProcessDirection.Fill(0.);
  this->fPostProcessDirection[0] = 1.;
  this->SetYieldingStress(1.);

  fPreStress.Resize(3);
  fPreStress[0] = preStressXX;
  fPreStress[1] = preStressYY;
  fPreStress[2] = preStressZZ;

//  this->fPlasticModel = NULL;
//  this->fDamageModel = NULL;
//  this->fCohesiveDamageModel = NULL;
} // method

TPZElasticity3DGD::TPZElasticity3DGD(int nummat) : TPZDiscontinuousGalerkin(nummat), fE(0.),
    fPoisson(0.), fForce(3, 0.), fPostProcessDirection(3, 0.), fFy(0.),
    fPreStress(3, 0.) /*, fPlasticModel(), fDamageModel(), fCohesiveDamageModel()*/ {
}

TPZElasticity3DGD::TPZElasticity3DGD() : TPZDiscontinuousGalerkin(), fE(0.), fPoisson(0.),
    fForce(3, 0.), fPostProcessDirection(3, 0.), fFy(0.), fPreStress(3, 0.) /*,
    fPlasticModel(), fDamageModel(), fCohesiveDamageModel()*/ {
}

TPZElasticity3DGD::~TPZElasticity3DGD() {
}

TPZElasticity3DGD::TPZElasticity3DGD(const TPZElasticity3DGD &cp)
    : TPZDiscontinuousGalerkin(cp), fE(cp.fE), fPoisson(cp.fPoisson), fForce(cp.fForce),
    fPostProcessDirection(cp.fPostProcessDirection), fFy(cp.fFy),
    fPreStress(cp.fPreStress) {
/*  if(cp.fPlasticModel){
    this->fPlasticModel = cp.fPlasticModel->Clone();
  }
  else{
    this->fPlasticModel = NULL;
  }

  if(cp.fDamageModel){
    this->fDamageModel = cp.fDamageModel->Clone();
  }
  else{
    this->fDamageModel = NULL;
  }

  if(cp.fCohesiveDamageModel){
    this->fCohesiveDamageModel = cp.fCohesiveDamageModel->Clone();
  }
  else{
    this->fCohesiveDamageModel = NULL;
  }
  */
}

void TPZElasticity3DGD::Print(std::ostream & out) {
  out << "\nTPZElasticity3DGD material:\n";
  out << "\tfE       = " << this->fE << std::endl;
  out << "\tfPoisson = " << this->fPoisson << std::endl;
  out << "\tfForce   = " << this->fForce << std::endl;
  out << "\tBase class print\n";
  TPZDiscontinuousGalerkin::Print(out);
  out << "End of TPZElasticity3DGD::Print\n";
}

void TPZElasticity3DGD::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){

	TPZFMatrix<STATE> &dphi = data.dphix;
	TPZFMatrix<STATE> &phi = data.phi;
  TPZManVector<REAL, 3>&x = data.x;

  const int phr = phi.Rows();
  if (this->fForcingFunction) {
    this->fForcingFunction->Execute(x, fForce);
  }

  //tensor de tensão
  TPZFNMatrix<9> StressTensor(3,3);
  bool yielding;
  REAL damage = 0.;
  this->ComputeStressTensor(StressTensor, yielding, damage, data.dsol[0], data.gelElId, data.intLocPtIndex, data.HSize);

  for (int in = 0; in < phr; in++) {
    for (int kd = 0; kd < 3; kd++) {
       ef(in * 3 + kd, 0) +=
          weight * (fForce[kd] * phi(in, 0) - fPreStress[kd] * dphi(kd, in)
           - (StressTensor(kd,0)*dphi(0,in)+StressTensor(kd,1)*dphi(1,in)+StressTensor(kd,2)*dphi(2,in))
           );
    } // kd
  } // in

}//Contribute residual

void TPZElasticity3DGD::ContributeOriginal(TPZMaterialData &data, REAL weight,
	TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
	TPZFMatrix<STATE> &dphi = data.dphix;
	TPZFMatrix<STATE> &phi = data.phi;
  TPZManVector<REAL, 3>&x = data.x;

  const int phr = phi.Rows();
  if (this->fForcingFunction) {
    this->fForcingFunction->Execute(x, fForce);
  }

  //tensor de tensão
  TPZFNMatrix<9> StressTensor(3,3);
  bool yielding;
  REAL damage = 0.;
  this->ComputeStressTensor(StressTensor, yielding, damage, data.dsol[0], data.gelElId, data.intLocPtIndex, data.HSize);

  TPZFNMatrix<9>Deriv(3,3,0.);

  const REAL damageFactor = (1. - damage);
  const REAL limite = 1./10.;// 1./10.; // Acho que esse limite deve depender da ponderacao do alpha do Contribute
  const REAL factor = (damageFactor > limite) ? damageFactor : limite;
  const REAL E = this->fE  * factor;
//  const REAL nu = (yielding == false) ? this->fPoisson : 0.49;
  const REAL nu = this->fPoisson;

  const REAL C1 = E / (2. + 2. * nu);
  const REAL C2 = E * nu / (-1. + nu + 2. * nu * nu);
  const REAL C3 = E * (nu - 1.) / (-1. + nu + 2. * nu * nu);

  int in;
  for (in = 0; in < phr; in++) {
    int kd;
    for (kd = 0; kd < 3; kd++) {
       ef(in * 3 + kd, 0) +=
          weight * (fForce[kd] * phi(in, 0) - fPreStress[kd] * dphi(kd, in)
           - (StressTensor(kd,0)*dphi(0,in)+StressTensor(kd,1)*dphi(1,in)+StressTensor(kd,2)*dphi(2,in))
           );
    } // kd
    REAL val;
    for (int jn = 0; jn < phr; jn++) {
      // Compute Deriv matrix
      for (int ud = 0; ud < 3; ud++) {
        for (int vd = 0; vd < 3; vd++) {
          Deriv(vd, ud) = dphi(vd, in) * dphi(ud, jn);
        } // ud
      } // vd

      // First equation Dot[Sigma1, gradV1]
      val = (Deriv(1, 1) + Deriv(2, 2)) * C1 + Deriv(0, 0) * C3;
      ek(in * 3 + 0, jn * 3 + 0) +=  weight * val;

      val = Deriv(1, 0) * C1 - Deriv(0, 1) * C2;
      ek(in * 3 + 0, jn * 3 + 1) +=  weight * val;

      val = Deriv(2, 0) * C1 - Deriv(0, 2) * C2;
      ek(in * 3 + 0, jn * 3 + 2) +=  weight * val;

      // Second equation Dot[Sigma2, gradV2]
      val = Deriv(0, 1) * C1 - Deriv(1, 0) * C2;
      ek(in * 3 + 1, jn * 3 + 0) +=  weight * val;

      val = (Deriv(0, 0) + Deriv(2, 2)) * C1 + Deriv(1, 1) * C3;
      ek(in * 3 + 1, jn * 3 + 1) +=  weight * val;

      val = Deriv(2, 1) * C1 - Deriv(1, 2) * C2;
      ek(in * 3 + 1, jn * 3 + 2) +=  weight * val;

      // Third equation Dot[Sigma3, gradV3]
      val = Deriv(0, 2) * C1 - Deriv(2, 0) * C2;
      ek(in * 3 + 2, jn * 3 + 0) +=  weight * val;

      val = Deriv(1, 2) * C1 - Deriv(2, 1) * C2;
      ek(in * 3 + 2, jn * 3 + 1) +=  weight * val;

      val = (Deriv(0, 0) + Deriv(1, 1)) * C1 + Deriv(2, 2) * C3;
      ek(in * 3 + 2, jn * 3 + 2) +=  weight * val;

    } // jn
  } // in

#ifdef DEBUG
  if (!ek.VerifySymmetry(1.e-8)){
    PZError << __PRETTY_FUNCTION__ << "\nERROR - NON SYMMETRIC MATRIX" << std::endl;
  }
#endif
} // method

//calculando jacobiana por diferenca finita
void TPZElasticity3DGD::EstimateJacobian(TPZFMatrix<STATE> &dsol, int elindex, int elIntPoint, const REAL elementCharacteristicSize,
                                       TPZVec< TPZVec<REAL> > &jac){
  TPZManVector<REAL,6> StrainVec(6), StressVec(6);
  this->ComputeStrainVector(StrainVec, dsol);
  bool yielding;
  REAL damage;
  this->ComputeStressVector(StressVec,yielding,damage,StrainVec,elindex,elIntPoint,elementCharacteristicSize);
  jac.Resize(6);
  for(int k = 0; k < 6; k++) jac[k].Resize(6);
  for(int i = 0; i < 6; i++){
    TPZManVector<REAL,6> lcStrainVec(6), lcStressVec(6);
    lcStrainVec = StrainVec;
    const REAL deltaEps = 1e-7;
    lcStrainVec[i] += deltaEps;
    this->ComputeStressVector(lcStressVec,yielding,damage,lcStrainVec,elindex,elIntPoint,elementCharacteristicSize);
    for(int j = 0; j < 6; j++){
      jac[j][i] = (lcStressVec[j] - StressVec[j])/deltaEps;
    }//j
  }

///Agora que estou combinando com a matriz elastica, nao precisa mais disso. Deixa dar zero que a elastica resolve
/*  const REAL min = this->fE * 1e-3;
  for(int i = 0; i < 3; i++){ //{Sx,Sy e Sz apenas}
    if(jac[i][i] < min){
     jac[i][i] = min;
    }
  } */

  return;

#ifdef DEBUG
  bool fail = false;
{
  REAL sum = 0.;
  for(int i = 0; i < 6; i++) for(int j = 0; j < 6; j++) sum += fabs(jac[j][i] - jac[i][j]);
  if (sum > 1.e-3){
    fail = true;
    {
      std::ofstream myfile("c:\\Temp\\estimatejacobian.txt");
      for(int i = 0; i < 6; i++) {
        for(int j = 0; j < 6; j++) myfile << jac[i][j] << ",\t";
        myfile << "\n";
      }
    }
   PZError << __PRETTY_FUNCTION__ << "\nERROR - NON SYMMETRIC MATRIX" << std::endl;
  }
}
#endif

  //a diferenca finita pode fazer com que jac nao fique simetrico. Faz sentido? Parece que faz: lembrar que dSigmaX/dExy é diferente de dTauXY/dExx por um fator 2
  // ou seja seriam iguais se fosse  dSigmaX/dGammaxy. Ainda nao tenho plena certeza disso, mas ta estranho. A simetria deveria voltar na jacobiana com respeito aos coeficientes
  //multiplicadores das funcoes de forma, eu acho
  //Outro argumento, que acho mais fragil: o Mohr-Coulomb nao associativo nao garante simetria.
  //simetrizando
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < i; j++){
      REAL newval = 0.5*(jac[i][j]+jac[j][i]);
      jac[i][j] = newval;
      jac[j][i] = newval;
    }
  }

  for(int i = 0; i < 3; i++) for(int j = 3; j < 6; j++){
    jac[i][j] = 0.;
    jac[j][i] = 0.;
  }

  for(int i = 3; i < 6; i++){
    for(int j = 3; j < 6; j++){
      if(i != j) jac[i][j] = 0.;
    }
  }

#ifdef DEBUG
  if(fail){
    REAL sum = 0.;
    for(int i = 0; i < 6; i++) for(int j = 0; j < 6; j++) sum += fabs(jac[j][i] - jac[i][j]);
    {
      std::ofstream myfile("c:\\Temp\\estimatejacobian2.txt");
      for(int i = 0; i < 6; i++) {
        for(int j = 0; j < 6; j++) myfile << jac[i][j] << ",\t";
        myfile << "\n";
      }
    }
  }
#endif

}

#ifdef _AUTODIFF
void TPZElasticity3DGD::ContributeDifFinita(TPZMaterialData &data, REAL weight,
                                      TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
  TPZFMatrix<REAL> &dphi = data.dphix;
  TPZFMatrix<REAL> &phi = data.phi;
  TPZManVector<REAL, 3>&x = data.x;

  const int phr = phi.Rows();
  if (this->fForcingFunction) {
    this->fForcingFunction->Execute(x, fForce);
  }

  TPZManVector< TPZVec<REAL>, 6 > jac(6);
  this->EstimateJacobian(data.dsol[0], data.gelElId, data.intLocPtIndex, data.HSize, jac);

  //tensor de tensão
  TPZFNMatrix<9> StressTensor(3,3);
  bool yielding;
  REAL damage = 0.;

  this->ComputeStressTensor(StressTensor, yielding, damage, data.dsol[0], data.gelElId, data.intLocPtIndex, data.HSize);

  TPZManVector<  TFad<6, REAL> , 6 > StressVecFad(6);
  //Sigma xx
  StressVecFad[0] = StressTensor(0,0);
  for(int j = 0; j < 6; j++){
    StressVecFad[0].fastAccessDx(j) = jac[0][j];
  }
  //Sigma yy
  StressVecFad[1] = StressTensor(1,1);
  for(int j = 0; j < 6; j++){
    StressVecFad[1].fastAccessDx(j) = jac[1][j];
  }
  //Sigma zz
  StressVecFad[2] = StressTensor(2,2);
  for(int j = 0; j < 6; j++){
    StressVecFad[2].fastAccessDx(j) = jac[2][j];
  }
  //T xy
  StressVecFad[3] = StressTensor(0,1);
  for(int j = 0; j < 6; j++){
    StressVecFad[3].fastAccessDx(j) = jac[3][j];
  }
  //T xz
  StressVecFad[4] = StressTensor(0,2);
  for(int j = 0; j < 6; j++){
    StressVecFad[4].fastAccessDx(j) = jac[4][j];
  }
  //T yz
  StressVecFad[5] = StressTensor(1,2);
  for(int j = 0; j < 6; j++){
    StressVecFad[5].fastAccessDx(j) = jac[5][j];
  }

  for (int in = 0; in < phr; in++) {

    TFad<6, REAL> Residual1 = StressVecFad[0]*dphi(0,in) //Sxx
                             +StressVecFad[3]*dphi(1,in) //Txy
                             +StressVecFad[4]*dphi(2,in);//Txz
    TFad<6, REAL> Residual2 = StressVecFad[3]*dphi(0,in) //Txy
                             +StressVecFad[1]*dphi(1,in) //Syy
                             +StressVecFad[5]*dphi(2,in);//Tyz
    TFad<6, REAL> Residual3 = StressVecFad[4]*dphi(0,in) //Txz
                             +StressVecFad[5]*dphi(1,in) //Tyz
                             +StressVecFad[2]*dphi(2,in);//Szz

    REAL lhs[3] = {Residual1.val(), Residual2.val(), Residual3.val()};
    for (int kd = 0; kd < 3; kd++) {
       ef(in * 3 + kd, 0) +=
          weight * (fForce[kd] * phi(in, 0) - fPreStress[kd] * dphi(kd, in)
           - (lhs[kd]) );
    } // kd
    REAL val;
    TFad<6, REAL> * Residual[3];
    Residual[0] = &Residual1;
    Residual[1] = &Residual2;
    Residual[2] = &Residual3;
    for (int jn = 0; jn < phr; jn++) {
      for(int kd = 0; kd < 3; kd++){
        // First equation Dot[Sigma1, gradV1]
        // Second equation Dot[Sigma2, gradV2]
        // Third equation Dot[Sigma3, gradV3]
        val = Residual[kd]->dx(0)*1.0*dphi(0,jn) +
              Residual[kd]->dx(3)*0.5*dphi(1,jn) +
              Residual[kd]->dx(4)*0.5*dphi(2,jn) ;//D[R,a_u]
        ek(in * 3 + kd, jn * 3 + 0) +=  weight * val;

        val = Residual[kd]->dx(3)*0.5*dphi(0,jn)+
              Residual[kd]->dx(1)*1.0*dphi(1,jn) +
              Residual[kd]->dx(5)*0.5*dphi(2,jn);   //D[R,a_v]
        ek(in * 3 + kd, jn * 3 + 1) +=  weight * val;

        val = Residual[kd]->dx(4)*0.5*dphi(0,jn)+
              Residual[kd]->dx(5)*0.5*dphi(1,jn)+
              Residual[kd]->dx(2)*1.0*dphi(2,jn);  //D[R,a_w]
        ek(in * 3 + kd, jn * 3 + 2) +=  weight * val;
      }

    } // jn
  } // in

#ifdef DEBUG
  if (!ek.VerifySymmetry(1.e-3)){
    {
      std::ofstream myfile("c:\\Temp\\pzelast3d_contribute.txt");
      ek.Print("ek=",myfile,EMathematicaInput);
    }
    {
      std::ofstream myfile("c:\\Temp\\pzelast3d_contribute_jac.txt");
      for(int i = 0; i < 6; i++) {
        for(int j = 0; j < 6; j++) myfile << jac[i][j] << ",\t";
        myfile << "\n";
      }
    }
    {
      std::ofstream myfile("c:\\Temp\\pzelast3d_contribute_StressVecFad.txt");
      for(int i = 0; i < 6; i++) {
        for(int j = 0; j < 6; j++) myfile << StressVecFad[i].dx(j) << ",\t";
        myfile << "\n";
      }
    }
    PZError << __PRETTY_FUNCTION__ << "\nERROR - NON SYMMETRIC MATRIX" << std::endl;
  }
#endif
} // method

#endif

bool TPZElasticity3DGD::UsesBigNumberForBC(int bcType) const{
  if(bcType == 0) return true;
  if(bcType == 1) return false;
  if(bcType == 2) return true;
  if(bcType == 3) return true;
  if(bcType == 4) return false;
  return false;
}

void TPZElasticity3DGD::ContributeBC(TPZMaterialData &data, REAL weight,
	TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
	TPZFMatrix<STATE> &phi = data.phi;

  const int phr = phi.Rows();
  REAL v2[3];
  v2[0] = bc.Val2()(0, 0);
  v2[1] = bc.Val2()(1, 0);
  v2[2] = bc.Val2()(2, 0);
  TPZFMatrix<STATE> &v1 = bc.Val1();

  switch (bc.Type()) {

  case 0: // Dirichlet condition
    for (int in = 0; in < phr; in++) {
      ef(3 * in + 0, 0) += gBigNumber * (v2[0]-data.sol[0][0]) * phi(in, 0) * weight;
      ef(3 * in + 1, 0) += gBigNumber * (v2[1]-data.sol[0][1]) * phi(in, 0) * weight;
      ef(3 * in + 2, 0) += gBigNumber * (v2[2]-data.sol[0][2]) * phi(in, 0) * weight;

      for (int jn = 0; jn < phr; jn++) {
        ek(3 * in + 0, 3 * jn + 0) +=
            gBigNumber * phi(in, 0) * phi(jn, 0) * weight;
        ek(3 * in + 1, 3 * jn + 1) +=
            gBigNumber * phi(in, 0) * phi(jn, 0) * weight;
        ek(3 * in + 2, 3 * jn + 2) +=
            gBigNumber * phi(in, 0) * phi(jn, 0) * weight;
      } // jn
    } // in
    break;

  case 1: // Neumann condition
    for (int in = 0; in < phi.Rows(); in++) {
      ef(3 * in + 0, 0) += v2[0] * phi(in, 0) * weight;
      ef(3 * in + 1, 0) += v2[1] * phi(in, 0) * weight;
      ef(3 * in + 2, 0) += v2[2] * phi(in, 0) * weight;
    } // in
    break;

  case 2: // Mixed condition

    //Tiago em 19 de outubro de 2015
    //a implementacao abaixo so funciona para Val1 diagonal. Se nao for, os termos de ef estão errados
    //eu constatei o erro, mas nao quis arrumar, pois não teria como testar. Acho que funcionava até hoje,
    //porque essa condicao so era usada para fazer apoios moveis alinhados com x,y e z
#ifdef DEBUG
    if( (bc.Val1().Rows() != 3) || (bc.Val1().Cols() != 3) ) DebugStop();
    for(int i = 0; i < 3; i++){
      for(int j = 0; j < 3; j++){
        if(i != j){
          if(IsZero(bc.Val1()(i,j)) == false){
            DebugStop();
          }
        }
      }
    }
#endif
    for (int in = 0; in < phr; in++) {
		ef(3 * in + 0, 0) += bc.Val1()(0, 0) * (v2[0] - data.sol[0][0]) * phi(in, 0) * weight;
		ef(3 * in + 1, 0) += bc.Val1()(1, 1) * (v2[1] - data.sol[0][1]) * phi(in, 0) * weight;
		ef(3 * in + 2, 0) += bc.Val1()(2, 2) * (v2[2] - data.sol[0][2]) * phi(in, 0) * weight;


      for (int jn = 0; jn < phr; jn++) {
        ek(3 * in + 0, 3 * jn + 0) += bc.Val1()(0,0) * phi(in, 0) * phi(jn, 0) * weight;
        ek(3 * in + 1, 3 * jn + 1) += bc.Val1()(1,1) * phi(in, 0) * phi(jn, 0) * weight;
        ek(3 * in + 2, 3 * jn + 2) += bc.Val1()(2,2) * phi(in, 0) * phi(jn, 0) * weight;
      } // jn
    } // in
    break;

  case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
    for (int in = 0; in < phr; in++) {

		ef(3 * in + 0, 0) += gBigNumber * (0. - data.sol[0][0]) * phi(in, 0) * v2[0] * weight;
		ef(3 * in + 1, 0) += gBigNumber * (0. - data.sol[0][0]) * phi(in, 0) * v2[1] * weight;
		ef(3 * in + 2, 0) += gBigNumber * (0. - data.sol[0][0]) * phi(in, 0) * v2[2] * weight;


      for (int jn = 0; jn < phr; jn++) {
        ek(3 * in + 0, 3 * jn + 0) +=
            gBigNumber * phi(in, 0) * phi(jn, 0) * weight * v2[0];
        ek(3 * in + 1, 3 * jn + 1) +=
            gBigNumber * phi(in, 0) * phi(jn, 0) * weight * v2[1];
        ek(3 * in + 2, 3 * jn + 2) +=
            gBigNumber * phi(in, 0) * phi(jn, 0) * weight * v2[2];
      } // jn
    } // in
    break;

  case 4: // stressField Neumann condition
    for (int in = 0; in < 3; in++){
      v2[in] = -(v1(in, 0) * data.normal[0] + v1(in, 1) * data.normal[1] + v1(in, 2) * data.normal[2]);
    }
    // The normal vector points towards the neighbour. The negative sign is there to
    // reflect the outward normal vector.
    for (int in = 0; in < phi.Rows(); in++) {
      ef(3 * in + 0, 0) += v2[0] * phi(in, 0) * weight;
      ef(3 * in + 1, 0) += v2[1] * phi(in, 0) * weight;
      ef(3 * in + 2, 0) += v2[2] * phi(in, 0) * weight;
    }
    break;

  default:
    {
      PZError << "TPZElastitity3D::ContributeBC error - Wrong boundary condition type" << std::endl;
      DebugStop();
    }
  } // switch
} // method

int TPZElasticity3DGD::VariableIndex(const std::string &name) {

  if (!strcmp("StressVector", name.c_str())) return TPZElasticity3DGD::EStressVector;
  if (!strcmp("Displacement", name.c_str())) return TPZElasticity3DGD::EDisplacement;
  if (!strcmp("state", name.c_str())) return TPZElasticity3DGD::EDisplacement;
  if (!strcmp("DisplacementX", name.c_str())) return TPZElasticity3DGD::EDisplacementX;
  if (!strcmp("DisplacementY", name.c_str())) return TPZElasticity3DGD::EDisplacementY;
  if (!strcmp("DisplacementZ", name.c_str())) return TPZElasticity3DGD::EDisplacementZ;
  if (!strcmp("PrincipalStress", name.c_str())) return TPZElasticity3DGD::EPrincipalStress;
  if (!strcmp("PrincipalStrain", name.c_str())) return TPZElasticity3DGD::EPrincipalStrain;
  if (!strcmp("PrincipalPlasticStrain", name.c_str())) return TPZElasticity3DGD::EPrincipalPlasticStrain;
  if (!strcmp("VonMises", name.c_str())) return TPZElasticity3DGD::EVonMisesStress;
  if (!strcmp("Stress", name.c_str())) return TPZElasticity3DGD::EStress;
  if (!strcmp("Strain", name.c_str())) return TPZElasticity3DGD::EStrain;
  if (!strcmp("Stress1", name.c_str())) return TPZElasticity3DGD::EStress1;
  if (!strcmp("Strain1", name.c_str())) return TPZElasticity3DGD::EStrain1;
  if (!strcmp("NormalStress", name.c_str())) return TPZElasticity3DGD::ENormalStress;
  if (!strcmp("NormalStrain", name.c_str())) return TPZElasticity3DGD::ENormalStrain;
  if (!strcmp("StressX", name.c_str())) return TPZElasticity3DGD::EStressX;
  if (!strcmp("StressY", name.c_str())) return TPZElasticity3DGD::EStressY;
  if (!strcmp("StressZ", name.c_str())) return TPZElasticity3DGD::EStressZ;
  if (!strcmp("I1Strain", name.c_str())) return TPZElasticity3DGD::EI1Strain;
  if (!strcmp("I1", name.c_str())) return TPZElasticity3DGD::EI1;
  if (!strcmp("I2", name.c_str())) return TPZElasticity3DGD::EI2;
  if (!strcmp("I3", name.c_str())) return TPZElasticity3DGD::EI3;
  if (!strcmp("StressTensor", name.c_str())) return TPZElasticity3DGD::EStressTensor;
  if (!strcmp("StrainTensor", name.c_str())) return TPZElasticity3DGD::EStrainTensor;
  if (!strcmp("DeviatoricStressTensor", name.c_str())) return TPZElasticity3DGD::EDeviatoricStressTensor;
  if (!strcmp("DeviatoricStrainTensor", name.c_str())) return TPZElasticity3DGD::EDeviatoricStrainTensor;
  if (!strcmp("MatId", name.c_str())) return TPZElasticity3DGD::EMatId;
  if (!strcmp("StressFlac", name.c_str())) return TPZElasticity3DGD::EStressFlac;
  if (!strcmp("yielding", name.c_str())) return TPZElasticity3DGD::EYielding;
  if (!strcmp("Damage", name.c_str())) return TPZElasticity3DGD::EDamage;
  if (!strcmp("StrainZ", name.c_str())) return TPZElasticity3DGD::EStrainZ;
  if (!strcmp("PlasticStrainZ", name.c_str())) return TPZElasticity3DGD::EPlasticStrainZ;
  if (!strcmp("PlasticFunctionVal", name.c_str())) return TPZElasticity3DGD::EPlasticFunctionVal;
  if (!strcmp("ElementSize", name.c_str())) return TPZElasticity3DGD::EElementSize;
  if (!strcmp("GeoElIndex", name.c_str())) return TPZElasticity3DGD::EGeoElIndex;

  //cmesh->ElementSolution()
  if (!strcmp("AvDamage", name.c_str())) return 100;
  if (!strcmp("AvStressX", name.c_str())) return 101;
  if (!strcmp("AvStressY", name.c_str())) return 102;
  if (!strcmp("AvStressZ", name.c_str())) return 103;

  // cout << "TPZElasticityMaterial::VariableIndex Error\n";
  return TPZDiscontinuousGalerkin::VariableIndex(name);
}

int TPZElasticity3DGD::NSolutionVariables(int var) {
  switch (var) {
  case TPZElasticity3DGD::EStressTensor:
  case TPZElasticity3DGD::EDeviatoricStressTensor:
  case TPZElasticity3DGD::EStrainTensor:
  case TPZElasticity3DGD::EDeviatoricStrainTensor:
  case TPZElasticity3DGD::EStressFlac:
    return 6;
  case TPZElasticity3DGD::EDisplacement:
  case TPZElasticity3DGD::EPrincipalStress:
  case TPZElasticity3DGD::EPrincipalStrain:
  case TPZElasticity3DGD::EPrincipalPlasticStrain:
  case TPZElasticity3DGD::EPrincipalDirection1:
  case TPZElasticity3DGD::EPrincipalDirection2:
  case TPZElasticity3DGD::EPrincipalDirection3:
  case TPZElasticity3DGD::EStress:
  case TPZElasticity3DGD::EStrain:
  case TPZElasticity3DGD::ENormalStress:
  case TPZElasticity3DGD::ENormalStrain:
    return 3;
  case TPZElasticity3DGD::EDisplacementX:
  case TPZElasticity3DGD::EDisplacementY:
  case TPZElasticity3DGD::EDisplacementZ:
  case TPZElasticity3DGD::EVonMisesStress:
  case TPZElasticity3DGD::EStrain1:
  case TPZElasticity3DGD::EStress1:
  case TPZElasticity3DGD::EStressX:
  case TPZElasticity3DGD::EStressY:
  case TPZElasticity3DGD::EStressZ:
  case TPZElasticity3DGD::EI1Strain:
  case TPZElasticity3DGD::EI1:
  case TPZElasticity3DGD::EI2:
  case TPZElasticity3DGD::EI3:
  case TPZElasticity3DGD::EMatId:
  case TPZElasticity3DGD::EYielding:
  case TPZElasticity3DGD::EDamage:
  case TPZElasticity3DGD::EStrainZ:
  case TPZElasticity3DGD::EPlasticStrainZ:
  case TPZElasticity3DGD::EPlasticFunctionVal:
  case TPZElasticity3DGD::EElementSize:
  case TPZElasticity3DGD::EGeoElIndex:
    return 1;
  case TPZElasticity3DGD::EStressVector: return 6;
  default:
    return TPZMaterial::NSolutionVariables(var);
  }
  return -1;
}

void TPZElasticity3DGD::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout) {

  if(var == TPZElasticity3DGD::EElementSize){
    Solout[0] = data.HSize;
    return;
  }

  if(var == TPZElasticity3DGD::EGeoElIndex){
 //   Solout[0] = data.fGeoElIndex;
	  Solout[0] = data.gelElId;
	  return;
  }

  TPZVec<REAL> &Sol = data.sol[0];
  TPZFMatrix<STATE> &DSol = data.dsol[0];
  TPZFMatrix<STATE> &axes = data.axes;
  bool yielding;

  //o pos processamento, normalmente nao eh no ponto de integracao. E agora ?
  //Vou extrapolar o valor do primeiro ponto de integracao pra todo o elemento.
  if(data.intLocPtIndex < 0){
    data.intLocPtIndex = 0;
  }

  if(var == TPZElasticity3DGD::EDamage){
//    TPZManVector<REAL,6> Strain(6);
//    this->ComputeStrainVector(Strain, DSol);
//    Solout[0] = this->GetDamage(/*Strain,*/ data.gelElId, data.intLocPtIndex);
    Solout[0] = 0.0;
    return;
  }

  if(var == TPZElasticity3DGD::EStrainZ){
    TPZManVector<REAL,6> Strain(6);
    this->ComputeStrainVector(Strain, DSol);
    Solout[0] = Strain[2];
    return;
  }

  if(var == TPZElasticity3DGD::EPlasticStrainZ){

//    if(this->fPlasticModel){
//      if(this->fPlasticModel->fLastPlasticDeformation.IsInitialized()  ){
//        TPZFNMatrix<9> plasticStrain;
//        this->fPlasticModel->fLastPlasticDeformation.Get(data.gelElId, data.intLocPtIndex, plasticStrain);
//        Solout[0] = plasticStrain(2,2);
//        return;
//      }
//    }
    Solout[0] = 0.;
    return;
  }

  if(var == TPZElasticity3DGD::EPrincipalPlasticStrain){

  //  if(this->fPlasticModel){
  //    if(this->fPlasticModel->fLastPlasticDeformation.IsInitialized()  ){
  //      TPZFNMatrix<9> plasticStrain(3,3);
  //      this->fPlasticModel->fLastPlasticDeformation.Get(data.gelElId, data.intLocPtIndex, plasticStrain);
  //      int64_t numiterations = 1000;
  //      REAL tol = TPZElasticity3DGD::gTolerance;
  //      TPZManVector<REAL, 3> eigv(3);
  //      bool result = plasticStrain.SolveEigenvaluesJacobi(numiterations, tol, &eigv);
  //      for (int i = 0; i < eigv.size(); i++) {
  //        Solout[i] = eigv[i];
  //      }
  //  #ifdef DEBUG
  //      if (result == false) {
  //        PZError << __PRETTY_FUNCTION__ <<
  //            " - ERROR! - result = false - numiterations = " << numiterations <<
  //            " - tol = " << tol << std::endl;
  //      }
  //  #endif
  //      return;
  //    }
  //  }
    Solout[0] = 0.;
    return;
  }

  if(var == TPZElasticity3DGD::EMatId){
    Solout[0] = this->Id();
    return;
  }
  if (var == TPZElasticity3DGD::EDisplacement) {
    int i;
    for (i = 0; i < 3; i++) {
      Solout[i] = Sol[i];
    } // for
    return;
  } // TPZElasticity3DGD::EDisplacement

  if (var == TPZElasticity3DGD::EDisplacementX) {
    Solout[0] = Sol[0];
    return;
  } // TPZElasticity3DGD::EDisplacementX

  if (var == TPZElasticity3DGD::EDisplacementY) {
    Solout[0] = Sol[1];
    return;
  } // TPZElasticity3DGD::EDisplacementY

  if (var == TPZElasticity3DGD::EDisplacementZ) {
    Solout[0] = Sol[2];
    return;
  } // TPZElasticity3DGD::EDisplacementZ

  if(var == TPZElasticity3DGD::EPlasticFunctionVal){
    REAL PHI = +1.;
  //  if(this->fPlasticModel){
  //    TPZFNMatrix<9> StressTensor(3, 3);
  //    REAL damage = 0.;
  //    this->ComputeStressTensor(StressTensor, yielding, damage, DSol, data.gelElId, data.intLocPtIndex, data.HSize);
  //    bool res = this->fPlasticModel->CheckPlasticFunction(StressTensor,damage,PHI);
  //  }
    Solout[0] = PHI;
    return;
  }//TPZElasticity3DGD::EPlasticFunctionVal

  if (var == TPZElasticity3DGD::EPrincipalStress) {
    TPZFNMatrix<9> StressTensor(3, 3);
    REAL damage = 0.;
    this->ComputeStressTensor(StressTensor, yielding, damage, DSol, data.gelElId, data.intLocPtIndex, data.HSize);
    int64_t numiterations = 1000;
    REAL tol = TPZElasticity3DGD::gTolerance;
    bool result = StressTensor.SolveEigenvaluesJacobi(numiterations, tol, &Solout);
#ifdef DEBUG
    if (result == false) {
      PZError << __PRETTY_FUNCTION__ <<
          " - ERROR! - result = false - numiterations = " << numiterations <<
          " - tol = " << tol << std::endl;
    }
#endif
    return;
  } // TPZElasticity3DGD::EPrincipalStress

  if (var == TPZElasticity3DGD::EStress1) {
    TPZFNMatrix<9>StressTensor(3, 3);
    TPZManVector<REAL, 3>PrincipalStress(3);
    REAL damage = 0.;
    this->ComputeStressTensor(StressTensor,yielding, damage, DSol, data.gelElId, data.intLocPtIndex, data.HSize);
    int64_t numiterations = 1000;
    REAL tol = TPZElasticity3DGD::gTolerance;
    bool result = StressTensor.SolveEigenvaluesJacobi(numiterations, tol, &PrincipalStress);
    Solout[0] = PrincipalStress[0];
#ifdef DEBUG
    if (result == false) {
      PZError << __PRETTY_FUNCTION__ <<
          " - ERROR! - result = false - numiterations = " << numiterations <<
          " - tol = " << tol << std::endl;
    }
#endif
    return;
  } // TPZElasticity3DGD::EStress1

  if (var == TPZElasticity3DGD::EPrincipalStrain) {
    TPZFNMatrix<9>StrainTensor(3, 3);
    this->ComputeStrainTensor(StrainTensor, DSol);
    int64_t numiterations = 1000;
    REAL tol = TPZElasticity3DGD::gTolerance;
    TPZManVector<REAL, 3>eigv;
    bool result;
    result = StrainTensor.SolveEigenvaluesJacobi(numiterations, tol, &eigv);
    for (int i = 0; i < eigv.size(); i++) {
      Solout[i] = eigv[i];
    }
#ifdef DEBUG
    if (result == false) {
      PZError << __PRETTY_FUNCTION__ <<
          " - ERROR! - result = false - numiterations = " << numiterations <<
          " - tol = " << tol << std::endl;
    }
#endif
    return;
  } // TPZElasticity3DGD::EPrincipalStrain

  if (var == TPZElasticity3DGD::EStrain1) {
    TPZFNMatrix<9>StrainTensor(3, 3);
    TPZManVector<REAL, 3>PrincipalStrain(3);
    this->ComputeStrainTensor(StrainTensor, DSol);
    int64_t numiterations = 1000;
    REAL tol = TPZElasticity3DGD::gTolerance;
    bool result;
    result = StrainTensor.SolveEigenvaluesJacobi(numiterations, tol, &PrincipalStrain);
    Solout[0] = PrincipalStrain[0];
#ifdef DEBUG
    if (result == false) {
      PZError << __PRETTY_FUNCTION__ <<
          " - ERROR! - result = false - numiterations = " << numiterations <<
          " - tol = " << tol << std::endl;
    }
#endif
    return;
  } // TPZElasticity3DGD::EPrincipalStrain

  TPZManVector<REAL, 3>eigvec;
  if (var == TPZElasticity3DGD::EPrincipalDirection1) {
    this->PrincipalDirection(DSol, eigvec, 0);
    for (int i = 0; i < eigvec.size(); i++) {
      Solout[i] = eigvec[i];
    }
    return;
  } // TPZElasticity3DGD::EPrincipalDirection1

  if (var == TPZElasticity3DGD::EPrincipalDirection2) {
    this->PrincipalDirection(DSol, eigvec, 1);
    for (int i = 0; i < eigvec.size(); i++) {
      Solout[i] = eigvec[i];
    }
    return;
  } // TPZElasticity3DGD::EPrincipalDirection2

  if (var == TPZElasticity3DGD::EPrincipalDirection3) {
    this->PrincipalDirection(DSol, eigvec, 2);
    for (int i = 0; i < eigvec.size(); i++) {
      Solout[i] = eigvec[i];
    }
    return;
  } // TPZElasticity3DGD::EPrincipalDirection3

  if (var == TPZElasticity3DGD::EVonMisesStress) {

    DebugStop();//nossa, que coisa eh esta ?
    TPZManVector<REAL, 3>PrincipalStress(3);
    TPZFNMatrix<9>StressTensor(3, 3);
    REAL damage = 0.;
    this->ComputeStressTensor(StressTensor, yielding, damage, DSol, data.gelElId, data.intLocPtIndex, data.HSize);
    int64_t numiterations = 1000;
    REAL tol = TPZElasticity3DGD::gTolerance;
    bool result;
    result = StressTensor.SolveEigenvaluesJacobi(numiterations, tol, &PrincipalStress);
#ifdef DEBUG
    if (result == false) {
      PZError << __PRETTY_FUNCTION__ <<
          " - ERROR! - result = false - numiterations = " << numiterations <<
          " - tol = " << tol << std::endl;
    }
#endif

    Solout.Resize(1);
    Solout[0] = (PrincipalStress[0] - PrincipalStress[1]) *
                (PrincipalStress[0] - PrincipalStress[1]) +
                (PrincipalStress[1] - PrincipalStress[2]) *
                (PrincipalStress[1] - PrincipalStress[2]) +
                (PrincipalStress[2] - PrincipalStress[0]) *
                (PrincipalStress[2] - PrincipalStress[0]);
    Solout[0] = Solout[0] / (2. * this->fFy*this->fFy);
    return;
  } // TPZElasticity3DGD::EVonMisesStress

  if (var == TPZElasticity3DGD::EStress) {
    TPZManVector<REAL,6> Stress(6);
    REAL damage = 0.;
    this->ComputeStressVector(Stress,yielding, damage,DSol, data.gelElId, data.intLocPtIndex, data.HSize);
    this->ApplyDirection(Stress, Solout);
    return;
  } // TPZElasticity3DGD::EStress

  if (var ==  TPZElasticity3DGD::EYielding) {
    TPZManVector<REAL,6> Stress(6);
    REAL damage = 0.;
    this->ComputeStressVector(Stress,yielding, damage,DSol, data.gelElId, data.intLocPtIndex, data.HSize);
    if( yielding ) Solout[0] = 1.;
    else Solout[0] = 0.;
    return;
  } //  TPZElasticity3DGD::EYielding


  if (var == TPZElasticity3DGD::EStressVector){
    REAL damage = 0.;
    this->ComputeStressVector(Solout, yielding,damage,DSol, data.gelElId, data.intLocPtIndex, data.HSize);
    return;
  } // TPZElasticity3DGD::EStressVector

  if (var == TPZElasticity3DGD::EStrain) {
    TPZManVector<REAL,6>Strain(6);
    this->ComputeStrainVector(Strain, DSol);
    this->ApplyDirection(Strain, Solout);
    return;
  } // TPZElasticity3DGD::EStrain

  if (var == TPZElasticity3DGD::ENormalStress) {
    TPZManVector<REAL,6> Stress(6);
    REAL damage = 0.;
    this->ComputeStressVector(Stress, yielding,damage,DSol, data.gelElId, data.intLocPtIndex, data.HSize);
    Solout[0] = Stress[0];
    Solout[1] = Stress[1];
    Solout[2] = Stress[2];
    return;
  } // TPZElasticity3DGD::ENormalStress

  if (var == TPZElasticity3DGD::ENormalStrain) {
    TPZManVector<REAL,6> Strain(6);
    this->ComputeStrainVector(Strain, DSol);
    Solout[0] = Strain[0];
    Solout[1] = Strain[1];
    Solout[2] = Strain[2];
    return;
  } // TPZElasticity3DGD::ENormalStrain

  REAL damage = 0.;
  if (var == TPZElasticity3DGD::EStressX) {
    TPZFNMatrix<9>Stress(3, 3);
    this->ComputeStressTensor(Stress, yielding, damage,DSol,data.gelElId, data.intLocPtIndex, data.HSize);
    Solout[0] = Stress(0, 0);
    return;
  } // TPZElasticity3DGD::EStressX

  if (var == TPZElasticity3DGD::EStressY) {
    TPZFNMatrix<9>Stress(3, 3);
    this->ComputeStressTensor(Stress, yielding, damage,DSol,data.gelElId, data.intLocPtIndex, data.HSize);
    Solout[0] = Stress(1, 1);
    return;
  } // TPZElasticity3DGD::EStressY

  if (var == TPZElasticity3DGD::EStressZ) {
    TPZFNMatrix<9>Stress(3, 3);
    this->ComputeStressTensor(Stress, yielding, damage,DSol,data.gelElId, data.intLocPtIndex, data.HSize);
    Solout[0] = Stress(2, 2);
    return;
  } // TPZElasticity3DGD::EStressZ

  if (var == TPZElasticity3DGD::EI1Strain) {
    TPZManVector<REAL,6> Strain(6);
    this->ComputeStrainVector(Strain, DSol);
    Solout[0] = Strain[0] + Strain[1] + Strain[2];
    return;
  } // I1Strain

  if (var == TPZElasticity3DGD::EI1) {
    TPZFNMatrix<9>Stress(3, 3);
    this->ComputeStressTensor(Stress, yielding, damage,DSol,data.gelElId, data.intLocPtIndex, data.HSize);
    Solout[0] = Stress(0, 0) + Stress(1, 1) + Stress(2, 2);
    // Calculado no Mathematica

    return;
  } // I1

  if (var == TPZElasticity3DGD::EI2) {
    TPZFNMatrix<9>Stress(3, 3);
    this->ComputeStressTensor(Stress, yielding, damage,DSol,data.gelElId, data.intLocPtIndex, data.HSize);

    Solout[0] = -(Stress(0, 1) * Stress(1, 0)) - Stress(0, 2) * Stress(2, 0) -
                  Stress(1, 2) * Stress(2, 1)  + Stress(1, 1) * Stress(2, 2) +
                  Stress(0, 0) * (Stress(1, 1) + Stress(2, 2));
    // Calculado no Mathematica

    return;
  } // I2

  if (var == TPZElasticity3DGD::EI3) {
    TPZFNMatrix<9>Stress(3, 3);
    this->ComputeStressTensor(Stress, yielding, damage,DSol,data.gelElId, data.intLocPtIndex, data.HSize);
    Solout[0] = - (Stress(0, 2) * Stress(1, 1) * Stress(2, 0))
                + (Stress(0, 1) * Stress(1, 2) * Stress(2, 0))
                + (Stress(0, 2) * Stress(1, 0) * Stress(2, 1))
                - (Stress(0, 0) * Stress(1, 2) * Stress(2, 1))
                - (Stress(0, 1) * Stress(1, 0) * Stress(2, 2))
                + (Stress(0, 0) * Stress(1, 1) * Stress(2, 2)); // Calculado no Mathematica

    return;
  } // I3

  if (var == TPZElasticity3DGD::EStressTensor) {
    Solout.Resize(6,0.);
    TPZFNMatrix<9>Stress(3,3);
    this->ComputeStressTensor(Stress,yielding,damage,DSol,data.gelElId, data.intLocPtIndex, data.HSize);
    Solout[0] = Stress(0,0);
    Solout[1] = Stress(1,1);
    Solout[2] = Stress(2,2);
    Solout[3] = Stress(0,1);
    Solout[4] = Stress(0,2);
    Solout[5] = Stress(1,2);
    return;
  }

  if (var == TPZElasticity3DGD::EStressFlac) {

    TPZFNMatrix<9>Stress(3,3);
    this->ComputeStressTensor(Stress,yielding,damage,DSol,data.gelElId, data.intLocPtIndex, data.HSize);
    Solout[0] = Stress(0,0);
    Solout[1] = Stress(1,0);
    Solout[2] = Stress(2,0);
    Solout[3] = Stress(3,0);
    Solout[4] = Stress(4,0);
    Solout[5] = Stress(5,0);
    return;
  }


  if (var == TPZElasticity3DGD::EStrainTensor) {
    Solout.Resize(6,0.);
    TPZFNMatrix<6>Strain(6, 1);
    this->ComputeStrainTensor(Strain, DSol);
    Solout[0] = Strain(0,0);
    Solout[1] = Strain(1,0);
    Solout[2] = Strain(2,0);
    Solout[3] = Strain(3,0);
    Solout[4] = Strain(4,0);
    Solout[5] = Strain(5,0);
    return;
  }


  if (var == TPZElasticity3DGD::EDeviatoricStrainTensor) {
    Solout.Resize(6,0.);
    TPZFNMatrix<6>Strain(6, 1);
    this->ComputeStrainTensor(Strain, DSol);
    const REAL I1s3 = (Strain(0,0) + Strain(1,0) + Strain(2,0))/3.;
    Solout[0] = Strain(0,0) - I1s3;
    Solout[1] = Strain(1,0) - I1s3;
    Solout[2] = Strain(2,0) - I1s3;
    Solout[3] = Strain(3,0);
    Solout[4] = Strain(4,0);
    Solout[5] = Strain(5,0);
    return;
  }


  if (var == TPZElasticity3DGD::EDeviatoricStressTensor) {
    Solout.Resize(6,0.);
    TPZFNMatrix<9>Stress(3,3);
    this->ComputeStressTensor(Stress,yielding,damage,DSol,data.gelElId, data.intLocPtIndex, data.HSize);
    const REAL I1s3 = (Stress(0, 0) + Stress(1, 1) + Stress(2, 2))/3.;
    Solout[0] = Stress(0,0) - I1s3;
    Solout[1] = Stress(1,1) - I1s3;
    Solout[2] = Stress(2,2) - I1s3;
    Solout[3] = Stress(0,1);
    Solout[4] = Stress(0,2);
    Solout[5] = Stress(1,2);
    return;
  }

  TPZMaterial::Solution(data, var, Solout);

} // Solution

void TPZElasticity3DGD::Errors(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix<STATE> &dudx,
	TPZFMatrix<STATE> &axes, TPZVec<REAL> &flux, TPZVec<REAL> &u_exact,
	TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &values) {
  int i, j;

  /** L2 norm */
  REAL L2 = 0.;
  for (i = 0; i < 3; i++)
    L2 += (u[i] - u_exact[i]) * (u[i] - u_exact[i]);

  /** H1 semi-norm */
  REAL SemiH1 = 0.;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      SemiH1 += (dudx(i, j) - du_exact(i, j)) * (dudx(i, j) - du_exact(i, j));

  /** H1 norm */
  REAL H1 = L2 + SemiH1;

  // values[1] : eror em norma L2
  values[1] = L2;

  // values[2] : erro em semi norma H1
  values[2] = SemiH1;

  // values[0] : erro em norma H1 <=> norma Energia
  values[0] = H1;

}

void TPZElasticity3DGD::ComputeStrainTensor(TPZFMatrix<STATE> &Strain, const TPZFMatrix<STATE> &DSol) const {
  Strain.Redim(3, 3);
  Strain(0, 0) = DSol.GetVal(0, 0);
  Strain(0, 1) = 0.5 * (DSol.GetVal(1, 0) + DSol.GetVal(0, 1));
  Strain(0, 2) = 0.5 * (DSol.GetVal(2, 0) + DSol.GetVal(0, 2));
  Strain(1, 0) = Strain(0, 1);
  Strain(1, 1) = DSol.GetVal(1, 1);
  Strain(1, 2) = 0.5 * (DSol.GetVal(2, 1) + DSol.GetVal(1, 2));
  Strain(2, 0) = Strain(0, 2);
  Strain(2, 1) = Strain(1, 2);
  Strain(2, 2) = DSol.GetVal(2, 2);
}

// REAL TPZElasticity3DGD::GetDamage( /*const TPZVec<REAL> &Strain, */ int elindex, int elIntPoint) const{
/*  REAL damage = 0.0;
  if(this->fDamageModel){
//    TPZFNMatrix<9> StrainTensor(3,3);//nao precisa mais do strain
//    this->Vector2Tensor(Strain, StrainTensor);
    this->fDamageModel->GetDamage(elindex, elIntPoint, damage);
  }
  return damage;
}*/

void TPZElasticity3DGD::ComputeStressTensor(TPZFMatrix<STATE> &Stress, bool &yielding, REAL &damage, const TPZFMatrix<STATE> &dsol,
int GeoElIndex, int IntPtIndex, const REAL elementCharacteristicSize) const {
  TPZManVector<REAL,6> Vec(6);
  this->ComputeStressVector(Vec, yielding, damage,dsol, GeoElIndex, IntPtIndex, elementCharacteristicSize);
  this->Vector2Tensor(Vec,Stress);
}

void TPZElasticity3DGD::ComputeStrainVector(TPZVec<REAL> &Strain, const TPZFMatrix<STATE> &DSol) const {
  Strain.Resize(6);
  Strain[0] = DSol.GetVal(0, 0);
  Strain[1] = DSol.GetVal(1, 1);
  Strain[2] = DSol.GetVal(2, 2);
  Strain[3] = 0.5 * (DSol.GetVal(1, 0) + DSol.GetVal(0, 1));
  Strain[4] = 0.5 * (DSol.GetVal(2, 0) + DSol.GetVal(0, 2));
  Strain[5] = 0.5 * (DSol.GetVal(2, 1) + DSol.GetVal(1, 2));
}

void TPZElasticity3DGD::ComputeStressVector(TPZVec<REAL> &Stress,bool &yielding, REAL &damage,
	const TPZFMatrix<STATE> &DSol, int elindex, int elIntPoint, const REAL elementCharacteristicSize) const {
  TPZManVector<REAL,6> StrainVec(6);
  this->ComputeStrainVector(StrainVec, DSol);
  this->ComputeStressVector(Stress,yielding,damage,StrainVec,elindex,elIntPoint,elementCharacteristicSize);
}

void TPZElasticity3DGD::ComputeStressVector(TPZVec<REAL> &StressVec, bool &yielding,
                                          REAL &damage,
                                          const TPZVec<REAL> &Strain,
                                          int elindex, int elIntPoint,
                                          const REAL elementCharacteristicSize) const {
  TPZFNMatrix<9> StrainTensor(3,3);
  this->Vector2Tensor(Strain, StrainTensor);

  TPZManVector<REAL,6> elasticStrain = Strain;
  TPZFNMatrix<9> plasticStrain(3,3,0.);
//  if(this->fPlasticModel){
//    if(this->fPlasticModel->fLastPlasticDeformation.IsInitialized()  ){
//      this->fPlasticModel->fLastPlasticDeformation.Get(elindex, elIntPoint, plasticStrain);
//      elasticStrain[0] -= plasticStrain.GetVal(0, 0);
//      elasticStrain[1] -= plasticStrain.GetVal(1, 1);
//      elasticStrain[2] -= plasticStrain.GetVal(2, 2);
//      elasticStrain[3] -= plasticStrain.GetVal(1, 0);
//      elasticStrain[4] -= plasticStrain.GetVal(0, 2);
//      elasticStrain[5] -= plasticStrain.GetVal(1, 2);
//    }
//    else{
      //primeiro passo de tempo, logo nao ha historico ainda
//    }
//  }

//  damage = this->GetDamage(/*Strain,*/ elindex, elIntPoint);
  damage = 0.0;
  const REAL E = this->fE*(1.-damage);
  const REAL ni = this->fPoisson;
  const REAL const1 = (1. + ni)*(1 - 2.*ni);
  const REAL c1 = E/const1;
  const REAL c2 = 1. - ni;
  const REAL c3 = 1. - 2.*ni;

  StressVec.Resize(6);
  StressVec[0] = c1 * (c2*elasticStrain[0] + ni * (elasticStrain[1] + elasticStrain[2])) + fPreStress[0];
  StressVec[1] = c1 * (c2*elasticStrain[1] + ni * (elasticStrain[0] + elasticStrain[2])) + fPreStress[1];
  StressVec[2] = c1 * (c2*elasticStrain[2] + ni * (elasticStrain[0] + elasticStrain[1])) + fPreStress[2];

  StressVec[3] = c1 * (c3 * elasticStrain[3]);
  StressVec[4] = c1 * (c3 * elasticStrain[4]);
  StressVec[5] = c1 * (c3 * elasticStrain[5]);

/*  if(this->fPlasticModel){
    TPZFNMatrix<9> StressTrial;
    this->Vector2Tensor(StressVec,StressTrial);

    TPZFNMatrix<9> PlasticModelStress, DeltaPlasticStrain;
    yielding = this->fPlasticModel->ApplyPlasticModel(StressTrial,damage,PlasticModelStress,DeltaPlasticStrain);
    this->fPlasticModel->fCurrentDeltaPlasticDeformation.Add(elindex, elIntPoint, DeltaPlasticStrain);
    this->Tensor2Vector(PlasticModelStress,StressVec);
    if(this->fDamageModel){
      //new plastic strain
      TPZFNMatrix<9> updatedPlasticStrain = plasticStrain;
      updatedPlasticStrain += DeltaPlasticStrain;

      //new elastic strain
      TPZFNMatrix<9> updatedElasticStrain(3,3);
      updatedElasticStrain = StrainTensor;
      updatedElasticStrain -= updatedPlasticStrain;

      this->fDamageModel->UpdateCurrentDamage(updatedElasticStrain, updatedPlasticStrain, elindex, elIntPoint, elementCharacteristicSize);
    }
  }
  else{//elastico apenas
    yielding = false;
    if(this->fDamageModel){
      plasticStrain.Redim(3,3);
      this->fDamageModel->UpdateCurrentDamage(StrainTensor,plasticStrain,elindex, elIntPoint, elementCharacteristicSize);
    }
  }*/
}

void TPZElasticity3DGD::ApplyDirection(const TPZVec<REAL> &StrVec, TPZVec<REAL> &Out) const {
  Out.Resize(3);
  const TPZVec<REAL>&Dir = this->fPostProcessDirection;
  Out[0] = Dir[0] * StrVec[0] + Dir[1] * StrVec[3] + Dir[2] * StrVec[4];
  Out[1] = Dir[1] * StrVec[1] + Dir[0] * StrVec[3] + Dir[2] * StrVec[5];
  Out[2] = Dir[2] * StrVec[2] + Dir[0] * StrVec[4] + Dir[1] * StrVec[5];
}

void TPZElasticity3DGD::PrincipalDirection(const TPZFMatrix<STATE> &DSol, TPZVec<REAL> &Solout, int direction) const {
  TPZFNMatrix<9>StrainTensor(3, 3);
  TPZManVector<REAL, 3>Eigenvalues;
  TPZFNMatrix<9>Eigenvectors(3, 3);

  this->ComputeStrainTensor(StrainTensor, DSol);
  int64_t numiterations = 1000;
  REAL tol = TPZElasticity3DGD::gTolerance;
  bool result;
  result = StrainTensor.SolveEigensystemJacobi(numiterations, tol, Solout, Eigenvectors);
  // Solout is used to store Eigenvaleus, but its values will be replaced below
#ifdef DEBUG
  if (result == false) {
    PZError << __PRETTY_FUNCTION__ <<
        " - ERROR! - result = false - numiterations = " << numiterations <<
        " - tol = " << tol << std::endl;
  }
#endif
  Solout.Resize(3);
  for (int i = 0; i < 3; i++) {
    Solout[i] = Eigenvectors(direction, i);
  }
}

/** Save the element data to a stream */
void TPZElasticity3DGD::Write(TPZStream &buf, int withclassid) const{
  TPZMaterial::Write(buf, withclassid);
  buf.Write(&fE, 1);
  buf.Write(&fForce[0], 3);
  buf.Write(&fFy, 1);
  buf.Write(&fPoisson, 1);
  buf.Write(&fPostProcessDirection[0], 3);
  DebugStop();//implementar fPlasticModel
}

/** Read the element data from a stream */
void TPZElasticity3DGD::Read(TPZStream &buf, void *context) {
  TPZMaterial::Read(buf, context);
  buf.Read(&fE, 1);
  fForce.Resize(3);
  buf.Read(&fForce[0], 3);
  buf.Read(&fFy, 1);
  buf.Read(&fPoisson, 1);
  fPostProcessDirection.Resize(3);
  buf.Read(&fPostProcessDirection[0], 3);
  DebugStop();//implementar fPlasticModel
}

int TPZElasticity3DGD::ClassId() const {
  return TPZELASTICITY3DMATERIALID;
}

void TPZElasticity3DGD::FillDataRequirements(TPZMaterialData &data) {
  TPZMaterial::FillDataRequirements(data);
  data.fNeedsSol = true;
  data.fNeedsHSize = true;
}

void TPZElasticity3DGD::SetIntegrationRule(TPZAutoPointer<TPZIntPoints> rule,
                                     int elPMaxOrder,
                                     int elDimension){

  TPZManVector<int,3> p2(elDimension,2*elPMaxOrder);
  rule->SetOrder(p2);
  if(this->HasForcingFunction()) {
    TPZManVector<int,3> order(elDimension,ForcingFunction()->PolynomialOrder());
    rule->SetOrder(order);
  }
}

void TPZElasticity3DGD::SetIntegrationRule(TPZAutoPointer<TPZIntPoints> rule,
                                      int elPMaxOrder,
                                      int elDimension,
                                      TPZBndCond &bc){
  int maxOrder = 0;
  if(bc.HasForcingFunction()){
    maxOrder = bc.ForcingFunction()->PolynomialOrder();
  }
  if(elPMaxOrder > maxOrder) maxOrder = elPMaxOrder;
  int integOrder = 2*maxOrder;
	TPZManVector<int,3> p2(elDimension,integOrder);
	rule->SetOrder(p2);
}

void TPZElasticity3DGD::Vector2Tensor(const TPZVec<REAL> &Vec, TPZFMatrix<STATE> &Tensor) const{
  Tensor.Redim(3, 3);
  Tensor(0, 0) = Vec[0];
  Tensor(0, 1) = Vec[3];
  Tensor(0, 2) = Vec[4];
  Tensor(1, 0) = Vec[3];
  Tensor(1, 1) = Vec[1];
  Tensor(1, 2) = Vec[5];
  Tensor(2, 0) = Vec[4];
  Tensor(2, 1) = Vec[5];
  Tensor(2, 2) = Vec[2];
}

void TPZElasticity3DGD::Tensor2Vector(const TPZFMatrix<STATE> &Tensor, TPZVec<REAL> &Vec) const{
  Vec.Resize(6);
  Vec[0] = Tensor.GetVal(0, 0);
  Vec[1] = Tensor.GetVal(1, 1);
  Vec[2] = Tensor.GetVal(2, 2);
  Vec[3] = Tensor.GetVal(1, 0);
  Vec[4] = Tensor.GetVal(0, 2);
  Vec[5] = Tensor.GetVal(1, 2);
}

void TPZElasticity3DGD::Contribute(TPZMaterialData &data, REAL weight,
	TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {

#ifdef _AUTODIFF
  REAL alpha = 1.;//1. significa apenas matriz jacobiana elastica
  const REAL eps = 1e-3;
  if( fabs(alpha - 0. ) > eps ){
    this->ContributeOriginal(data,alpha*weight,ek,ef);
  }
  if( fabs(alpha -1. ) > eps ){
    this->ContributeDifFinita(data,(1.-alpha)*weight,ek,ef);
  }
#else
  this->ContributeOriginal(data,weight,ek,ef);
#endif

}

void TPZElasticity3DGD::FillDataRequirementsInterface(TPZMaterialData &facedata){
  facedata.SetAllRequirements(true);
  facedata.fNeedsSol = true;
	facedata.fNeedsHSize = true;
}

void TPZElasticity3DGD::InterfaceGamma(const int elindex, const int elIntPoint,
                                     REAL h, REAL normalOpening,
                                     REAL &gamma, REAL &maxTensionStress, bool &limitedTensionStress) const{
  const REAL alpha = 1e-2;//tamarindo 1e-5;
  const REAL gammaVirginState = fE/(alpha*2.*h);

 /* if(this->fCohesiveDamageModel){
    gamma = 0.;
    this->fCohesiveDamageModel->GetDamageSlope(elindex, elIntPoint, normalOpening, gammaVirginState, gamma, maxTensionStress);
    limitedTensionStress = true;
    this->fCohesiveDamageModel->UpdateCurrentDamage(normalOpening, elindex, elIntPoint);
  }
  else{ */
    limitedTensionStress = false;
    maxTensionStress = 0.;
    gamma = gammaVirginState;
 // }
}

void TPZElasticity3DGD::ContributeInterface(TPZMaterialData &facedata,
                                   TPZMaterialData &leftdata, TPZMaterialData &rightdata,
								   REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){

	TPZFMatrix<STATE> &phiL = leftdata.phi;
	TPZFMatrix<STATE> &phiR = rightdata.phi;
  TPZManVector<REAL,3> &normal = facedata.normal;
  const int nphiL = phiL.Rows();
  const int nphiR = phiR.Rows();

  //eu posso usar a maxima abertura ou a abertura normal. Assim, ao inves de descontar a cisalhante (abertura normal), eu a incluiria
  REAL normalOpening = 0.;
  for(int i = 0; i < 3; i++){
    normalOpening += (rightdata.sol[0][i]-leftdata.sol[0][i]) * normal[i];
  }
// abertura total, mas prefiro a normal acima. Imagine uma face apenas cisalhando. Não quero dano de fratura modo I nela.
//  for(int i = 0; i < 3; i++){
//    normalOpening += pow(leftdata.sol[i]-rightdata.sol[i],2);
//  }
//  normalOpening = sqrt(normalOpening);


  const REAL h = (leftdata.HSize < rightdata.HSize) ? leftdata.HSize : rightdata.HSize;
  REAL gamma, maxTensionStress;
  bool limitedTensionStress;
  this->InterfaceGamma(facedata.gelElId, facedata.intLocPtIndex, h, normalOpening, gamma, maxTensionStress, limitedTensionStress );
  TPZManVector<REAL,3> stress(3,0.);
  for(int i = 0; i < 3; i++){
    stress[i] = gamma*(rightdata.sol[0][i]-leftdata.sol[0][i]);
  }
  bool tamarindo = false;
  if(limitedTensionStress && tamarindo){//para limitar a tensao a maxima admissivel ao atual estado de dano
    REAL normalStress = 0.;
    for(int i = 0; i < 3; i++) normalStress += stress[i]*normal[i];
    if(normalStress > 1e-6){//tem que ser positivo (tracao) e ser diferente de zero
      if(normalStress > maxTensionStress){
        const REAL factor = maxTensionStress/normalStress;
        for(int i = 0; i < 3; i++) stress[i] *= factor;
        gamma *= factor;
      }
    }
  }

  //left test functions {du, dv, dw} left
  for(int il = 0; il < nphiL; il++){

    for(int idim = 0; idim < 3; idim++){
      ef(il * 3 + idim, 0 ) += -1.*weight * (-stress[idim]) * phiL(il,0);
    }

    //left trial functions {u, v, w} left
    for(int jl = 0; jl < nphiL; jl++){
      //(dul,ul)
      ek(il * 3 + 0, jl * 3 + 0) += weight * gamma * phiL(il,0) * phiL(jl,0);
      //(dvl,vl)
      ek(il * 3 + 1, jl * 3 + 1) += weight * gamma * phiL(il,0) * phiL(jl,0);
      //(dwl,wl)
      ek(il * 3 + 2, jl * 3 + 2) += weight * gamma * phiL(il,0) * phiL(jl,0);
    }

    //right trial functions {u, v, w} right
    for(int jr = 0; jr < nphiR; jr++){
      //(dul,ur)
      ek(il * 3 + 0, 3*nphiL + jr * 3 + 0) += weight * (-gamma) * phiL(il,0) * phiR(jr,0);
      //(dvl,vr)
      ek(il * 3 + 1, 3*nphiL + jr * 3 + 1) += weight * (-gamma) * phiL(il,0) * phiR(jr,0);
      //(dwl,wr)
      ek(il * 3 + 2, 3*nphiL + jr * 3 + 2) += weight * (-gamma) * phiL(il,0) * phiR(jr,0);
    }
  }

  //right test functions {du, dv, dw} right
  for(int ir = 0; ir < nphiR; ir++){

    for(int idim = 0; idim < 3; idim++){
      ef(3*nphiL + ir * 3 + idim, 0 ) += -1.*weight * (+stress[idim]) * phiR(ir,0);
    }

    //left trial functions {u, v, w} left
    for(int jl = 0; jl < nphiL; jl++){
      //(dur,ul)
      ek(3*nphiL + ir * 3 + 0, jl * 3 + 0) += weight * (-gamma) * phiR(ir,0) * phiL(jl,0);
      //(dvr,vl)
      ek(3*nphiL + ir * 3 + 1, jl * 3 + 1) += weight * (-gamma) * phiR(ir,0) * phiL(jl,0);
      //(dwr,wl)
      ek(3*nphiL + ir * 3 + 2, jl * 3 + 2) += weight * (-gamma) * phiR(ir,0) * phiL(jl,0);
    }

    //right trial functions {u, v, w} right
    for(int jr = 0; jr < nphiR; jr++){
      //(dur,ur)
      ek(3*nphiL + ir * 3 + 0, 3*nphiL + jr * 3 + 0) += weight * gamma * phiR(ir,0) * phiR(jr,0);
      //(dvr,vr)
      ek(3*nphiL + ir * 3 + 1, 3*nphiL + jr * 3 + 1) += weight * gamma * phiR(ir,0) * phiR(jr,0);
      //(dwr,wr)
      ek(3*nphiL + ir * 3 + 2, 3*nphiL + jr * 3 + 2) += weight * gamma * phiR(ir,0) * phiR(jr,0);
    }
  }

#ifdef DEBUG
  if ( !ek.VerifySymmetry( 1.e-3 ) ){
    std::cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << std::endl;
  }
#endif

}///void

void TPZElasticity3DGD::ContributeBCInterface(TPZMaterialData & /*facedata*/ ,
                                     TPZMaterialData &leftdata,
									 REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                                     TPZBndCond &bc){
  //implementando BC de formulacao continua classica
  this->ContributeBC(leftdata,weight,ek,ef,bc);
}



#ifndef BORLAND

template class TPZRestoreClass<TPZElasticity3DGD, TPZELASTICITY3DMATERIALID>;
#endif

