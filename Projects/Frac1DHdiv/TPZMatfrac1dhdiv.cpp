#include "pzlog.h"
#include "TPZMatfrac1dhdiv.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"

#include <iostream>

#ifdef PZ_LOG
static PZLogger logger("pz.multiphase");
#endif

#ifdef PZ_LOG
static PZLogger logdata("pz.material.multiphase.data");
#endif

TPZMatfrac1dhdiv::TPZMatfrac1dhdiv(): TPZMatWithMem<TPZFMatrix<REAL>, TPZMaterial >()
{
  fDim = 1;
  TPZFNMatrix<3,REAL> Vl(1,1,0.);
  this->SetDefaultMem(Vl);
}

TPZMatfrac1dhdiv::TPZMatfrac1dhdiv(int matid): TPZMatWithMem<TPZFMatrix<REAL>, TPZMaterial >(matid)
{
  fDim = 1;
  TPZFNMatrix<3,REAL> Vl(1,1,0.);
  this->SetDefaultMem(Vl);
}

TPZMatfrac1dhdiv::~TPZMatfrac1dhdiv()
{
  
}

int TPZMatfrac1dhdiv::Dimension() const {return fDim;};


void TPZMatfrac1dhdiv::Print(std::ostream &out) {
  out << "name of material : " << Name() << "\n";
  out << "Coeficient which multiplies the gradient operator "<< "my var" << std::endl;
  out << "Base Class properties :";
  TPZMaterial::Print(out);
  out << "\n";
}

REAL TPZMatfrac1dhdiv::Getw(REAL p)
{
  DebugStop(); // deprecated
  return 0.817 * (1. - fData->Poisson()) * fData->Hf() / fData->G() * (p - fData->SigmaConf());
}

REAL TPZMatfrac1dhdiv::Getdwdp()
{
  DebugStop(); // deprecated
  return 0.817 * (1. - fData->Poisson()) * fData->Hf() / fData->G();
}

void TPZMatfrac1dhdiv::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  
#ifdef PZDEBUG
  int nref =  datavec.size();
  if (nref != 2 )
  {
    std::cout << " Erro. The size of the datavec is different from 2 \n";
    DebugStop();
  }
#endif
  
  if(fUpdateMem){
    this->UpdateMemory(datavec);
    return;
  }
  
  // Simulation Data
  const REAL mu = fData->Viscosity();
  const REAL DeltaT = fData->TimeStep();
  
  // Getting shape functions
  TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
  TPZFMatrix<REAL>  &phiP = datavec[1].phi;
  TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
  
  // number of test functions for each state variable
  const int phrQ = phiQ.Rows(), phrP = phiP.Rows();
  
  // blocks
  const int FirstQ  = 0;
  const int FirstP  = phrQ+FirstQ;
  
  // Solutions and derivate of solutions
  const REAL Q = datavec[0].sol[0][0];
  const REAL dQdx = datavec[0].dsol[0](0,0);
  const REAL p = datavec[1].sol[0][0];
  
  // Abertura
  REAL w = fData->GetW(p);
  const REAL dwdp = fData->GetDwDp();
  
  
  // Leak off
  
  const int intGlobPtIndex = datavec[0].intGlobPtIndex;
  TPZFMatrix<REAL> Vl = this->MemItem(intGlobPtIndex);
  const REAL ql = 2.*fData->QlFVl(Vl(0,0),p);
  const REAL dqldp = 2.*fData->dQlFVl(Vl(0,0), p);

  //  Contribution of domain integrals for Jacobian matrix and Residual vector
  //  n + 1 time step
  if(!fData->IsLastState())
  {
    REAL cte = (12.*mu)/(w*w*w);
    for (int iq = 0; iq < phrQ; iq++) { // Constitutive Law
      const REAL res = cte * Q * phiQ(iq,0) - dphiQ(0,iq) * p;
      ef(FirstQ+iq,0) += weight * res;
      for (int jq = 0; jq < phrQ; jq++) {
        const REAL jac = cte * phiQ(iq,0) * phiQ(jq,0);
        ek(FirstQ+iq,FirstQ+jq) += weight * jac;
      }
      for (int jp = 0; jp < phrP; jp++) {
        const REAL jac = - 3 / (w*w*w*w) * dwdp * phiP(jp,0) * 12.*mu * Q * phiQ(iq,0) - dphiQ(0,iq) * phiP(jp,0);
        ek(FirstQ+iq,FirstP+jp) += weight * jac;
      }
    }
    for (int ip = 0; ip < phrP ; ip++) { // Conservation of Mass
      const REAL res = - dQdx * phiP(ip,0) - 1./DeltaT * w * phiP(ip,0) - ql * phiP(ip,0);
      ef(FirstP+ip,0) += weight * res;
      for (int jq = 0; jq < phrQ; jq++) {
        const REAL jac = - dphiQ(0,jq) * phiP(ip,0);
        ek(FirstP+ip,FirstQ+jq) += weight * jac;
      }
      for (int jp = 0; jp < phrP; jp++) {
        const REAL jac = - 1./DeltaT * dwdp * phiP(jp,0) * phiP(ip,0) - dqldp * phiP(jp,0) * phiP(ip,0);
        ek(FirstP+ip,FirstP+jp) += weight * jac;
      }
    }
  }
  
  //  n time step
  //  This values are constant in Newton iteration
  if(fData->IsLastState())
  {
    for (int ip = 0; ip < phrP; ip++) {
      const REAL res = + 1./DeltaT * w * phiP(ip,0);
      ef(FirstP+ip,0) += weight * res;
    }
  }
}

void TPZMatfrac1dhdiv::UpdateMemory(TPZVec<TPZMaterialData> &datavec)
{
  const int intGlobPtIndex = datavec[0].intGlobPtIndex;
  TPZFMatrix<REAL> Vl = this->MemItem(intGlobPtIndex);
  const STATE pfrac = datavec[1].sol[0][0];
  const REAL deltaT = fData->TimeStep();
  
  REAL tStar = fData->FictitiousTime(Vl(0,0), pfrac);
  REAL Vlnext = fData->VlFtau(pfrac, tStar + deltaT);
  
  Vl(0,0) = Vlnext;
  this->MemItem(intGlobPtIndex) = Vl;
}

void TPZMatfrac1dhdiv::UpdateMemory(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec)
{
  const int intGlobPtIndex = data.intGlobPtIndex;
  TPZFMatrix<REAL> Vl = this->MemItem(intGlobPtIndex);
  const STATE pfrac = datavec[1].sol[0][0];
  const REAL deltaT = fData->TimeStep();
  
  REAL tStar = fData->FictitiousTime(Vl(0,0), pfrac);
  REAL Vlnext = fData->VlFtau(pfrac, tStar + deltaT);
  
  Vl(0,0) = Vlnext;
  this->MemItem(intGlobPtIndex) = Vl;
}

void TPZMatfrac1dhdiv::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
  
  if(fUpdateMem){
    this->UpdateMemory(data, dataright);
    return;
  }
  
  // ---------------------- Getting Data for Darcy Flow BC ----------------------
  TPZFMatrix<REAL> &phiQL = dataleft[0].phi;
  TPZFMatrix<REAL> &phiPL = dataleft[1].phi;
  
  TPZManVector<REAL,3> &normal = data.normal;
  
  REAL n1 = normal[0];
  REAL n2 = normal[1];
  
  TPZManVector<STATE,3> sol_qL =dataleft[0].sol[0];
  TPZManVector<STATE,3> sol_pL =dataleft[1].sol[0];
  
  
  //  Getting Q solution for left and right side
  STATE qxL = sol_qL[0];
  STATE qyL = sol_qL[1];
  STATE qnL = (qxL*n1) + (qyL*n2);
  
  //  Getting P solution for left and right side
  STATE PressureL = sol_pL[0];
  
  //  Getting another required data
  REAL TimeStep = fData->TimeStep();
  REAL Theta = fData->Theta();
  
  int QRowsleft = dataleft[0].fVecShapeIndex.NElements();
  int PRowsleft = phiPL.Rows();
  
  int FirstQL = 0;
  int FirstPL = QRowsleft + FirstQL;
  
  // ---------------------- Getting Data for Leak Off ----------------------
  // Getting shape functions
  TPZFMatrix<REAL>  &phiQR = dataright[0].phi;
  TPZFMatrix<REAL>  &phiPR = dataright[1].phi;
//  TPZFMatrix<REAL> &dphiQR = dataright[0].dphix;
  
  // number of test functions for each state variable
  const int phrQR = phiQR.Rows(), phrPR = phiPR.Rows();
  
  // blocks
  const int FirstQR  = 0;
  const int FirstPR  = phrQR+FirstQR;
  
  // Solutions and derivate of solutions
  const REAL p = dataright[1].sol[0][0];
  
  // Leak off
  const int intGlobPtIndex = data.intGlobPtIndex;
  TPZFMatrix<REAL> Vl = this->MemItem(intGlobPtIndex);
  REAL ql, dqldp, dqldpomar;
  if (fData->IsCoupled()) {
    ql = 2.*fData->QlFVl(Vl(0,0), p, PressureL);
    dqldp = 2.*fData->dQlFVl(Vl(0,0), p, PressureL);
    dqldpomar = 2.*fData->dQlFVlPoros(Vl(0,0), p, PressureL);
  }
  else{
    ql = 2.*fData->QlFVl(Vl(0,0), p);
    dqldp = 2.*fData->dQlFVl(Vl(0,0), p);
  }
  
  // ************************ Leak Off Contribution for Integral *************************************
  
  int iblock = QRowsleft + PRowsleft;
  int jblock = QRowsleft + PRowsleft;
  
  //  n + 1 time step
  if(!fData->IsLastState())
  {
    for (int ip = 0; ip < phrPR ; ip++) { // Conservation of Mass
      const REAL res = - ql * phiPR(ip,0);
      ef(FirstPR+ip+ iblock,0) += weight * res;
      for (int jp = 0; jp < phrPR; jp++) {
        const REAL jac = - dqldp * phiPR(jp,0) * phiPR(ip,0);
        ek(FirstPR+ip+iblock,FirstPR+jp+jblock) += weight * jac;
      }
      if (fData->IsCoupled()) {
        for (int jpomar = 0; jpomar < PRowsleft; jpomar++) {
          const REAL jac = - dqldpomar * phiPR(ip,0) * phiPL(jpomar,0);
          ek(FirstPR+ip+iblock, jpomar + FirstPL) += weight * jac;
        }
      }
    }
  }
  
  //  ////////////////////////// Residual Vector ///////////////////////////////////
  //  Contribution of contour integrals for Residual Vector
  //  Time step n
  
  if(fData->IsLastState())
  {
    
    //  Second Block (Equation Two) Bulk flux  equation
    // Integrate[L dot(v, n), Gamma_{e}]    (Equation Two) Left-Left Part
    for (int ip=0; ip < PRowsleft; ip++)
    {
      ef(ip + FirstPL) += (-0.0) * (1.0-Theta) * (TimeStep) * weight * phiPL(ip,0) * qnL;
    }
    
    return;
  }
  
  //  ////////////////////////// Residual Vector ///////////////////////////////////
  //  End of contribution of contour integrals for Residual Vector
  
  
  
  
  //  Second Block (Equation Two) Bulk flux  equation
  // Integrate[L dot(v, n), Gamma_{e}]    (Equation Two) Left-Left Part
  for (int ip=0; ip < PRowsleft; ip++)
  {
    
    ef(ip + FirstPL) += (-1.0) * (Theta) * (TimeStep) * weight * qnL * phiPL(ip,0);
    
    for (int jq=0; jq<QRowsleft; jq++)
    {
      int jvectorindex    = dataleft[0].fVecShapeIndex[jq].first;
      int jshapeindex     = dataleft[0].fVecShapeIndex[jq].second;
      
      REAL vnL =
      (n1) * (phiQL(jshapeindex,0)*dataleft[0].fDeformedDirections(0,jvectorindex)) +
      (n2) * (phiQL(jshapeindex,0)*dataleft[0].fDeformedDirections(1,jvectorindex)) ;
      
      ek(ip + FirstPL,jq + FirstQL) += (-1.0) * (Theta) * (TimeStep) * weight * phiPL(ip,0) * vnL;
      
    }
  }
  
  REAL qN = -ql;  // HERE -> Normal Flux
  
  
  for(int iq=0; iq < QRowsleft; iq++)
  {
    int iLvectorindex       = dataleft[0].fVecShapeIndex[iq].first;
    int iLshapeindex        = dataleft[0].fVecShapeIndex[iq].second;
    
    REAL vni    =   (phiQL(iLshapeindex,0)*dataleft[0].fDeformedDirections(0,iLvectorindex)*n1)+(phiQL(iLshapeindex,0)*dataleft[0].fDeformedDirections(1,iLvectorindex)*n2);
    ef(iq + FirstQL) += weight * ( (gBigNumber * ( qnL - qN ) * vni ) );
    
    for (int jq=0; jq < QRowsleft; jq++)
    {
      int jLvectorindex       = dataleft[0].fVecShapeIndex[jq].first;
      int jLshapeindex        = dataleft[0].fVecShapeIndex[jq].second;
      
      REAL vnj    =   (phiQL(jLshapeindex,0)*dataleft[0].fDeformedDirections(0,jLvectorindex)*n1)+(phiQL(jLshapeindex,0)*dataleft[0].fDeformedDirections(1,jLvectorindex)*n2);
      ek(iq + FirstQL,jq + FirstQL) += weight * ( (gBigNumber * ( vnj ) * vni ) );
    }
  }
  
  
  
  //ApplyQnD(data, dataleft, weight, ek,ef);
}

void TPZMatfrac1dhdiv::ApplyQnD(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef)
{
}


void TPZMatfrac1dhdiv::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  
  if(fData->IsLastState() || this->fUpdateMem){
    return;
  }
  
#ifdef PZDEBUG
  int nref =  datavec.size();
  if (nref != 2 )
  {
    std::cout << " Erro. The size of the datavec is different from 2 \n";
    DebugStop();
  }
#endif
  
  //BigNumber
  const REAL BIGNUMBER = TPZMaterial::gBigNumber;
  
  // Getting shape functions
  TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
  
  // number of test functions for each state variable
  const int phrQ = phiQ.Rows();
  
  // blocks
  const int FirstQ  = 0;
  
  // Solutions and derivate of solutions
  const REAL Q = datavec[0].sol[0][0];
  
  // Valores impostos
  REAL v2 = bc.Val2()(0,0);
  
  switch (bc.Type()) {
    case 0: // Flux strong imposition (v2 = imposed flux)
    {
      const REAL difQ = Q - v2;
      for(int iq = 0 ; iq < phrQ ; iq++) {
        ef(FirstQ+iq,0) += BIGNUMBER * difQ * phiQ(iq,0) * weight; // forced v2 flux
        for (int jq = 0 ; jq < phrQ ; jq++)
        {
          ek(FirstQ+iq,FirstQ+jq) += BIGNUMBER * phiQ(iq,0) * phiQ(jq,0) * weight;
        }
      }
    }
      break;
    case 1: // Pressure weak imposition (v2 = imposed pressure)
    {
      for (int iq = 0; iq < phrQ; iq++) {
        ef(FirstQ+iq,0) += v2 * phiQ(iq,0);
      }
    }
      break;
      
  }
  
}

/** Returns the variable index associated with the name */
int TPZMatfrac1dhdiv::VariableIndex(const std::string &name){
  if(!strcmp("Pressure",name.c_str())) return  0;
  if(!strcmp("State",name.c_str())) return  0;
  if(!strcmp("Flow",name.c_str())) return  1;
  if(!strcmp("Flux",name.c_str())) return  1;
  if(!strcmp("Opening",name.c_str())) return  2;
  
  return TPZMaterial::VariableIndex(name);
}

int TPZMatfrac1dhdiv::NSolutionVariables(int var){
  if(var == 0) return 1;
  if(var == 1) return 1;
  if(var == 2) return 1;
  
  return TPZMaterial::NSolutionVariables(var);
}

void TPZMatfrac1dhdiv::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
  
  Solout.Resize(this->NSolutionVariables(var));
  const REAL SolQ = datavec[0].sol[0][0];
  const REAL SolP = datavec[1].sol[0][0];
  
  if(var == 0){ //function (state variable Q)
    Solout[0] = SolP;
    return;
  }
  
  if(var == 1){
    Solout[0] = SolQ;//function (state variable p)
    return;
  }
  
  if (var == 2) {
    const REAL w = fData->GetW(SolP);
    Solout[0] = w;
  }
}


void TPZMatfrac1dhdiv::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
  int nref = datavec.size();
  for(int i = 0; i<nref; i++ )
  {
    datavec[i].SetAllRequirements(false);
    datavec[i].fNeedsSol = true;
  }
  
}

void TPZMatfrac1dhdiv::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){
  int nref = datavec.size();
  for(int i = 0; i<nref; i++)
  {
    datavec[i].SetAllRequirements(false);
    datavec[i].fNeedsSol = true;
  }
}
