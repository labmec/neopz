#include "pzlog.h"
#include "TPZMatfrac1dhdiv.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"

#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphase"));
#endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.multiphase.data"));
#endif

TPZMatfrac1dhdiv::TPZMatfrac1dhdiv(): TPZMaterial()
{
  fDim = 1;
}

TPZMatfrac1dhdiv::TPZMatfrac1dhdiv(int matid): TPZMaterial(matid)
{
  fDim = 1;
}


TPZMatfrac1dhdiv::~TPZMatfrac1dhdiv()
{
  
}

int TPZMatfrac1dhdiv::Dimension() const {return fDim;};

int TPZMatfrac1dhdiv::NStateVariables() {return 1;}

void TPZMatfrac1dhdiv::Print(std::ostream &out) {
  out << "name of material : " << Name() << "\n";
  out << "Coeficient which multiplies the gradient operator "<< "my var" << std::endl;
  out << "Base Class properties :";
  TPZMaterial::Print(out);
  out << "\n";
}


void TPZMatfrac1dhdiv::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  
#ifdef DEBUG
  int nref =  datavec.size();
  if (nref != 2 )
  {
    std::cout << " Erro. The size of the datavec is different from 2 \n";
    DebugStop();
  }
#endif
  
  // Getting shape functions
  TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
  TPZFMatrix<REAL>  &phiP = datavec[1].phi;
  TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
  TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
  
  // number of test functions for each state variable
  const int phrQ = phiQ.Rows(), phrP = phiP.Rows();
  
  // blocks
  const int FirstQ  = 0;
  const int FirstP  = phrQ+FirstQ;
  
  // Solutions and derivate of solutions
  const REAL Q = datavec[0].sol[0][0];
  const REAL dQ = datavec[0].dsol[0](0,0);
  const REAL p = datavec[1].sol[0][0];
  
  // Abertura
  const REAL w = 1.; // sera uma funcao de p
  const REAL dwdp = 0.;
  
  // Leak Off
  const REAL ql = datavec[0].x[0];
  const REAL dqldp = 0.;
  
  // Simulation Data
  const REAL mu = fData->Viscosity();
  const REAL DeltaT = fData->TimeStep();
  
  //  Contribution of domain integrals for Jacobian matrix and Residual vector
  //  n + 1 time step
  if(!fData->IsLastState())
  {
    REAL cte = w*w*w/(12.*mu);
    for (int iq = 0; iq < phrQ; iq++) {
      const REAL res = Q * phiQ(iq,0) - cte * dphiQ(0,iq) * p;
      ef(FirstQ+iq,0) += weight * res;
      for (int jq = 0; jq < phrQ; jq++) {
        const REAL jac = phiQ(iq,0) * phiQ(jq,0);
        ek(FirstQ+iq,FirstQ+jq) += weight * jac;
      }
      for (int jp = 0; jp < phrP; jp++) {
        const REAL jac = - 3 *  w * w / (12. * mu) * dwdp * phiP(jp,0) * dphiQ(0,iq) * p - cte * dphiQ(0,iq) * phiP(jp,0);
        ek(FirstQ+iq,FirstP+jp) += weight * jac;
      }
    }
    for (int ip = 0; ip < phrP ; ip++) {
      const REAL res = dQ * phiP(ip,0) + 1./DeltaT * w * phiP(ip,0) + ql * phiP(ip,0);
      ef(FirstP+ip,0) += weight * res;
      for (int jq = 0; jq < phrQ; jq++) {
        const REAL jac = dphiQ(0,jq) * phiP(ip,0);
        ek(FirstP+ip,FirstQ+jq) += weight * jac;
      }
      for (int jp = 0; jp < phrP; jp++) {
        const REAL jac = 1./DeltaT * dwdp * phiP(jp,0) * phiP(ip,0) + dqldp * phiP(jp,0) * phiP(ip,0);
        ek(FirstP+ip,FirstP+jp) += weight * jac;
      }
    }
  }
  
  //  n time step
  //  This values are constant in Newton iteration
  if(fData->IsLastState())
  {
    for (int ip = 0; ip < phrP; ip++) {
      const REAL res = - 1./DeltaT * w * phiP(ip,0);
      ef(FirstP+ip,0) += weight * res;
    }
  }
  
}

void TPZMatfrac1dhdiv::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{

  if(fData->IsLastState()){
    return;
  }
  
#ifdef DEBUG
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
  //TPZFMatrix<REAL>  &phiP = datavec[1].phi;
  
  // number of test functions for each state variable
  const int phrQ = phiQ.Rows();
  
  // blocks
  const int FirstQ  = 0;

  // Solutions and derivate of solutions
  const REAL Q = datavec[0].sol[0][0];
  //const REAL p = datavec[1].sol[0][0];
  
  // Abertura
  const REAL w = 1.; // sera uma funcao de p
  //const REAL dwdp = 0.;
  
  // Dados da simulacao
  const REAL mu = fData->Viscosity();
  
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
        ef(FirstQ+iq,0) += w*w*w/(12.*mu) * v2 * phiQ(iq,0);
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
    const REAL w = 1.; // AQUINATHAN Trocar pela formula do w
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
