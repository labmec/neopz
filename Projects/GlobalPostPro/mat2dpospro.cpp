/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2014  <copyright holder> <email>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "pzmaterialdata.h"
#include <math.h>
#include "mat2dpospro.h"
#include "pzmatrix.h"
#include "pzaxestools.h"

Mat2Dpospro::Mat2Dpospro() : TPZMaterial(), fDim(0){ //void constructor

  fAlpha = 1.;
  fDelta = 1.; 
}

Mat2Dpospro::Mat2Dpospro(int nummat, int dim) : TPZMaterial(nummat), fDim(dim){ //constructor

  fAlpha = 1.;
  fDelta = 1.; 
}

Mat2Dpospro::~Mat2Dpospro()
{
  
}

Mat2Dpospro::Mat2Dpospro(const Mat2Dpospro& other) : TPZMaterial(other)
{
  fAlpha = other.fAlpha;
  fDelta = other.fDelta;
}

Mat2Dpospro & Mat2Dpospro::operator=(const Mat2Dpospro &copy)
{
    TPZMaterial::operator=(copy);
    fAlpha = copy.fAlpha;
    fDelta = copy.fDelta;
    return *this;
}


int Mat2Dpospro::Dimension() const {return fDim;};


void Mat2Dpospro::Print(std::ostream &out) {
  out << "name of material : " << Name() << "\n";
  out << "Coeficient which multiplies the gradient operator "<< "my var" << std::endl;
  out << "Base Class properties :";
  TPZMaterial::Print(out);
  out << "\n";
}

void Mat2Dpospro::SetParameters(REAL alpha, REAL delta) {
  fAlpha = alpha;
  fDelta = delta;
}

void Mat2Dpospro::GetParameters(REAL &alpha, REAL &delta) {
    alpha = fAlpha;
    delta = fDelta;
}


void Mat2Dpospro::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	//Find u in V c Hdiv, tal a_alpha(u,v) = L_alpha(v) para todo v in V
	//where a_alpha(u,v) = integral(u dot v) +  co^alpha integral(div u div v)
	//e L_alpha(v) = 
	//(delta h)^alpha = co^alpha
    REAL h=data.HSize;
    REAL co = fDelta*h;
    REAL coef = pow(co,fAlpha);
    
    TPZFMatrix<REAL> &phiTau = data.phi;
    TPZFMatrix<REAL> &dphiTau = data.dphix;
    TPZVec<REAL> divphi(fDim);//zerar
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = data.fVecShapeIndex.NElements();
    
    REAL F = 0, U = 0;
    TPZFMatrix<STATE> GradU(fDim,1,0.0);
    
    
    if(fForcingFunction) {
        TPZVec<STATE> solExact(1);
        fForcingFunction->Execute(data.x, solExact);
        F = solExact[0];
    }
    if(fTimeDependentForcingFunction)
    {
        TPZVec<STATE> solExact(1);
	REAL time =0.0;
        fTimeDependentForcingFunction->Execute(data.x,time, solExact,GradU);
	U = solExact[0];
    }
    
    
    // Matrix A
    for(int iq=0; iq<phr; iq++)
    {
	int ivecind = data.fVecShapeIndex[iq].first;
	int ishapeind = data.fVecShapeIndex[iq].second;
	
	TPZFMatrix<REAL> Vi(3,1);
	Vi(0,0) = data.fDeformedDirections(0,ivecind);
	Vi(1,0) = data.fDeformedDirections(1,ivecind);
	Vi(2,0) = data.fDeformedDirections(2,ivecind);
	TPZFMatrix<REAL> axesveci(3,1);
	data.axes.Multiply(Vi,axesveci);
	
	REAL divphitau_i = 0.;
	for(int iloc=0; iloc<fDim; iloc++)
	{
	    divphitau_i += dphiTau(iloc,ishapeind)*axesveci(iloc,0);
	}	
	
	REAL GradudotphiTau=0.0;
	REAL fdivphitau_i = F * divphitau_i;
	
	GradudotphiTau = GradU[0]*phiTau(ishapeind,0)*Vi(0,0) + GradU[1]*phiTau(ishapeind,0)*Vi(1,0);
	
	ef(iq) += weight*(GradudotphiTau - coef*fdivphitau_i);
	
	for (int jq=0; jq<phr; jq++)
	{
	    int jvecind = data.fVecShapeIndex[jq].first;
	    int jshapeind = data.fVecShapeIndex[jq].second;
	    
	    TPZFMatrix<REAL> Vj(3,1);
	    Vj(0,0) = data.fDeformedDirections(0,jvecind);
	    Vj(1,0) = data.fDeformedDirections(1,jvecind);
	    Vj(2,0) = data.fDeformedDirections(2,jvecind);
	    TPZFMatrix<REAL> axesvecj(3,1);
	    data.axes.Multiply(Vj,axesvecj);
	    
	    REAL divphitau_j = 0.;
	    for(int jloc=0; jloc<fDim; jloc++)
	    {
		divphitau_j += dphiTau(jloc,jshapeind)*axesvecj(jloc,0);
	    }		    
	    
	    REAL dotprod = 
	    phiTau(ishapeind,0)*Vi(0,0) * phiTau(jshapeind,0)*Vj(0,0)+
	    phiTau(ishapeind,0)*Vi(1,0) * phiTau(jshapeind,0)*Vj(1,0);
	    
	    //dot product between phi_i and phi_j
	    ek(iq,jq) += weight*(dotprod + coef*(divphitau_i*divphitau_j));
	}
    }    
}

/** Returns the variable index associated with the name */
int Mat2Dpospro::VariableIndex(const std::string &name){
  if(!strcmp("U",name.c_str()))		return  1;
  if(!strcmp("F",name.c_str()))		return  2;
  if(!strcmp("Sigma",name.c_str()))		return  3;
  if(!strcmp("SigmaExact",name.c_str()))	return  4;
  if(!strcmp("GradU",name.c_str()))		return  5;
  
  return TPZMaterial::VariableIndex(name);
}

int Mat2Dpospro::NSolutionVariables(int var){
  if(var == 1) return 1;
  if(var == 2) return 1;
  if(var == 3) return fDim;
  if(var == 4) return fDim;
  if(var == 5) return fDim; 
  
  return TPZMaterial::NSolutionVariables(var);
}

void Mat2Dpospro::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    
    
    
    TPZVec<STATE> SolSigma;
    SolSigma = data.sol[0];
    
    if(var == 1){ //function (U)
        TPZVec<STATE> solExact(1);
        TPZFMatrix<STATE> flux(fDim,1);
	REAL time =0.0;
        fTimeDependentForcingFunction->Execute(data.x,time, solExact,flux);
        Solout[0] = solExact[0];
        return;
    }
    
    if(var == 2){//function (F)
        TPZVec<STATE> solExact(1);
        fForcingFunction->Execute(data.x, solExact);
        Solout[0] = solExact[0];
        return;
    }

    if(var == 3){ //function (state variable Q anal)
        Solout[0] = SolSigma[0];
        Solout[1] = SolSigma[1];
        return;
    }
    
    if(var == 4){//function (state variable p anal)
        TPZVec<STATE> solExact(fDim);
        TPZFMatrix<STATE> flux(fDim,1);
        fExactSol->Execute(data.x, solExact,flux);
        Solout[0] = solExact[0];
        Solout[1] = solExact[1];
        return;
    }
    if(var == 5){//function (state variable p anal)
        TPZVec<STATE> solExact(fDim);
        TPZFMatrix<STATE> flux(fDim,1);
	REAL time =0.0;
        fTimeDependentForcingFunction->Execute(data.x,time, solExact,flux);
        Solout[0] = flux(0,0);
        Solout[1] = flux(1,0);
        return;
    }    
    
}


void Mat2Dpospro::FillDataRequirements(TPZMaterialData  &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = false;
    data.fNeedsHSize = true;
    data.fNeedsNeighborSol = false;
    data.fNeedsNormal = false;
    data.fNeedsNeighborCenter = false;
}

void Mat2Dpospro::FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data){

    data.fNeedsSol = false;
    data.fNeedsNormal = false;
    data.fNeedsNeighborSol = false;
}

