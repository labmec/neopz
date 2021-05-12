

#include "TPZLinearConvecDiff.h"
#include "pzaxestools.h"
#include "pzbndcond.h"

TPZLinearConvecDiff::TPZLinearConvecDiff(int nummat, REAL k, const TPZVec<REAL> &conv, REAL f, REAL SD)
:TPZRegisterClassId(&TPZLinearConvecDiff::ClassId), TPZMaterial(nummat){
  fK = k;
  fConvDir[0] = conv[0];
  fConvDir[1] = conv[1];
  fXf = f;
  fSD = SD;
}

TPZLinearConvecDiff::TPZLinearConvecDiff(int matid)
:TPZRegisterClassId(&TPZLinearConvecDiff::ClassId),  TPZMaterial(matid), fXf(0.), fK(0.), fSD(0.){
  fConvDir[0] = 0.;
  fConvDir[1] = 0.;
}

TPZLinearConvecDiff::TPZLinearConvecDiff():
TPZRegisterClassId(&TPZLinearConvecDiff::ClassId), TPZMaterial(), fXf(0.), fK(0.), fSD(0.){
  fConvDir[0] = 0.;
  fConvDir[1] = 0.;
}

TPZLinearConvecDiff::TPZLinearConvecDiff(const TPZLinearConvecDiff &c)
:TPZRegisterClassId(&TPZLinearConvecDiff::ClassId), TPZMaterial(c){
  fK = c.fK;
  fConvDir[0] = c.fConvDir[0];
  fConvDir[1] = c.fConvDir[1];
  fXf = c.fXf;
  fSD = c.fSD;
}

TPZLinearConvecDiff::~TPZLinearConvecDiff(){
///nothing here
}

void TPZLinearConvecDiff::Print(std::ostream & out){
	out << "name of material : " << Name() << "\n";
	out << "fK "<< fK << std::endl;
	out << "Convection vector " << fConvDir[0] << "\t" << fConvDir[1] << std::endl;
	out << "fXf " << fXf << std::endl;
  out << "fSD " << fSD << std::endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}

void TPZLinearConvecDiff::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){

  TPZFMatrix<REAL>  &phi = data.phi;
  TPZFMatrix<REAL> &dphidaxes = data.dphix;
  TPZFNMatrix<200,REAL> dphi(dphidaxes.Rows(),dphidaxes.Cols(),0.);
  TPZAxesTools<REAL>::Axes2XYZ(dphidaxes, dphi, data.axes);

  TPZVec<REAL>  &x = data.x;
//  TPZFMatrix<REAL> &axes = data.axes;
//  TPZFMatrix<REAL> &jacinv = data.jacinv;
  const int nshape = phi.Rows();

  STATE FVal = this->fXf;
  if(fForcingFunction) {
      TPZManVector<STATE,1> res(1);
      TPZFMatrix<STATE> dres(Dimension(),1);
      fForcingFunction->Execute(x,res,dres);
      FVal = res[0];
  }
  const REAL normaConveccao = sqrt(fConvDir[0]*fConvDir[0]+fConvDir[1]*fConvDir[1]);

  const REAL h = data.HSize;
  for( int in = 0; in < nshape; in++ ) {
    const REAL gradVBeta = this->fConvDir[0] * dphi(0,in) + this->fConvDir[1] * dphi(1,in);
    ef(in, 0) += weight * FVal * ( phi(in,0) + this->fSD*(0.5*h/normaConveccao)*gradVBeta );
    for( int jn = 0; jn < nshape; jn++ ) {
      ek(in,jn) += weight * (
        +fK * ( dphi(0,in) * dphi(0,jn) + dphi(1,in) * dphi(1,jn) )
        - ( (dphi(0,in) * phi(jn)) * fConvDir[0] + (dphi(1,in) * phi(jn)) * fConvDir[1] )
        + this->fSD*(0.5*h/normaConveccao)*( 
                                             (fConvDir[0]*dphi(0,jn))*(fConvDir[0]*dphi(0,in)) +
                                             (fConvDir[1]*dphi(1,jn))*(fConvDir[1]*dphi(1,in))   )
        );
    }
  }
}///void

void TPZLinearConvecDiff::ContributeBC(TPZMaterialData &data,REAL weight,
							              TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
	TPZFMatrix<REAL>  &phi = data.phi;
	TPZFMatrix<REAL> &axes = data.axes;
	const int nshape = phi.Rows();
	STATE v2 = bc.Val2()(0,0);

	if(bc.HasForcingFunction()) {
		TPZManVector<STATE> res(1);
		bc.ForcingFunction()->Execute(data.x,res);
		v2 = res[0];
	}

	switch (bc.Type()) {
		case 0 :			// Dirichlet condition
			for(int in = 0 ; in < nshape; in++) {
				ef(in,0) += weight * gBigNumber* phi(in,0) * v2;
				for (int jn = 0 ; jn < nshape; jn++) {
					ek(in,jn) += weight * gBigNumber * phi(in,0) * phi(jn,0);
				}
			}
			break;
		case 1 :			// Neumann condition
			for(int in = 0 ; in < nshape; in++) {
				ef(in,0) += weight * v2 * phi(in,0);
			}
			break;
		case 3: // outflow condition
			REAL normal[2];
			normal[0] = axes(0,1);
			normal[1] = axes(1,1);

			REAL ConvNormal = 0.;
			for(int id = 0; id < 2; id++) ConvNormal += fConvDir[id]*normal[id];
			if(ConvNormal > 0.){
				for(int il = 0; il < nshape; il++){
					for(int jl = 0; jl < nshape; jl++){
						ek(il,jl) += weight * ConvNormal * phi(il)*phi(jl);
					}
				}
			}
			else{
				if (ConvNormal < 0.) std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
			}
			break;
	}

}

int TPZLinearConvecDiff::VariableIndex(const std::string &name){
	if(!strcmp("Solution",name.c_str()))        return  1;
	if(!strcmp("Derivative",name.c_str()))      return  2;
	return TPZMaterial::VariableIndex(name);
}

int TPZLinearConvecDiff::NSolutionVariables(int var){
	if(var == 1) return 1;
	if(var == 2) return 2;
	return TPZMaterial::NSolutionVariables(var);
}

void TPZLinearConvecDiff::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){

	Solout.Resize( this->NSolutionVariables( var ) );

	if(var == 1){
		Solout[0] = data.sol[0][0];//solution - escalar
		return;
	}
	if(var == 2) {
    TPZFNMatrix<9,STATE> dsoldx;
    TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
		for(int id = 0 ; id < 2; id++) {
			Solout[id] = dsoldx(id,0);//derivative - vetorial
		}
		return;
	}//var == 2

  return TPZMaterial::Solution(data,var,Solout);
}

void TPZLinearConvecDiff::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values){

	values.Resize(3);
	///L2 norm
	values[1] = (u[0] - u_exact[0])*(u[0] - u_exact[0]);
	///semi norma de H1
	values[2] = 0.;
	for(int i = 0; i < 2; i++){
		values[2] += (dudx(i,0) - du_exact(i,0))*(dudx(i,0) - du_exact(i,0));
	}
	///H1 norm
	values[0] = values[1]+values[2];

}

int TPZLinearConvecDiff::ClassId() const{
    return Hash("TPZLinearConvecDiff") ^ TPZMaterial::ClassId() << 1;
}
