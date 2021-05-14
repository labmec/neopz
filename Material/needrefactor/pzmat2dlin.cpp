/**
 * \file
 * @brief Contains implementations of the TPZMat2dLin methods.
 */

#include "pzmat2dlin.h"
#include "TPZMaterial.h"
#include "pzconnect.h"
#include "pzbndcond.h"
#include <math.h>
#include "pzvec.h"
#include "pzerror.h"
#include "TPZStream.h"

#include "pzpoisson3d.h"

using namespace std;

void TPZMat2dLin::Contribute( TPZMaterialData &data,REAL weight,
							 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef )
{
	// this method adds the contribution of the material to the stiffness
	// matrix and right hand side
	
	// check on the validity of the arguments
	//rows x cols
	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
	TPZManVector<REAL,3> &x = data.x;
	TPZFMatrix<REAL> &axes=data.axes;
	
	if(phi.Cols() != 1 || dphi.Rows() != 2 || phi.Rows() != dphi.Cols())
	{
		PZError << "TPZMat2dLin.contr, inconsistent input data : phi.Cols() = "
		<< phi.Cols() << " dphi.Cols + " << dphi.Cols()
		<< " phi.Rows = " << phi.Rows() << " dphi.Rows = " << dphi.Rows() << "\n";
        DebugStop();
	}
	if(fForcingFunction)
	{
		TPZManVector<STATE> xfloat(fXf.Rows());
		fForcingFunction->Execute(x,xfloat);//fXf = xfloat
		int i;
		for(i=0; i<fXf.Rows(); i++) fXf(i,0) = xfloat[i];
	}
	int r = fKxx.Rows();
	int ic,jc;
	for(int in=0 ; in < phi.Rows() ; ++in)
	{
		for(ic=0; ic< r; ic++)
		{
			ef(in*r+ic, 0) +=  (STATE)(weight*phi(in,0))*fXf(ic,0);
		}
		for(int jn=0 ; jn<phi.Rows() ; ++jn)
		{
			for(ic=0; ic<r; ic++) for(jc=0; jc<r; jc++)
			{
				REAL dphix = dphi(0,in)*axes(0,0)+axes(1,0)*dphi(1,in);
				REAL dphiy = dphi(0,in)*axes(0,1)+axes(1,1)*dphi(1,in);
				REAL dphjx = dphi(0,jn)*axes(0,0)+axes(1,0)*dphi(1,jn);
				REAL dphjy = dphi(0,jn)*axes(0,1)+axes(1,1)*dphi(1,jn);
				ek(in*r+ic,jn*r+jc) += (STATE)(weight*phi(in,0)*phi(jn,0))*fK00(ic,jc) +
				((STATE)(dphix*dphjx)*fKxx(ic,jc) +
				 (STATE)(dphiy*dphjy)*fKyy(ic,jc) +
				 (STATE)(dphix*phi(jn,0))*fKx0(ic,jc) +
				 (STATE)(dphiy*phi(jn,0))*fKy0(ic,jc) +
				 (STATE)(phi(in,0)*dphjx)*fK0x(ic,jc) +
				 (STATE)(phi(in,0)*dphjy)*fK0y(ic,jc) ) * (STATE)weight;
			}
		}
	}
}

void TPZMat2dLin::ContributeBC(TPZMaterialData &data,REAL weight,
							   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
	TPZFMatrix<REAL> &phi = data.phi;
	TPZFMatrix<REAL> &axes=data.axes;
	
	if(bc.Material() != this){
		PZError << "TPZMat1dLin.apply_bc warning : this material didn't create the boundary condition!\n";
	}
	
	if(bc.Type() < 0 && bc.Type() > 2){
		PZError << "TPZMat1dLin.aplybc, unknown boundary condition type :"  <<
		bc.Type() << " boundary condition ignored\n";
	}
	
	int r = fKxx.Rows();
	int numnod = ek.Rows()/r;
	
	int idf,jdf,in,jn;
	switch(bc.Type()){
			
		case 0:
			for(in=0 ; in<numnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					(ef)(in*r+idf,0) += (STATE)(weight*gBigNumber*phi(in,0))*bc.Val2()(idf,0);
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
					(ef)(in*r+idf,0) += (STATE)(weight*phi(in,0))*bc.Val2()(idf,0);
				}
			}
			break;
			
		case 2:
			for(in=0 ; in<numnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					(ef)(in*r+idf,0) += (STATE)(weight*phi(in,0))*bc.Val2()(idf,0);
				}
				for(jn=0 ; jn<numnod ; ++jn) {
					for(idf = 0;idf<r;idf++) {
						for(jdf = 0;jdf<r;jdf++) {
							ek(in*r+idf,jn*r+jdf) += bc.Val1()(idf,jdf)*(STATE)(phi(in,0)*phi(jn,0)*weight);
						}
					}
				}
			}
			break;
		case 3:
			for(in=0; in<numnod; in++) {
				for(idf=0; idf<r; idf++) {
					for(jdf=0; jdf<r; jdf++) {
						ef(in*r+idf,0) += (fKx0(idf,jdf)*(STATE)axes(0,1)-fKy0(idf,jdf)*(STATE)axes(0,0))*(STATE)(phi(in,0)*weight);
					}
				}
			}
			break;
	}//fim switch
}

void TPZMat2dLin::Print(std::ostream & out) {
	out << "Material type TPZMat2dLin -- number = " << Id() << "\n";
	out << "Matrix Kxx ->  "; fKxx.Print("fKxx",out);
	out << "Matrix Kyy ->  "; fKyy.Print("fKyy",out);
	out << "Matrix Kxy ->  "; fKxy.Print("fKxy",out);
	out << "Matrix Kyx ->  "; fKyx.Print("fKyx",out);
	out << "Matrix Kx0 ->  "; fKx0.Print("fKx0",out);
	out << "Matrix K0x ->  "; fK0x.Print("fK0x",out);
	out << "Matrix Ky0 ->  "; fKy0.Print("fKy0",out);
	out << "Matrix K0y ->  "; fK0y.Print("fK0y",out);
	out << "Matrix K00 ->  "; fK00.Print("fK00",out);
	out << "Matrix xf ->  "; fXf.Print("fXf",out);
}

int TPZMat2dLin::VariableIndex(const std::string &name) {
	if(!strcmp(name.c_str(),"displacement")) return 1;
	if(!strcmp(name.c_str(),"Solution")) return 1;
	if(!strcmp("Derivate",name.c_str())) return 2;
    if(!strcmp("Pressure",name.c_str())) return 3;
	return TPZMaterial::VariableIndex(name);
}

int TPZMat2dLin::NSolutionVariables(int index) {
	if(index == 1) return 3;
	if(index == 2) return 2;
    if(index == 3) return 1;
	return TPZMaterial::NSolutionVariables(index);
}


void TPZMat2dLin::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes, int var,TPZVec<STATE> &Solout) {
	
	if(var == 0) {
#ifndef STATE_COMPLEX
        Solout.Resize(Sol.size());
        for (int i = 0; i<Sol.size(); i++) {
            Solout[i] = Sol[i];
        }
#else
        for (int i=0; i<Sol.size(); i++) {
            Solout[i] = Sol[i].real();
        }
#endif
	} else if(var == 1) {
		Solout.Resize(3,0.);
#ifndef STATE_COMPLEX
		Solout[0] = Sol[0];
#else
        Solout[0] = Sol[0].real();
#endif
	} else if(var == 2) {
#ifndef STATE_COMPLEX
      	Solout[0] = DSol(0,0);//derivate
		Solout[1] = DSol(1,0);//derivate
#else
      	Solout[0] = DSol(0,0).real();//derivate
		Solout[1] = DSol(1,0).real();//derivate
#endif
		return;
	} else if(var == 3) {
		Solout.Resize(1,0.);
#ifndef STATE_COMPLEX
		Solout[0] = Sol[0];
#else
        Solout[0] = Sol[0].real();
#endif
	}
    
	else TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

void TPZMat2dLin::Errors(TPZVec<REAL> &/*x*/,TPZVec<STATE> &u,TPZFMatrix<STATE> &dudx,TPZFMatrix<REAL> &axes,
						 TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
	
	STATE dx = dudx(0,0)*(STATE)axes(0,0)+dudx(1,0)*(STATE)axes(1,0);
	STATE dy = dudx(0,0)*(STATE)axes(0,1)+dudx(1,0)*(STATE)axes(1,1);
	//Norma Energia
	STATE parc1 = dx-du_exact(0,0) ;
	STATE parc2 = dy-du_exact(1,0) ;
	values[0] = abs(parc1*parc1 + parc2*parc2);/*pow(parc1,2.)+pow(parc2,2.);*/
	//Norma L2
	values[1] = pow((REAL)fabs(u[0] - u_exact[0]),(REAL)2.0);
	values[2] = 0.;
}

TPZMaterial * TPZMat2dLin::NewMaterial() {
	return new TPZMat2dLin(*this);
}

void TPZMat2dLin::ConvectionDiffusion(REAL angle,REAL diff){
	fKxx.Redim(1,1);
	fKxy.Redim(1,1);
	fKyx.Redim(1,1);
	fKyy.Redim(1,1);
	fKx0.Redim(1,1);
	fK0x.Redim(1,1);
	fKy0.Redim(1,1);
	fK0y.Redim(1,1);
	fK00.Redim(1,1);
	fXf.Redim(1,1);
	REAL cosa = cos(angle);
	REAL sina = sin(angle);
	fKxx(0,0) = cosa*cosa;
	fKxy(0,0) = cosa*sina;
	fKyy(0,0) = sina*sina;
	fKyx(0,0) = cosa*sina;
	fKx0(0,0) = -cosa;
	fKy0(0,0) = -sina;
}

TPZBndCond *TPZMat2dLin::OutflowFlux(TPZMaterial * &reference, int bc){
	
	int nstate = fKxx.Rows();
	TPZFMatrix<STATE> val1(nstate,nstate),val2(nstate,1);
	return TPZMaterial::CreateBC(reference,bc,3,val1,val2);
}

/** TPZSavable methods ***/

/** returns the unique identifier for reading/writing objects to streams */
int TPZMat2dLin::ClassId() const{
    return Hash("TPZMat2dLin") ^ TPZMaterial::ClassId() << 1;
}

/** Save the element data to a stream */
void TPZMat2dLin::Write(TPZStream &buf, int withclassid) const
{
#ifdef PZDEBUG2
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " before write material ";
		LOGPZ_DEBUG( logger,sout.str().c_str() );
	}
#endif
	TPZMaterial::Write(buf,withclassid);
#ifdef PZDEBUG2
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " after write material ";
		LOGPZ_DEBUG( logger,sout.str().c_str() );
	}
#endif
	
	fKxx.Write(buf,0);
	fKxy.Write(buf,0);
	fKyx.Write(buf,0);
	fKyy.Write(buf,0);
	fKx0.Write(buf,0);
	fK0x.Write(buf,0);
	fKy0.Write(buf,0);
	fK0y.Write(buf,0);
	fK00.Write(buf,0);
	fXf.Write(buf,0);
}

/** Read the element data from a stream */
void TPZMat2dLin::Read(TPZStream &buf, void *context)
{
	TPZMaterial::Read(buf,context);
	fKxx.Read(buf,0);
	fKxy.Read(buf,0);
	fKyx.Read(buf,0);
	fKyy.Read(buf,0);
	fKx0.Read(buf,0);
	fK0x.Read(buf,0);
	fKy0.Read(buf,0);
	fK0y.Read(buf,0);
	fK00.Read(buf,0);
	fXf.Read(buf,0);
}

#ifndef BORLAND
template class TPZRestoreClass<TPZMat2dLin>;
#endif
